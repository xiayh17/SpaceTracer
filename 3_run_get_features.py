# spatial features 
# phase combine 
# phase no combine
# extract features
import argparse
from functools import partial
import multiprocessing
import os
import shutil
import subprocess
import uuid
from module.extract_features import extract_feature_perline, get_mutation_bed_file, handle_h5ad_file, handle_mappbablity_file
from module.phase_combine import handel_candidate_informative_SNP_site, phase_combine_get_candidate_germline, short_bed_file, short_gencode_file, short_mosaic_sites_and_combine_with_germ
from module.phase_no_combine import phase_no_combine_get_candidate_germline
from module.read_file import handle_barcode, handle_spot_geno_to_get_mutation_list, read_spot_posterior_file
import pandas as pd
from module.spatial_features import handle_per_line, parallel_merge_and_save
import numpy as np
from utils import barcode_cell_mapping, check_dir, check_input, check_output, str2bool
from tqdm import tqdm


def main():
    check_dir(args.outdir)
    check_input([args.ind_genotype,args.spot_genotype,args.fasta,args.germline,args.artifact_signature])
    in_name=args.outprefix
    indgenotype_file=args.ind_genotype
    spotgenotype_file=args.spot_genotype
    germline_file=args.germline
    ref_fasta=args.fasta
    filter_bam=args.filter_bam
    raw_bam=args.raw_bam
    thread=int(args.thread)
    artifact_signature=args.artifact_signature

    barcode_dir=os.path.join(args.outdir,"barcode_dir");check_dir(barcode_dir)
    spatial_feature_file=os.path.join(args.outdir,in_name+".spatial_feature.txt")
    # phase_beforeUMIcombination_file=os.path.join(args.outdir,in_name+".phase_beforeUMIcombination.txt")
    phase_afterUMIcombination_file=os.path.join(args.outdir,in_name+".phase_afterUMIcombination.txt")
    feature_file=os.path.join(args.outdir,in_name+".features.txt")
    add_hFR_feature=os.path.join(args.outdir,"all_RF_filter_hFDR_less_08_feature_remove_AtoG.txt")
    
    def get_spatial_test_result():
        mutation_identifier_list=handle_spot_geno_to_get_mutation_list(indgenotype_file)
        barcode_dict=handle_barcode(args.barcodes)
        # get the spatial scatter figure range
        barcode_df = pd.read_csv(args.barcodes)
        min_array_row, max_array_row = barcode_df['array_row'].agg(['min', 'max'])
        min_array_col, max_array_col = barcode_df['array_col'].agg(['min', 'max'])

        array_size_row = max_array_row - min_array_row  # range of array_row
        array_size_col = max_array_col - min_array_col  # range of array_col
        array_area = array_size_row * array_size_col

        if args.point_size is not None:
            point_size = args.point_size
        else:
            num_points = len(barcode_df)
            base_point_size = 1.5
            point_density = num_points / array_area if array_area != 0 else 1
            point_size = round((base_point_size / point_density)**2, 1)
            min_point_size = 5
            max_point_size = 50
            point_size = np.clip(point_size, min_point_size, max_point_size)

        spot_geno_file=read_spot_posterior_file(spotgenotype_file)
        plot_dir=os.path.join(args.outdir,"sf_plot")
        check_dir(plot_dir)

        spot_barcodes = list(barcode_dict.keys())

        partial_func=partial(handle_per_line,barcode_dir,plot_dir,in_name,spot_geno_file,barcode_dict,args.alpha, \
                            args.thr_r2,args.thr_prob,args.thr_likelihood,args.thr_vaf, \
                            args.plot,args.plot_supp,args.fig_size, \
                            min_array_row, max_array_row, min_array_col, max_array_col, point_size, \
                            args.method,args.num_directions,args.n_quantile,args.output_info)
        with multiprocessing.Pool(args.thread) as pool:
            results = list(tqdm(pool.imap(partial_func, mutation_identifier_list, chunksize=10), total=len(mutation_identifier_list)))
        

        out_file=open(spatial_feature_file,"w")    
        out_file.write(f'#chr\tpos\tref\talt\ttest_sig\tearly_mutation\tlate_mutation\tverylate_mutation\t'
                    'ks_stat\tks_pvalue\tmoranI_stat\tmoranI_pvalue\tmutant_rate\tmutant_rate_prob\tmutant_rate_likelihood\tmutant_rate_vaf\t'
                    'mean_vaf\tmax_vaf\tr_squared\twilcoxon_stat\twilcoxon_pvalue\toutlier_clusters\toutlier_vaf\toutlier_moranI_stat\toutlier_moranI_pvalue\n')

        for values in results:
            if values!=[]:
                out_lines = values[0]
                if out_lines != []:
                    out_file.write("\t".join(out_lines) + "\n")
        out_file.close()

        if args.output_info:
            barcode_mutation_mutation_prob_df,barcode_mutation_mutation_likelihood_df,barcode_mutation_mutation_vaf_df=pd.DataFrame(index=spot_barcodes),pd.DataFrame(index=spot_barcodes),pd.DataFrame(index=spot_barcodes)
            barcode_mutation_mutation_prob_sig_df,barcode_mutation_mutation_likelihood_sig_df,barcode_mutation_mutation_vaf_sig_df=pd.DataFrame(index=spot_barcodes),pd.DataFrame(index=spot_barcodes),pd.DataFrame(index=spot_barcodes)
            barcode_mutation_prob_list = [barcode_mutation_mutation_prob_df]
            barcode_mutation_likelihood_list = [barcode_mutation_mutation_likelihood_df]
            barcode_mutation_vaf_list = [barcode_mutation_mutation_vaf_df]
            barcode_mutation_prob_list_sig = [barcode_mutation_mutation_prob_sig_df]
            barcode_mutation_likelihood_list_sig = [barcode_mutation_mutation_likelihood_sig_df]
            barcode_mutation_vaf_list_sig = [barcode_mutation_mutation_vaf_sig_df]

            for values in results:
                if values!=[]:
                    if not values[1].empty:
                        barcode_mutation_prob_list.append(values[1])
                    if not values[2].empty:
                        barcode_mutation_likelihood_list.append(values[2])
                    if not values[3].empty:    
                        barcode_mutation_vaf_list.append(values[3])
                    if values[4]:
                        barcode_mutation_prob_list_sig.append(values[1])
                        barcode_mutation_likelihood_list_sig.append(values[2])
                        barcode_mutation_vaf_list_sig.append(values[3])

            df_lists = [barcode_mutation_prob_list, barcode_mutation_likelihood_list, barcode_mutation_vaf_list, \
                        barcode_mutation_prob_list_sig, barcode_mutation_likelihood_list_sig, barcode_mutation_vaf_list_sig]
            filenames = [
                os.path.join(args.outdir, in_name+ ".barcode_mutation_prob.txt"),
                os.path.join(args.outdir, in_name+ ".barcode_mutation_likelihood.txt"),
                os.path.join(args.outdir, in_name+ ".barcode_mutation_mutantallele.txt"),
                os.path.join(args.outdir, in_name+ ".barcode_mutation_prob_sig.txt"),
                os.path.join(args.outdir, in_name+ ".barcode_mutation_likelihood_sig.txt"),
                os.path.join(args.outdir, in_name+ ".barcode_mutation_mutantallele_sig.txt")
            ]
            parallel_merge_and_save(df_lists, filenames)

    def phase_combine():
        tmp_dir=os.path.join(str(args.outdir),"tmp_phase_combine")
        check_dir(tmp_dir)

        germ_site_bed_file=handel_candidate_informative_SNP_site(args.species,germline_file, tmp_dir)
        short_ind_file, combine_ind_germ_file=short_mosaic_sites_and_combine_with_germ(indgenotype_file, germ_site_bed_file, tmp_dir)

        if args.gtexGene:
            short_gencode_file_name=short_bed_file(args.gtexGene, short_ind_file, tmp_dir)
            gene_name_index=3
        if args.gencode:
            short_gencode_file_name=short_gencode_file(args.gencode,short_ind_file,tmp_dir,rerun=True)
            gene_name_index=4

        if args.gender in ["M","male"]:
            gender="male"
        else:
            gender="female"

        partial_func=partial(phase_combine_get_candidate_germline,ref_fasta,combine_ind_germ_file, filter_bam,args.flanking,args.minprior,gender,gene_name_index,args.run_type,args.bins)
        short_gene_review=open(short_gencode_file_name,"r")
        lines=short_gene_review.readlines()

        with multiprocessing.Pool(thread) as pool:
            results=pool.map(partial_func, lines,chunksize=1) # type: ignore # per result is [[site_info]]


        out_file=open(phase_afterUMIcombination_file,"w")
        res=[]
        for lists in results:
            if lists not in res:
                res.append(lists)
                for line in lists:
                    out_text="\t".join([str(i) for i in line])
                    out_file.write(f'{out_text}\n')

        out_file.close()
        shutil.rmtree(tmp_dir)  

    # def phase_no_combine():
    #     tmp_dir=os.path.join(str(args.outdir),"tmp_phase_no_combine")
    #     check_dir(tmp_dir)

    #     germ_site_bed_file=handel_candidate_informative_SNP_site(args.species,args.germline, tmp_dir)
    #     short_ind_file, combine_ind_germ_file=short_mosaic_sites_and_combine_with_germ(args.ind_genotype, germ_site_bed_file, tmp_dir)
    #     if args.gtexGene:
    #         short_gencode_file_name=short_bed_file(args.gtexGene, short_ind_file, tmp_dir)
    #         gene_name_index=3
    #     if args.gencode:
    #         short_gencode_file_name=short_gencode_file(args.gencode,short_ind_file,tmp_dir,rerun=True)
    #         gene_name_index=4

    #     if args.gender in ["M","male"]:
    #         gender="male"
    #     else:
    #         gender="female"

    #     partial_func=partial(phase_no_combine_get_candidate_germline,ref_fasta,combine_ind_germ_file, filter_bam,args.flanking,args.minprior,gender,gene_name_index)
    #     short_gene_review=open(short_gencode_file_name,"r")
    #     lines=short_gene_review.readlines()

    #     with multiprocessing.Pool(thread) as pool:
    #         results=pool.map(partial_func,lines,chunksize=1) # type: ignore # per result is [[site_info]]


    #     out_file=open(phase_beforeUMIcombination_file,"w")
    #     res=[]
    #     for lists in results:
    #         if lists not in res:
    #             res.append(lists)
    #             for line in lists:
    #                 out_text="\t".join([str(i) for i in line])
    #                 out_file.write(f'{out_text}\n')
            
    #     out_file.close()
    #     shutil.rmtree(tmp_dir)  

    def extract_feature():
        tmp_name=str(uuid.uuid4())
        tmpdir=os.path.join(str(args.outdir),"tmp"+tmp_name)
        check_dir(tmpdir)
        current_directory = os.path.dirname(os.path.abspath(__file__))
        compare_pl_path=os.path.join(current_directory,"others/compare_files.pl")

        gene_count_file=""
        gff3_file=args.gff3_file
        knownGene_file=""
        combine_phase_file=phase_afterUMIcombination_file
        no_combine_phase_file=""
        # annovar_annotaion_file=args.annovar_annotaion_file
        annovar_annotaion_file=""
        thread=int(args.thread)

        reference_fasta=args.fasta
        run_type=args.run_type
        bins=args.bins
        bam_file=raw_bam
        ind_count_file=args.ind_count_file
        ind_geno_file=args.ind_genotype
        h5ad_path=args.h5ad
        empty_df,adata=handle_h5ad_file(h5ad_path)
       
        vaf_cluster_file=args.vaf_cluster_file
        sample=args.outprefix
        cell_info=args.cell_info
        if cell_info :
            cell_dict=barcode_cell_mapping(cell_info)
        else:
            cell_dict={}
        
        downsample=args.downsample
        targe_dp=args.downsample_dp
        # umi_downsample_dp=args.umi_downsample_dp
        seed=args.seed
        file=open(spatial_feature_file,"r")
        lines=[]
        for line in file:
            if line[0]!="#":
                lines.append(line)
        if len(lines)==0:
            raise ValueError(f'No mosaic sites were caught!')
        
        sep="\t";grep_sf_info=1
        mode="normal"
        tmp_mutation_bed=get_mutation_bed_file(tmpdir,spatial_feature_file,sep)

        mappbablity_file=args.mappbablity_file
        used_tmp_mappbablity_file=handle_mappbablity_file(tmpdir,mappbablity_file,tmp_mutation_bed)
        
        vaf_cluster_file=args.vaf_cluster_file
        readLen=args.readLen
        prior=args.prior
        total_gene_count=""
  
        partial_func=partial(extract_feature_perline,
                            sample,
                            run_type,
                            bins,
                            tmpdir,
                            mode,
                            compare_pl_path,
                            reference_fasta,
                            grep_sf_info,
                            spatial_feature_file,
                            bam_file,
                            ind_count_file,
                            ind_geno_file,
                            gene_count_file,
                            gff3_file,
                            knownGene_file,
                            total_gene_count,
                            combine_phase_file,
                            no_combine_phase_file,
                            empty_df,
                            adata,
                            used_tmp_mappbablity_file,
                            annovar_annotaion_file,
                            vaf_cluster_file,
                            readLen,
                            barcode_dir,
                            prior,
                            cell_dict,
                            downsample,
                            targe_dp,
                            seed
                            # umi_downsample_dp
                            )
    
    
        with multiprocessing.Pool(thread) as pool:
            results=pool.map(partial_func, lines,chunksize=1) # type: ignore # per result is [[site_info]]
        

        outfile=open(feature_file,"w")
        header=""

        for _,mutation_dict in results: 
            if mutation_dict is not None and mutation_dict!={}:
                if header!="" or mutation_dict=={}:
                    pass
                else:
                    header="\t".join(list(mutation_dict.keys()))
                    outfile.write(f'#{header}\n')
                
                per_line="\t".join([str(i) for i in list(mutation_dict.values())])
                outfile.write(f'{per_line}\n')

        del results
        outfile.close()                
        
        shutil.rmtree(tmpdir)  

    
    def add_hFDR():
        if os.path.exists(artifact_signature):
            current_directory = os.path.dirname(os.path.abspath(__file__))
            script_path=os.path.join(current_directory,"others/statistic_by_mutation_signature.R")
            
            command=f"Rscript {script_path} {args.outprefix} {feature_file} {args.outdir} {args.artifact_signature} {args.reference_error_profile}"
            result=subprocess.run(command,shell=True,check=True)


    ############# Run #############
    if args.rerun:
        # get_spatial_test_result()
        phase_combine()
        extract_feature()
        add_hFDR()

    elif not check_output(spatial_feature_file):
        get_spatial_test_result()
        phase_combine()
        extract_feature()
        add_hFDR()

    elif not check_output(phase_afterUMIcombination_file):
        phase_combine()
        extract_feature()
        add_hFDR()

    elif not check_output(feature_file):
        extract_feature()
        add_hFDR()

    elif not check_output(add_hFR_feature):
        add_hFDR()
        
    else:
        print("All features have been extracted before, if you want to rerun them, please use --rerun")

    # check_output(add_hFR_feature,print_log=True)



## parameters
parser = argparse.ArgumentParser()
# required and recommend
parser.add_argument("--gender",required=False,default="female",choices=["F","M","female","male"],help="The gender of input sample, if male, chrX would not be considered")
parser.add_argument("--outdir",required=True,help="output dir")
parser.add_argument("--outprefix", help="Custom prefix for the output files.", required=False, default="sample",type=str)
parser.add_argument("--thread",required=False,default=2, type=int,help="thread")
parser.add_argument("--species",required=False,default="human", type=str,help="human species would filter hSNP, while other not")
parser.add_argument("--fasta", required=True,help="fasta file")
parser.add_argument("--filter_bam", required=False, help="bam file")
parser.add_argument("--raw_bam", required=False, help="bam file")
parser.add_argument("--germline",required=True, help="germline file")
parser.add_argument("--ind_genotype",required=True,help="individual genotype file")
parser.add_argument("--spot_genotype",required=True,help="spot genotype file")
parser.add_argument("--barcodes", required=False, help="Input a barcodes file (e.g., tissue_positions.csv from Spaceranger results). \
                                                        Example line: (barcode, in_tissue_or_not, simplified_location_x, simplified_location_y, real_location_x, real_location_y). \
                                                        'in_tissue_or_not' column: 0 means not in tissue, 1 means in tissue. This file is optional. If provided, the output count file will only include barcodes that are in tissue.")


# features related parameters
parser.add_argument("--readLen", required=False,default=150,type=int, help="read length")
parser.add_argument("--prior",required=False,default="", help="prior file")
parser.add_argument("--h5ad",required=False,default="", help="The h5ad file")
parser.add_argument("--spaceranger_result_dir",required=False,default="", help="The directory storing spaceranger result, such as demo/outs")
parser.add_argument("--ind_count_file",required=False,default="", help="ind_count_file")
parser.add_argument("--mappbablity_file",required=False,default="", help="mappbablity_file")
parser.add_argument("--gff3_file",required=False,default="", help="gff3_file")
parser.add_argument("--vaf_cluster_file",required=False,default="", help="vaf_cluster_file")
parser.add_argument("--artifact_signature",required=False,default="", help="artifact_signature")
parser.add_argument("--reference_error_profile",required=False,default="", help="the reference error profile (we provided)")


#choice one
parser.add_argument("--gtexGene",required=False, help="gene bed file downloaded from uscs goldpath")
parser.add_argument("--gencode",required=False, help="gtf file downloaded from gencode")

# optional
parser.add_argument("--downsample", required=False, default=False, action="store_true", help="downsample or not (defalut: False)")
parser.add_argument("--downsample_dp", required=False, type=int, default=2000, help="if downsample, the max dp allowed (default: 2000)")
parser.add_argument("--seed", required=False, type=int, default=42,help="ramdom seed (default: 42)")
parser.add_argument("--rerun", required=False, action="store_true", help="If set, forces the program to rerun even if the output already exists.")
parser.add_argument("--type", dest='run_type',default="visium",choices=["visium","stereo","ST","HD"],type=str, required=False, help="Your input sequence type")
parser.add_argument("--cell_info",default="",type=str, required=False, help="If you want to extract the features by cell level, please offer the relationship between barcode and cells")
parser.add_argument("--bins", required=False,default=100,type=int, help="only work for stereo-seq, if you want to combine the UMI in a bin level")
parser.add_argument("--flanking",required=False,type=int,default=120, help="read length or the flanking length")
parser.add_argument("--minprior",required=False,type=float,default=0.01, help="the min prior for filtering informative SNP")
parser.add_argument("--alpha", required=False, default=0.05, type=float, help="the significance level for the tests (default=0.05)")
parser.add_argument("--thr_r2", required=False, default=0.5, type=float, help="the threshold of R-square for distinguish early-embryonic mutation (default=0.5)")
parser.add_argument("--thr_prob", required=False, default=0.5, type=float, help="the threshold of mosaic mutation probability (default=0.5)")
parser.add_argument("--thr_likelihood", required=False, default=0.5, type=float, help="the threshold of mosaic likelihood (default=0.5)")
parser.add_argument("--thr_vaf", required=False, default=0, type=float, help="the threshold of mutant allele frequency (default=0)")
parser.add_argument("--plot", required=False, default=False, type=str2bool, help="Boolean variable of whether plotting the scatter plot of the spots in the original axes")
parser.add_argument("--plot_supp", required=False, default=False, type=str2bool, help="Boolean variable of whether plotting the supplementary plots")
parser.add_argument("--fig_size", required=False, default=5, type=int, help="Size of the output figure")
parser.add_argument("--point_size", required=False, default=None, type=int, help="Point size used in the spatial scatter plot (default=None)")
parser.add_argument("--method","-m", required=False, default="LDA", choices=["LDA", "OPT", "PCA"], help="the method used to do the dimension-reduction on which we perform the KS-test (default=LDA)")
parser.add_argument("--num_directions", required=False, default=8, type=int, help="the number of directions to test in the optimization method (OPT) (default=8)")
parser.add_argument("--n_quantile", required=False, default=100, type=int, help="the number of quantiles plotted in the QQ-plot (default=100)")
parser.add_argument("--output_info", required=False, default=True, type=str2bool, help="Boolean variable of whether output the spot mutation info file")

args = parser.parse_args()


if __name__ == '__main__':
    main()  