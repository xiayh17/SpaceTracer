import argparse
import multiprocessing
import os
import pandas as pd
from module.cluster_vaf_calculation import calculate_percluster
from module.individual_genotyper import individual_genotype
from module.read_file import read_cell_number_file, read_cluster_file, read_cluster_vaf_file, read_count_file, read_count_file_df, read_ind_posterior_file, read_prior_file
from module.spot_genotyper import spot_genotype
from utils import check_dir, check_input, check_output, str2bool


def main():
    check_dir(args.outdir)
    check_input([args.spot_count,args.cluster_count,args.ind_count,args.cluster])
    spot_count_file=args.spot_count
    cluster_count_filter_file =args.cluster_count
    ind_count_filter_file = args.ind_count

    ind_geno_file = os.path.join(args.outdir, args.outprefix +".ind_genotype.out")
    ind_geno_filter_file = os.path.join(args.outdir, args.outprefix +".ind_genotype_filter.out")
    ind_germ_file = os.path.join(args.outdir, args.outprefix +".germ_genotype.out")
    cluster_vaf_file =  os.path.join(args.outdir, args.outprefix + ".cluster_vaf.out")
    spot_genotype_file = os.path.join(args.outdir,args.outprefix +  ".spot_genotype.out")

    def get_posterior_ind():
        """get the posterior for cluster or individual"""

        prior_file=read_prior_file(args.prior)
        count_file=read_count_file(ind_count_filter_file)
        out_colname = "\t".join(["#chrom","site","ID",'germline','mutant',"cluster","spot_number","consensus_read_count", \
                                'genotype','p_mosaic','Gi','vaf'])
        out_file=open(ind_geno_file, "w")
        out_file.write(f'{out_colname}\n')
        out_germ_colname = "\t".join(["#chrom","site","genotype","allele","count","prior"])
        out_germ_file=open(ind_germ_file, "w")
        out_germ_file.write(f'{out_germ_colname}\n')

        if len(count_file)>1:
        # calculate individual and germline genotype
            for lines in count_file:
                if lines[0][0]=="#":
                    continue
                
                identifier=lines[0]+str("_")+lines[1]
                if prior_file.empty:
                    f_info=[1,1,1,1]
                try:
                    prior_info=prior_file.loc[identifier]
                    f_info=list(prior_info)[3:]
                except:
                    f_info=[0,0,0,0]

                if lines[3]=="N" or len(lines[3])!=1:
                    continue
                
                out_line=lines+f_info
                genotype_list, germline_list = individual_genotype(out_line, mu=args.mu, thr_dp=args.max_dp, pop_vaf=args.vaf, \
                                                                filter_oneallele=args.filter_oneallele)
                for genotype in genotype_list:
                    if genotype[8]!="NA":
                        out_geno="\t".join([str(i) for i in genotype])
                        out_file.write(f'{out_geno}\n')

                        out_germline="\t".join([str(i) for i in germline_list])
                        out_germ_file.write(f'{out_germline}\n')

        else:
            pass
        out_file.close()
        out_germ_file.close()
        
    def filter_pind():
        ind_raw=pd.read_csv(ind_geno_file,sep="\t")
        ind_raw['spot_number']=ind_raw['spot_number'].astype(int)
        ind_raw['vaf']=ind_raw['vaf'].astype(float)
        filter_mask = (ind_raw['vaf'] <= args.max_vaf) & (ind_raw['spot_number'] >= args.min_dp) & (ind_raw['genotype']=="mosaic")
        ind_raw[filter_mask].to_csv(ind_geno_filter_file,index=False,sep="\t")

    def get_cluster_vaf():
        """get cluster mutant allele frequency from individual genotype and cluster count"""
        # read input files
        ind_geno_file = read_ind_posterior_file(ind_geno_filter_file)
        count_file = read_count_file(cluster_count_filter_file)

        # format output file
        out_colname = "\t".join(["#chrom", "site", "ID", "germline", "mutant", "cluster", "spot_number", "consensus_read_count", "vaf"])
        out_cluster_vaf_file = open(cluster_vaf_file, "w") 
        out_cluster_vaf_file.write(f'{out_colname}\n')
            
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
            results = pool.starmap(calculate_percluster, [(cluster_count, ind_geno_file) for cluster_count in count_file])
            
        for result in results:
            if result: 
                out_cluster_vaf_file.write("\n".join(result) + "\n")

    def get_posterior_spot():
        """get spot genotype from cluster genotype"""
        # read in inputs
        ind_geno_file = read_ind_posterior_file(ind_geno_filter_file, index=False)
        ind_geno_file = ind_geno_file.drop(["cluster","spot_number","consensus_read_count"], axis=1)
        count_file = read_count_file_df(spot_count_file)
        cluster_file = read_cluster_file(args.cluster, index=False)
        cluster_vaf = read_cluster_vaf_file(cluster_vaf_file, index=False)
        cluster_vaf = cluster_vaf.drop(["spot_number","consensus_read_count"], axis=1)

        # join the data
        count_cluster_df = pd.merge(count_file, cluster_file, on='spot_barcode', how='left')
        count_geno_join = pd.merge(count_cluster_df, ind_geno_file, on=['chr','pos','ID'], how='inner')
        count_geno_join = count_geno_join.rename(columns={'vaf': 'ind_vaf'})
        count_geno_vaf_join = pd.merge(count_geno_join, cluster_vaf, on=['chr','pos','ID','germline','mutant','cluster'], how='left')
        count_geno_vaf_join = count_geno_vaf_join.rename(columns={'vaf': 'cluster_vaf'})

        # input the cell number for each spot
        if args.cellnum_file is not None:
            cell_num_file = read_cell_number_file(args.cellnum_file)
            cell_num_file = cell_num_file.drop(["cluster", "nUMI", "nREAD"], axis=1)
            count_geno_vaf_join = pd.merge(count_geno_vaf_join, cell_num_file, on='spot_barcode', how='left')

        out_colname = "\t".join(["#chrom","pos","ID","germline","mutant","cluster","spot_barcode","consensus_read_count", \
                                "l_germline", "l_mosaic", "max_spot_geno","G_spot_max","depth","vaf","p_mosaic"])
        out_file = open(spot_genotype_file, "w")
        out_file.write(f'{out_colname}\n')
        
        # calculate spot genotype
        for i in range(count_geno_vaf_join.shape[0]):
            count_geno = list(count_geno_vaf_join.iloc[i])
            # only consider the mosaic individual genotype
            if count_geno[14] != "mosaic":
                continue
            # get the cell number
            if args.cellnum_file is not None:
                cell_num = count_geno[-1]
                count_geno = count_geno[:-1]
            else:
                cell_num = args.cell_num
            # calculate spot genotype
            spot_geno_info=spot_genotype(join_info=count_geno, cell_num=cell_num, \
                                        epsQ=args.epsQ, thr_dp=args.max_dp, pop_vaf=args.min_vaf)
            if spot_geno_info[10]!="NA":
                out_spot_geno="\t".join([str(i) for i in spot_geno_info])
                out_file.write(f'{out_spot_geno}\n')
        out_file.close()

    ############# Run #############
    if args.rerun:
        get_posterior_ind()
        filter_pind()
        get_cluster_vaf()
        get_posterior_spot()
    elif not check_output(ind_geno_filter_file):
        get_posterior_ind()
        filter_pind()
        get_cluster_vaf()
        get_posterior_spot()
    elif not check_output(cluster_vaf_file):
        get_cluster_vaf()
        get_posterior_spot()
    elif not check_output(spot_genotype_file):
        get_posterior_spot()
    else:
        print("All genotyper files (spot/cell, cluster and ind) have been finished before, if you want to rerun them, please use --rerun")

    check_output(spot_genotype_file,print_log=True)
## parameters
parser = argparse.ArgumentParser()
# required and recommend
parser.add_argument("--spot_count", required=True, help="UMI counts file at spot level")
parser.add_argument("--cluster_count", required=True, help="UMI counts file at cluster level")
parser.add_argument("--ind_count", required=True, help="UMI counts file at pseudo bulk level")
parser.add_argument("--cluster", required=True, help="file contains the cluster info")
parser.add_argument("--outprefix", help="Custom prefix for the output files.", required=False, default="sample",type=str)
parser.add_argument("--outdir",required=False, default="./", help="output dir")
parser.add_argument("--type", dest='run_type',default="visium",choices=["visium","stereo","ST"],type=str, required=False, help="Your input sequence type")
parser.add_argument("--prior",required=False, default="",help="input prior file. You are not required to provide this file, and all priors of alleles would be regared as 1")
parser.add_argument("--cellnum_file",required=False, default=None, help="The file path to the file which contains the estimated cell numbers per spot")

# optional
parser.add_argument("--rerun", required=False, action="store_true", help="If set, forces the program to rerun even if the output already exists.")
parser.add_argument("--epsQ",required=False, default=20, help="the threshold for the consensus read quality (default=20)")
parser.add_argument("--alpha", required=False, default=0.05, help="the significance level for the tests (default=0.05)")
parser.add_argument("--epsAF",required=False, default=0.003,type=float, help="the threshold for the alternative allele frequency, i.e. the background error (default=0.03)")
parser.add_argument("--mu",required=False, default=1e-7, type=float, help="the population mutation rate prior (default=1e-7)")
parser.add_argument("--max_dp",required=False, default=1000,type=int, help="the max threshold for the read depth, if depth larger than the threshold, the allele numbers would be downsampled (default=1000)")
parser.add_argument("--vaf",required=False, default=1e-5, type=float, help="avoid 0 allele frequency")
parser.add_argument("--filter_oneallele",required=False, default=True, type=str2bool, help="true for delete one allele case")
parser.add_argument("--min_dp",required=False, default=30,type=int, help="the min threshold for the read depth, if depth smaller than the threshold, the allele would be removed in following analysis")
parser.add_argument("--max_vaf",required=False, default=0.3,type=float, help="the max threshold for the vaf, if vaf hihgher than the threshold, the allele would be removed in the following analysis")
parser.add_argument("--min_vaf",required=False, default=1e-5,type=float, help="the min threshold for the vaf, to avoid 0 allele frequency")
parser.add_argument("--cell_num",required=False, default=20,type=int, help="The estimated cell numbers per spot. If not given cellnum_file")


args = parser.parse_args()


if __name__ == '__main__':
    main()
