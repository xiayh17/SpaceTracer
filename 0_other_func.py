import argparse
import datetime
from functools import partial
import multiprocessing
import os
import shutil
from module.read_file import get_freq_for_mutation_simple_from_annovar, handel_and_split_pos_for_bed, handle_prior_bed
from others.domain_segmentation import cluster2domain, cluster_spot_GraphST, cluster_spot_SpaGCN
from utils import check_dir, check_input, check_output, get_chrom_list_from_file, intersect_pos_gnomad, str2bool


def get_prior(args):
    vcf2bed_split_dir=os.path.join(args.outdir,"split_vcf")
    short_gnomad_dir=os.path.join(args.outdir,"short_gnomad")
    check_dir(args.outdir);check_dir(vcf2bed_split_dir);check_dir(short_gnomad_dir)
    
    check_input([args.posfile,args.annovar])
    pos_bed_file_dict=handel_and_split_pos_for_bed(args.posfile,vcf2bed_split_dir)
    in_name=args.outprefix
    thread=args.thread
    out_prior=open(os.path.join(args.outdir, in_name + ".prior.out"), "w")
    gnomad_dict=get_chrom_list_from_file(args.annovar)
    # print(pos_bed_file_dict)
    # print(gnomad_dict)
    short_gnomad_dict=intersect_pos_gnomad(pos_bed_file_dict, gnomad_dict, short_gnomad_dir)
    # print(short_gnomad_dict)
    partial_func=partial(get_freq_for_mutation_simple_from_annovar,short_gnomad_dict)

    for chrom in pos_bed_file_dict.keys():
        bed_file=pos_bed_file_dict[chrom]
        new_list=handle_prior_bed(bed_file)
        print("read_vcf:"+chrom,datetime.datetime.now())
        with multiprocessing.Pool(thread) as pool:
            result=pool.map(partial_func, new_list, chunksize=1000) 

        for line in result:
            out="\t".join([str(i) for i in line])
            out_prior.write(f'{out}\n')

    out_prior.close()
    check_output(os.path.join(args.outdir, in_name + ".prior.out"))
    shutil.rmtree(vcf2bed_split_dir);shutil.rmtree(short_gnomad_dir)  
    

def ase(args):
        
    def binomial_test(row):
        h1, h2 = map(int, row['count'].split(','))
        p_value = binom_test(h1, h1 + h2, p=0.5, alternative='two-sided')
        return p_value

    current_directory = os.path.dirname(os.path.abspath(__file__))

    germline_sites=pd.read_csv(args.germline,sep="\t",header=None,names=["#chrom","site","genotype","allele","count","prior"])
    germline_sites=germline_sites[germline_sites["genotype"]=="het"]
    germline_sites['count_sum'] = germline_sites['count'].apply(lambda x: sum(map(int, x.split(','))))
    germline_sites['prior_valid'] = germline_sites['prior'].apply(lambda x: all(float(val) > args.prior_threshold for val in x.split(',')))
    filtered_df = germline_sites[(germline_sites['count_sum'] > args.count_threshold) & (germline_sites['prior_valid'])]
    filtered_df['binomial_p_value'] = filtered_df.apply(binomial_test, axis=1)
    filtered_df=filtered_df[filtered_df['binomial_p_value']<=args.p_threshold]

    filtered_df.rename(columns={"#chrom": "chrom"}, inplace=True)
    filtered_df.rename(columns={"site": "pos"}, inplace=True)
    germ_result_file=os.path.join(args.outdir,"temp_filter_germ.txt")
    filtered_df[['ref', 'alt']] = filtered_df['allele'].str.split(',', expand=True)
    filtered_df[['chrom','pos','ref', 'alt']].to_csv(germ_result_file,sep="\t",index=False)

    compare_pl_path=os.path.join(os.path.dirname(current_directory),"others/compare_files.pl")
    vcf_file=args.dbSNP
    if vcf_file !="":
        index_in_query="0 1 3 4"
        filter_germ_file=germ_result_file
        command=f"grep -v '#' {vcf_file} | perl {compare_pl_path} {filter_germ_file} - {index_in_query}"
        result=subprocess.run(command,shell=True,check=True, stdout=subprocess.PIPE, text=True)
        df = pd.read_csv(StringIO(result.stdout), delim_whitespace=True,header=None,names=["chrom","pos","dbSNPID","ref","alt","info","filter","detail"])
        filtered_df['pos'] = filtered_df['pos'].astype('int64')
        df['pos'] = df['pos'].astype('int64')
        merged_df = pd.merge(filtered_df, df[['chrom', 'pos', 'ref', 'alt', 'dbSNPID']], 
                            on=['chrom', 'pos', 'ref', 'alt'], how='left')

        merged_df=merged_df[pd.notna(merged_df['dbSNPID'])]

    else:
        merged_df["dbSNPID"]="NA"

    gtex_bed=BedTool(args.gtexGene)
    def gene_detect(row):
        gene_name=[];start_list=[];end_list=[]
        chrom = row['chrom']
        site = row['pos']
        gene_line="\t".join([chrom,str(site),str(site)])
        gene_bed=BedTool(gene_line,from_string=True)
        pos_intersect = gtex_bed.intersect(gene_bed,u=True)
        region_chrom="NA"
        for line in pos_intersect:
            new_sline=line.fields
            gene_name.append(new_sline[3])
            region_chrom=new_sline[0]
            start_list.append(new_sline[1])
            end_list.append(new_sline[2])

        if start_list!=[]:
            region_start=min([int(i) for i in start_list])
        else:
            region_start="NA"
        
        if end_list!=[]:
            region_end=max([int(i) for i in end_list])
        else:
            region_end="NA"

        gene=",".join(gene_name)

        if gene=="":
            gene="NA"
            region_chrom=chrom
            region_start=site-150
            region_end=site+150

        return pd.Series([gene, region_chrom, region_start, region_end])

    merged_df[['gene', 'region_chrom', 'region_start', 'region_end']] = merged_df.apply(gene_detect, axis=1)
    short_df=merged_df[['chrom', 'pos', 'ref', 'alt','count', 'prior','binomial_p_value','dbSNPID','gene','region_chrom','region_start','region_end']]
    short_df.to_csv(os.path.join(args.outdir,args.outname),sep="\t",index=False)


    try:
        os.remove(germ_result_file)
    except:
        pass

def cluster_spot(args):
    """
    Cluster the spots
    """
    check_dir(args.outdir)
    if args.method == "SpaGCN":
        # input files from indir
        h5_file_name = "raw_feature_bc_matrix.h5"
        spatial_file_name = "spatial/tissue_positions.csv"
        image_file_name = "spatial/tissue_hires_image.png"
        # run SpaGCN
        cluster_spot_SpaGCN(args.indir, h5_file_name, spatial_file_name, image_file_name, args.outdir, args.sample, args.ncluster, \
                            plot=args.plot, init_cluster=args.init_method, s=args.weight_histology, b=args.b, histology=args.histology, \
                            p=args.pertentage, seed=args.seed, tol=args.tol, lr=args.lr, max_epochs_run=args.max_epochs_run, \
                            l_start=0.01, l_end=1000, l_tol=0.01, l_max_run=100, \
                            res_start=0.7, res_step=0.1, res_tol=5e-3, res_lr=0.05, res_max_epochs=20)
        # divide seperate cluster into domains
        cluster2domain(args.outdir, args.sample, plot=args.plot, distance_threshold=args.distance_threshold, \
                       min_samples=args.min_samples, num_threshold=args.num_threshold, shape="hexagon", keep=False)
        
    if args.method == "GraphST":
        # run GraphST
        cluster_spot_GraphST(args.indir, args.outdir, args.sample, args.ncluster, type=args.type, \
                            h5_file_name=args.h5_file_name, \
                            R_HOME=args.R_HOME, radius=args.radius, tool=args.tool, refinement=args.refinement, plot=args.plot)
        if args.type == "Visium":
            # divide seperate cluster into domains
            cluster2domain(args.outdir, args.sample, plot=args.plot, distance_threshold=args.distance_threshold, \
                        min_samples=args.min_samples, num_threshold=args.num_threshold, shape="hexagon", keep=False)
     

if __name__ == "__main__":   
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    # get prior
    parser_prior = subparsers.add_parser('prior', help='get prior from gnomad file')
    parser_prior.add_argument("--outprefix", help="Custom integer value for outprefix. The prefix of your output", required=True, type=str)
    parser_prior.add_argument("--posfile",required=False, help="input mutation information POS file; chr pos ref alt")
    parser_prior.add_argument("--outdir",required=False, default="./", help="output dir")
    parser_prior.add_argument("--annovar",required=False, help="gnomad list file from annovar containing all gnomad file, or one specific file")
    parser_prior.add_argument("--thread",required=False, default=2,type=int, help="the thread you want to use, please make sure thread is equal to the cpu numbers you used" )
    parser_prior.set_defaults(func=get_prior)

    # cluster spots
    parser_cluster = subparsers.add_parser('cluster', help='cluster spots') 
    parser_cluster.add_argument("--indir","-d", required=True, help="the directory of data input")
    parser_cluster.add_argument("--outdir",required=True, help="output dir")
    parser_cluster.add_argument("--sample","-s", required=True, help="sample name")
    parser_cluster.add_argument("--ncluster","-n", required=True, default=5, type=int, help="the number of clusters")
    parser_cluster.add_argument("--plot","-p", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for plotting the cluster results")
    parser_cluster.add_argument("--method", required=False, default="SpaGCN", choices=['SpaGCN', 'GraphST'], help="clustering method (SpaGCN or GraphST, default = SpaGCN)")
    parser_cluster.add_argument("--init_method", required=False, default="louvain", choices=["louvain", "kmeans"], help="initial clustering method (louvain or kmeans)")
    parser_cluster.add_argument("--distance_threshold", required=False, default=2, type=int, help="the maximum threshold for the distance between two spots which could be included in one domain (default = 2)")
    parser_cluster.add_argument("--min_samples", required=False, default=1, type=int, help="the minimum number of spots in one domain when seperation (default = 1)")
    parser_cluster.add_argument("--num_threshold", required=False, default=30, type=int, help="the minimum number of spots in one domain that would be left (default = 30)")
    parser_cluster.add_argument("--weight_histology", required=False, default=1, type=int, help="the weight given to histology when calculating Euclidean distance between every two spots. (default = 1). \
                                value 1 means that the histology pixel intensity value has the same scale variance as the (x,y) coordinates, whereas higher value indicates higher weight to histology")
    parser_cluster.add_argument("--b", required=False, default=49, type=int, help="the area of each spot when extracting color intensity (default = 49)")
    parser_cluster.add_argument("--histology", required=False, default=True, type=str2bool, choices=[True, False], help="whether use the histology image for calculating the adjacency matrix (default = True)")
    parser_cluster.add_argument("--pertentage", required=False, default=0.5, type=float, help="Percentage of total expression contributed by neighborhoods (default = 0.5)")
    parser_cluster.add_argument("--seed", required=False, default=100, type=int, help="the random seed (default = 100)")
    parser_cluster.add_argument("--tol", required=False, default=5e-3, type=float, help="the tolerance (default = 5e-3)")
    parser_cluster.add_argument("--lr", required=False, default=0.05, type=float, help="the learning rate (default = 0.05)")
    parser_cluster.add_argument("--max_epochs_run", required=False, default=200, type=int, help="the maximum number of epochs when running SpaGCN (default = 200)")
    parser_cluster.add_argument("--type", required=False, default='Visium', choices=['Visium', 'Stereo'], help="the data type ('Visium' or 'Stereo', default = 'Visium')")
    parser_cluster.add_argument("--h5_file_name", required=False, default='filtered_feature_bc_matrix.h5', help="the h5 file name (default = 'filtered_feature_bc_matrix.h5', not need to modify for 10x Visium data)")
    parser_cluster.add_argument("--R_HOME", required=False, default=None, help="R installation path (default = None)")
    parser_cluster.add_argument("--radius", required=False, default=6, type=int, help="the number of neighbors considered during refinement for GraphST (default = 6)")
    parser_cluster.add_argument("--tool", required=False, default='louvain', help="clustering method for GraphST ('mclust', 'leiden', and 'louvain', default = 'louvain')")
    parser_cluster.add_argument("--refinement", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable of whether refine the clustering result for GraphST (refinement = True)")
    parser_cluster.set_defaults(func=cluster_spot)

    # ase
    parser_ase = subparsers.add_parser('ase', help='get prior from gnomad file')
    parser_ase.add_argument("--outprefix", help="Custom integer value for outprefix. The prefix of your output", required=True, type=str)
    parser_ase.add_argument("--posfile",required=False, help="input mutation information POS file; chr pos ref alt")
    parser_ase.add_argument("--outdir",required=False, default="./", help="output dir")
    parser_ase.add_argument("--annovar",required=False, help="gnomad list file from annovar containing all gnomad file, or one specific file")
    parser_ase.add_argument("--thread",required=False, default=2,type=int, help="the thread you want to use, please make sure thread is equal to the cpu numbers you used" )
    parser_ase.set_defaults(func=ase)

    args = parser.parse_args()
    args.func(args)

