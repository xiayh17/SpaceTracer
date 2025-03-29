from functools import partial
import multiprocessing
from module.UMI_combine import * # type: ignore
from module.cluster_count import *
from module.cluster_count_filter import *
from module.individual_count_filter import *
from module.read_file import handle_vcf, handle_barcode
from utils import check_dir,check_input,check_output
import os
import argparse


def main():
    ############# prepare input #############
    check_input([args.bam,args.posfile])
    check_dir(args.outdir)
    mutation_identifier_list=handle_vcf(args.posfile)
    if args.barcodes!='':
        check_input(args.barcodes)
        barcode_dict=handle_barcode(args.barcodes)
    else:
        barcode_dict={}

    tmp_count_dir=os.path.join(args.outdir,"count_files");check_dir(tmp_count_dir)

    spot_count_file=os.path.join(tmp_count_dir, args.outprefix+".spot_count.out")
    cell_count_file = os.path.join(tmp_count_dir, args.outprefix +".cell.count.out") ## only for stereo-seq
    cluster_count_file = os.path.join(tmp_count_dir, args.outprefix+ ".cluster.count.out")
    cluster_count_filter_file = os.path.join(tmp_count_dir, args.outprefix + ".cluster_filter.count.out")
    ind_count_filter_file = os.path.join(tmp_count_dir, args.outprefix + ".ind_filter.count.out")
 
    ############# Function: spot count #############
    def spot_count():
        partial_func=partial(combine_all,args.bam,barcode_dict,args.platform)
        with multiprocessing.Pool(args.thread) as pool:
            result=pool.map(partial_func, mutation_identifier_list, chunksize=5)

        out_spot_count_file=open(spot_count_file, "w")
        for list_all in result:
            spot_list=list_all[0]
            if spot_list==[]:
                continue
            for one_barcode_list in spot_list:
                out_spot_count="\t".join([str(i) for i in one_barcode_list])
                out_spot_count_file.write(f'{out_spot_count}\n')

        out_spot_count_file.close()
    
    def cell_count():
        if args.platform=="stereo":
            UMI_count_cell(spot_count_file, cell_count_file, type="cell", cluster_file=args.cellpos, epsQ=args.epsQ)

    ############# Function: cluster and ind count #############
    def cluster_count():
        if args.platform=="stereo":
            UMI_count_cluster(cell_count_file, cluster_count_file, type="cluster", cluster_file=args.cellcluster, epsQ=args.epsQ)
        else:
            UMI_count_cluster(spot_count_file, cluster_count_file, type="cluster", cluster_file=args.cellcluster, epsQ=args.epsQ)
        cluster_allele_filter(cluster_count_file, cluster_count_filter_file, alpha=args.alpha, epsAF=args.epsAF)
        UMI_count_ind_from_cluster(cluster_count_filter_file, ind_count_filter_file, epsQ=args.epsQ)

    ############# Run #############
    if args.rerun:
        spot_count()
        cell_count()
        cluster_count()

    elif not check_output(spot_count_file):
        spot_count()
        cell_count()
        cluster_count()

    elif not check_output(cluster_count_file):

        cell_count()
        cluster_count()
    else:
        print("All count files (spot/cell, cluster and ind) have been finished before, if you want to rerun them, please use --rerun")

    check_output(ind_count_filter_file,print_log=True)

## parameters
parser = argparse.ArgumentParser()
# required and recommend
parser.add_argument("--posfile", help="Input a position file. Either a VCF file (e.g., chr1\t1000\t.\tA\t...) or a position file (e.g., chr1\t1000\tA\t...) are allowed. Each col is separated by tab.", required=True, type=str)
parser.add_argument("--outprefix", help="Custom prefix for the output files.", required=False, default="sample",type=str)
parser.add_argument("--bam", required=True, help="Input BAM file.")
parser.add_argument("--outdir",required=False, default="./", help="Directory where the output files will be saved.")
parser.add_argument("--cellpos", required=False, help="Required only for stereo-seq. Provide a file containing the relationship between positions and cells, e.g., '1000_1000\tcell1'.")
parser.add_argument("--cellcluster", required=False, help="Provide a file containing the relationship between spots/cells and cluster, e.g., 'cell1\tclusterA'.")

# optional
parser.add_argument("--rerun", required=False, action="store_true", help="If set, forces the program to rerun even if the output already exists.")
parser.add_argument("--barcodes", required=False, help="Input a barcodes file (e.g., tissue_positions.csv from Spaceranger results). \
                                                        Example line: (barcode, in_tissue_or_not, simplified_location_x, simplified_location_y, real_location_x, real_location_y). \
                                                        'in_tissue_or_not' column: 0 means not in tissue, 1 means in tissue. This file is optional. If provided, the output count file will only include barcodes that are in tissue.")
parser.add_argument("--platform",default="visium",choices=["visium","stereo","ST"],type=str, required=False, help="Specify the platform or sequencing type used for input data. Choices are: 'visium', 'stereo', 'ST'.")
parser.add_argument("--thread",required=False, default=2,type=int, help="The number of threads to use. Ensure this number matches the number of CPU cores available." )
parser.add_argument("--epsQ",required=False, default=20, help="Threshold for consensus read quality. Reads with quality below this value will be excluded (default=20).")
parser.add_argument("--alpha", required=False, default=0.05, help="The significance level (alpha) for statistical tests. Default is 0.05.")
parser.add_argument("--epsAF",required=False, default=0.01,type=float, help="Threshold for alternative allele frequency (background error rate). Variants with frequency below this value will be ignored (default=0.01).")
args = parser.parse_args()

if __name__ == '__main__':
    main()


