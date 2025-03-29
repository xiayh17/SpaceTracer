import pysam
from collections import defaultdict
import argparse
import os
import pandas as pd
import numpy as np

def get_cb_ub(bam_file,outname,run_type):
    # bam_file=args.bam
    # bam_file = '/storage/douyanmeiLab/yangzhirui/01.Data_download/06.skin/04.Analysis/03.bamQC/split_bam/P6_ST_vis_rep2_tissue/IN_filter.bam'  # 替换为你的 BAM 文件路径
    samfile = pysam.AlignmentFile(bam_file, "rb")

    cb_ub_stats = defaultdict(lambda: {'read_count': 0, 'ub_set': set()})

    for read in samfile:
        try:
            if run_type=="visium":
                CB=read.get_tag("CB").strip()
                UB=read.get_tag("UB").strip()

                barcode_name=str(CB)
                UMI_name=str(UB)

            elif run_type=="stereo":
                Cx=str(read.get_tag("Cx"))
                Cy=str(read.get_tag("Cy"))
                UR=read.get_tag("UR").strip()

                barcode_name=Cx+"_"+Cy
                UMI_name=str(UR)

            elif run_type=="ST":
                CB=str(read.get_tag("B0"))
                UB=str(read.get_tag("B3"))

                barcode_name=str(CB)
                UMI_name=str(UB)
            else:
                print("type",run_type)

        except:
            barcode_name=None
            UMI_name=None

        if barcode_name is not None and UMI_name is not None:
            cb_ub_stats[barcode_name]['read_count'] += 1
            cb_ub_stats[barcode_name]['ub_set'].add(UMI_name)

    samfile.close()
    
    outfile=open(outname,"w")
    outfile.write(f'barcode\tnUMI\tnREAD\n')
    for cb, stats in cb_ub_stats.items():
        outfile.write(f'{cb}\t{len(stats["ub_set"])}\t{stats["read_count"]}\n')
    outfile.close()


def main():
    if args.run:
        count_file=os.path.join(args.outdir,"raw_umi_read_count.txt")
        get_cb_ub(args.bam,count_file,args.run_type)
    else:
        if os.path.exists(args.count):
            count_file=args.count
        else:
            print("your input file not exist! Please check", args.count )

    cb_ub_df = pd.read_csv(count_file, sep='\t', header=0)
                            
    cluster_file = pd.read_csv(args.cluster, sep='\t', header=None, names=['barcode', 'cluster'])

    file_merged = pd.merge(cluster_file, cb_ub_df, on='barcode')

    cluster_sums = file_merged.groupby('cluster')['nUMI'].sum()
    max_cluster = cluster_sums.idxmax()
    max_cluster_data = file_merged[file_merged['cluster'] == max_cluster]
    median_nUMI = max_cluster_data['nUMI'].median()
    file_merged['calculated_cell_num'] = np.ceil(20 * file_merged['nUMI'] / median_nUMI)
    file_merged['refined_cell_num'] = np.where(file_merged['calculated_cell_num'] > 25, 25, file_merged['calculated_cell_num'])
    
    final_output_df = file_merged[['barcode', 'cluster', 'nUMI', 'nREAD', 'refined_cell_num']]

    output_file=os.path.join(args.outdir,"refined_umi_read_cellNum.txt")
    final_output_df.to_csv(output_file, index=False,sep="\t")

## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--bam",'-b', required=True,help="bam file")
parser.add_argument("--run", required=False,action='store_true',help="do you want to run umi count? If not, please give ")
parser.add_argument("--cluster",'-c', required=True,help="cluster file")
parser.add_argument("--type", dest='run_type',default="visium",choices=["visium","stereo","ST"],type=str, required=False, help="Your input sequence type")
parser.add_argument("--outdir",'-o',required=False,default="./", help="Output directory")
args = parser.parse_args()
    
if __name__ == '__main__':
    main()



