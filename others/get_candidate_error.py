

from functools import partial
import multiprocessing
import argparse
import multiprocessing
from module.read_file import handle_vcf, handle_barcode, read_bam
from module.UMI_combine import *
from utils import split_bam_by_chrom,check_dir, split_vcf
import os
import argparse
import gc
import pandas as pd

def combine_bulk_for_errors(bam,barcode_dict,threshold,run_type,mv_chr,identifier):
    #spot_info=combine_UMI_spot_tidy(bam_dict,barcode_dict,identifier)
    # print("type:",run_type)
    
    if identifier.split("_")[0] in mv_chr:
        return []

    else:
        print(identifier)
        error_info=combine_UMI_bulk_for_errors_tidy(bam,barcode_dict,threshold,run_type,identifier)
        return [error_info]


def combine_UMI_bulk_for_errors_tidy(bam_file,barcode_dict,threshold,run_type,identifier):
    '''
    identifier:"chr1_1000_A_alt"
    '''
    identifier_list=identifier.split("_")
    chr=identifier_list[0]
    
    pos=int(identifier_list[1])
    ref=identifier_list[2]
    # print("\n======",chr,bam_file)
    handel_bam=read_bam(bam_file)
    site_barcode_UMI_dict,strand=handle_reads_per_pos_read_count_and_strand(handel_bam,chr,pos,run_type)
    if barcode_dict=={}:
        filter_site_barcode_UMI_dict=site_barcode_UMI_dict
    else:
        filter_site_barcode_UMI_dict={}
        for key in site_barcode_UMI_dict.keys():
            if key.split(".")[0] in list(barcode_dict.keys()):
                filter_site_barcode_UMI_dict[key]=site_barcode_UMI_dict[key]

    error_info=combine_UMI_bulk_for_errors(filter_site_barcode_UMI_dict,chr,pos,ref,threshold) + [strand]

    return error_info


def combine_UMI_bulk_for_errors(site_barcode_UMI_dict,chrom,pos,ref,threshold):
    #out_df = pd.DataFrame(columns=["chr","pos","ref","alt","spot_number","consensus_read_count","read_quality_A","read_quality_T","read_quality_C","read_quality_G"])
    #spot_number=len(site_barcode_UMI_dict.keys())
    pcr_errors,lysis_errors=0,0
    pcr_alts,lysis_alts=[],[]
    UMI_dp=0
    dp=0
    # lysis_UMI_count=defaultdict(int)
    for barcode in site_barcode_UMI_dict.keys():
        UMI_dp+=len(site_barcode_UMI_dict[barcode].keys())
        for UMI in site_barcode_UMI_dict[barcode]:
            count_dict=site_barcode_UMI_dict[barcode][UMI]["count"]
            dp+=sum(count_dict.values())
            #quality_dict=site_barcode_UMI_dict[barcode][UMI]["quality"]
            pcr_error, lysis_error,pcr_alt,lysis_alt=check_errors(count_dict,ref,threshold)
            pcr_errors += pcr_error; lysis_errors += lysis_error
            pcr_alts += pcr_alt; lysis_alts += lysis_alt

    if lysis_errors >1 or lysis_errors==0:
        lysis_alt_out="."   
    else:
        lysis_alt_out=lysis_alts[0]

    pcr_alt_out=",".join(list(set(pcr_alts)))
    # lysis_alt_out=",".join(list(set(lysis_alts)))
    if pcr_alt_out =="":
        pcr_alt_out="." 
    # if lysis_alt_out =="":
    #     lysis_alt_out="."     

    error_info=[chrom, pos, dp,UMI_dp,ref, pcr_alt_out, pcr_errors, lysis_alt_out,lysis_errors]

    return error_info


def handle_reads_per_pos_read_count_and_strand(bam_handle,chrom,pos,run_type):
    '''
    input:
    bam_handel: the initial bam file handeled by pysam
    
    output:
    dict: {barcode_name: {UMI_name: {"count": 1, "quality": {"A":{30:9, 10:1}, "T":{10:1}}}}}
    '''
    reads=bam_handle.fetch(chrom,pos-1,pos,multiple_iterators=True)
    pos_index = pos-1

    site_barcode_UMI_dict={}
    reverse_dp=0
    forward_dp=0
    for item in reads:
        # all_reads_count += 1
        try:
            # a part of reads didn't have the information of "CB", bacause the "CR" didn't pass QC
            if run_type=="visium":
                CB=item.get_tag("CB").strip()
                UB=item.get_tag("UB").strip()

                barcode_name=str(CB)
                UMI_name=str(UB)

            elif run_type=="stereo":
                Cx=str(item.get_tag("Cx"))
                Cy=str(item.get_tag("Cy"))
                UR=item.get_tag("UR").strip()

                barcode_name=Cx+"_"+Cy
                UMI_name=str(UR)
            elif run_type=="ST":
                CB=str(item.get_tag("B0"))
                UB=str(item.get_tag("B3"))

                barcode_name=str(CB)
                UMI_name=str(UB)
            else:
                # print("type",run_type)
                continue
                
        except:
            # barcode_low_quality += 1
            # print("no CB UB for ", item)
            continue
        
        try:
            item.get_reference_positions().index(pos_index)
        except:
            # skip_reads += 1
            # print("no pos_index for ", item)
            continue
        # mapping_reads += 1

        if item.is_reverse in [True,"TRUE","true","True"]:
            reverse_dp+=1
        else:
            forward_dp+=1        
        
        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)

        if pos_index in cut_pos:
            # effective_DP += 1
            geno = cut_seq[cut_pos.index(pos_index)]
            if geno not in "ATCG":
                print("not ATCG ", item)
                continue

            raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
            try:
                qualities=item.get_forward_qualities()
                quality=qualities[raw_index]
            except:
                # qualities=None
                print("not qualities ", item)
                # i+=1
                continue

            if barcode_name not in site_barcode_UMI_dict.keys():
                site_barcode_UMI_dict[barcode_name]=defaultdict(dict)

            if UMI_name not in site_barcode_UMI_dict[barcode_name].keys():
                site_barcode_UMI_dict[barcode_name][UMI_name]["count"]=defaultdict(int)
                site_barcode_UMI_dict[barcode_name][UMI_name]["quality"]={"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
 
            site_barcode_UMI_dict[barcode_name][UMI_name]["count"][geno]+=1
            site_barcode_UMI_dict[barcode_name][UMI_name]["quality"][geno][quality]+=1
    
    # print(i)
    #stat_info = [all_reads_count, effective_DP, barcode_low_quality, mapping_reads, skip_reads]
    #print(stat_info)
    if reverse_dp>=forward_dp:
        major_read_strand="-"
    elif reverse_dp<forward_dp:
        major_read_strand="+"
    else:
        major_read_strand="unknown"
        
    return site_barcode_UMI_dict, major_read_strand


def UMI_combination_for_errors(args):
    # print("yes")
    work_dir=os.getcwd()
    check_dir(args.outdir)
    out_dir=os.path.join(args.outdir, "split_bam_dir")
    split_vcf_dir=os.path.join(args.outdir, "split_vcf_dir")
    # check_dir(out_dir); check_dir(split_vcf_dir)
    if not os.path.exists(args.bam):
        print("The bam file " , args.bam, "does not exist!")
    #split_bam_dict=split_bam_by_chrom(args.bam,out_dir,work_dir,run_index=False) ## if you do not need to run the split bam, please input no
    #vcf_file_dict=handel_vcf_for_bed(args.vcffile,split_vcf_dir)
    mutation_identifier_list=handle_vcf(args.vcffile)
    #in_name=args.vcffile.split("/")[-1].split(".")[0]
    in_name=args.sample 
    thread=int(args.thread)
    # print("stop")

    print(f"The following steps will be runned in {thread} thread.")
    new_list=[]
    for mutation in mutation_identifier_list:
        # print(mutation)
        new_list.append(mutation)

    if args.outsuffix=="":
        outsuffix="."
    else:
        outsuffix="."+args.outsuffix+"."

    if args.barcodes:
        barcode_dict=handle_barcode(args.barcodes)
    else:
        barcode_dict={}
    mv_chrom=["chrM","MT"]
    bam=args.bam
    partial_func=partial(combine_bulk_for_errors,bam,barcode_dict,args.threshold,args.run_type,mv_chrom)
    # with multiprocessing.Pool(thread) as pool:
    #     # result_spot=pool.map(partial_func_spot, new_list, chunksize=1000) # type: ignore # per result is [[site_info]]
    #     # result_ind=pool.map(partial_func_bulk, new_list, chunksize=1000)
    #     result=pool.map(partial_func, new_list, chunksize=5)

    COLUMNS = ['#chr', 'pos', 'dp','UMI_dp','ref', 'pcr_alt_out', 'pcr_errors', 'lysis_alt_out','lysis_errors','strand']
    out_error_file=os.path.join(args.outdir, in_name + outsuffix+ "error.count.out")
    with multiprocessing.Pool(thread) as pool, open(out_error_file, "w") as f:
        pd.DataFrame(columns=COLUMNS).to_csv(f, header=True, index=False,sep='\t')
        
        for result in pool.imap(partial_func, new_list, chunksize=1000):
            if result: 
                for per_result in result: 
                    chunk_df = pd.DataFrame([per_result], columns=COLUMNS)
                    chunk_df.to_csv(f, header=False, mode="a", index=False,sep='\t')
            del result
            gc.collect()


if __name__ == '__main__':
    ## parameters
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')
    
    # combine_UMI for error in bulk level
    parser_umi_error = subparsers.add_parser('error', help='combine UMI for ind ')
    parser_umi_error.add_argument("--bam","-b", required=True, help="input bam file")
    parser_umi_error.add_argument("--sample","-s", required=True, help="sample name")
    parser_umi_error.add_argument("--barcodes", required=False, help="inut gz barcodes file; if yes, the output count file will only contain barcodes in given file")
    parser_umi_error.add_argument("--posfile",required=False, help="input mutation information POS file; chr pos ref alt")
    parser_umi_error.add_argument("--vcffile",required=False, help="input mutation information VCF file; chr pos ref alt")
    parser_umi_error.add_argument("--type", dest='run_type',default="visium",choices=["visium","stereo","ST"],type=str, required=False, help="Your input sequence type")
    parser_umi_error.add_argument("--outdir",required=False, default="./", help="output dir")
    parser_umi_error.add_argument("--outsuffix",required=False, default="", help="output name suffix")
    parser_umi_error.add_argument("--threshold",required=False, default=3,type=int, help="threshold to filter errors" )
    parser_umi_error.add_argument("--thread",required=False, default=2,type=int, help="the thread you want to use, please make sure thread is equal to the cpu numbers you used" )
    parser_umi_error.set_defaults(func=UMI_combination_for_errors)

    args = parser.parse_args()
    args.func(args)