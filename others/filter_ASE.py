import argparse
from collections import defaultdict
from functools import reduce
from math import ceil, log10
import os
import subprocess
import pysam
import pandas as pd
from pybedtools import BedTool
import numpy as np
from module.UMI_combine import calculate_UMI_combine_phred, get_most_candidate_allele
from utils import check_dir
from scipy.stats import binom
import statsmodels.stats.multitest as smm

def get_intergration_from_identifier_and_file(compare_pl_path,identifier,query_file,index_in_query="0 1 2 3"):
    """
    This function is used to get the intergration result from identifier and one special query file.
    Note: The top 4 columns of query_file must be "chrom pos ref alt"
    
    """
    perl_script=compare_pl_path
    identifier_line="\t".join(identifier.strip().split("_"))

    command=f"echo -e \"{identifier_line}\" |perl {perl_script} - {query_file} {index_in_query}|sort -u"
    try:
        result=subprocess.check_output(command,text=True,shell=True)
        return result
    except:
        print(f"Something wrong when run the command: {command}")
        return ""


def handle_barcode(barcode_file,pos=0,in_tissue_choose=0):
    '''
    input: (barcode,in_tissue_or_not,simplified_location_x,simplified_location_y,real_location_x,real_location_y) #in_tissue_or_not:0 means not in tissue, 1 means in tissue
    ACGCCTGACACGCGCT-1,0,0,0,323,308
    TACCGATCCAACACTT-1,0,1,1,334,326
    
    argrument:
    pos means which kind of position you want to use, if you given a int not equal to 0, the simplified location will be used.

    output:
    dict: {barcode: (in_tissue, pos1, pos2),...}
    '''
    barcode_dict={}
    f=open(barcode_file,"r")
    for line in f.readlines():
        s = line.strip().split(",")
        barcode=s[0]
        if barcode!="barcode":
            in_tissue=int(s[1])
            if pos==0:
                pos1= int(s[2]); pos2=int(s[3])
            else:
                pos1= int(s[4]); pos2=int(s[5])
            
            if in_tissue_choose==0 and in_tissue==1:
                barcode_dict[barcode]= (in_tissue, pos1, pos2)
            elif in_tissue_choose==1:
                barcode_dict[barcode]= (in_tissue, pos1, pos2)
            else:pass
        else: pass
    return barcode_dict


def add_prior(ref,alt,prior_info):
    infos=prior_info.strip().split("\t")
    print(infos)
    # print("ATCG".index(self.ref))
    if "," not in ref:
        fref=str(infos[3+"ATCG".index(ref)])
    else:
        fref=infos[3+"ATCG".index(ref[0])] + "," + infos[3+"ATCG".index(ref[-1])]
    falt=infos[3+"ATCG".index(alt)]
    return fref,falt

def main():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    compare_pl_path=os.path.join(os.path.dirname(current_directory),"others/compare_files.pl")

    sig_df=pd.read_csv(args.features,sep="\t")
    sig_df['identifier']
    identifier_list=sig_df['identifier'].tolist()
    sig_df=sig_df.set_index('identifier')

    outdir=args.outdir
    check_dir(outdir)

    if args.prefix!="":
        prefix=args.prefix+"."
    else:
        prefix=""

    tmp_identifier_file=os.path.join(outdir,prefix+"tmp_identifier.txt")
    with open(tmp_identifier_file,"w") as f :
        f.write("\n".join(identifier_list))
    f.close()

    # bam_path="/storage/douyanmeiLab/yangzhirui/01.Data_download/01.Skin_Cancer/04.Analysis/01.cellranger_scrna_hg38"

    bed_file=os.path.join(outdir,prefix+"tmp.bed")
    run_shell="awk -F'_' '{print $1, $2, $2, $3, $4}' OFS='\t'  %s > %s" % (tmp_identifier_file, bed_file)
    result=subprocess.run(run_shell,shell=True)
    if result.returncode!=0:
        print(f'Something wrong when identifier tmp bed file for {bed_file}.')

    ase_file=args.ase_file
    # ase_file="/storage/douyanmeiLab/yangzhirui/01.Data_download/07.brain/04.Analysis/04.mutations/151673/candidate_ASE_sites.txt"
    ase_tmp_file=os.path.join(outdir,prefix+prefix+"ase.tmp.bed")
    run_shell = "sed '1d' %s |awk '{print $10, $11, $12, $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' - | sort -u | \
                 bedtools intersect -a - -b %s -wa > %s" % (ase_file,bed_file, ase_tmp_file)
    result=subprocess.run(run_shell,shell=True)
    print("running...")
    if result.returncode!=0:
        print(f'Something wrong when ase tmp file for {ase_tmp_file}.')
    print("finishe ase tmp bed")

    if args.RNA_editing!="":
        RNA_editing=args.RNA_editing
        RNA_editing_tmp_file=os.path.join(outdir,prefix+"RNA_editing.tmp.bed")
        run_shell = "bedtools intersect -a %s -b %s -wa > %s" % (RNA_editing,bed_file, RNA_editing_tmp_file)
        result=subprocess.run(run_shell,shell=True)
        if result.returncode!=0:
            print(f'Something wrong when RNA editing tmp file for {RNA_editing_tmp_file}.')
        
        print("editing bed")

    threshold=0.05
    sig_df["ase"]="no"
    sig_df["editing"]="no"

    if 'fref' not in sig_df.columns:
        sig_df["fref"]="no"
        sig_df["falt"]="no"


    for identifier in identifier_list:
        chrom,pos,ref,alt=identifier.split("_")
        print(identifier)
        #countsum=int(sig_df.loc[[identifier]]["consensus_UMI_count"].iloc[0])
        #alt_count=int(sig_df.loc[[identifier]]["consensus_alt_allele_count"].iloc[0])
        countsum=int(sig_df.loc[[identifier]]["consensus_UMI_count"].iloc[0])
        alt_count=int(sig_df.loc[[identifier]]["consensus_alt_allele_count"].iloc[0])

        identifier_bed=BedTool("\t".join([chrom,pos,pos,ref,alt]),from_string=True)
        ase_tmp_bed=BedTool(ase_tmp_file)

        pos_intersect = identifier_bed.intersect(ase_tmp_bed, wb=True)
        print(pos_intersect)
        if pos_intersect.count()==0:
            ase="Unknown"
        else:
            ase_list=[]
            directions=[]
            for line in pos_intersect:
                new_sline=line.fields
                # print(new_sline)
                count1,count2=new_sline[12].split(",") 
                germ_vaf=int(count1)/(int(count1)+int(count2))
                ase_p_value=float(new_sline[14])
                if germ_vaf>=0.35 and germ_vaf<=0.65:
                    directions.append("balance")

                mosaic_p_value =  binom.cdf(int(alt_count), countsum, p=min(germ_vaf,(1-germ_vaf)))
                ase_list.append(mosaic_p_value)

            balance_count = directions.count("balance")

            #if balance_count == len(directions):
            #ase="False"
            #else:
            _, p_adj, _, _ = smm.multipletests(ase_list, method='fdr_bh')
                
            if max(p_adj) >= threshold:
                ase="True"
            else:
                ase="False"
            print(identifier,ase)
            # if result in ["ref","alt"] and sum(ase_list)>=1:
            #     ase="True"
            # elif result in ["balance"] or sum(ase_list)==0:
            #     ase="False"
            # else:
            #     ase="Unclear"

        if args.RNA_editing!="":
            RNA_editing_bed=BedTool(RNA_editing_tmp_file)
            editing_intersect = identifier_bed.intersect(RNA_editing_bed, wb=True)
            if editing_intersect.count()!=0:
                editing="True"
            else:
                editing="False"

            sig_df.loc[[identifier], 'ase'] = ase
            sig_df.loc[[identifier], 'editing'] = editing

        if args.prior!="" and sig_df.loc[[identifier]]['fref'].iloc[0]=="no":
            only_pos_identifier="\t".join([chrom,str(pos),chrom,str(pos)])
            prior=args.prior
            prior_info=get_intergration_from_identifier_and_file(compare_pl_path,only_pos_identifier,prior,index_in_query="0 1 0 1")
            print(prior_info)
            if prior_info!="":
                fref,falt=add_prior(ref,alt,prior_info)
                sig_df.loc[[identifier], 'fref'] = fref
                sig_df.loc[[identifier], 'falt'] = falt

    sig_df.to_csv(outdir+"/"+prefix+"features_add_ase_editing.txt",index=True,sep="\t")
    
    sig_df['ase']=sig_df['ase'].astype(str)
    sig_df['editing']=sig_df['editing'].astype(str)
    filtered_df2 = sig_df[
        ((sig_df['ase'] == 'False') | (sig_df['ase'] == 'Unknown')) & 
        (sig_df['editing'] == 'False')
    ]
    filtered_df2.index.to_frame(index=False).to_csv(outdir+"/"+prefix+"filter_ase.filter_editing.identifier.txt",index=False,header=False,sep="\t")
    filtered_df2.to_csv(outdir+"/"+prefix+"filter_ase.filter_editing.feature.txt",index=True,sep="\t")


## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--features", required=True, help="feature csv file")
parser.add_argument("--RNA_editing", default="",required=False, help="RNA editing bed file")
parser.add_argument("--ase_file", default="", required=True, help="ase file")
parser.add_argument("--prior",required=False,default="", help="prior file")
parser.add_argument("--outdir","-o",required=False, default="IGV_Plot", help="output dir")
parser.add_argument("--prefix",required=False, default="", help="output file prefix")

args = parser.parse_args()


if __name__ == '__main__':
    main()
