import collections
import multiprocessing
import subprocess
from pybedtools import BedTool
import os
from pickle import TRUE
import re
import sys
from math import exp, floor, log, log10,ceil,inf
import argparse
from collections import Counter, defaultdict
import numpy as np

######### check #########
def check_input(files):
    if type(files)==list:
        for file in files:
            if file is None:
                raise ValueError(f"Wrong input! Please offer '{file}'")
            elif not os.path.exists(file):
                raise FileNotFoundError(f"Wrong input! '{file}' not found!")
    elif type(files)=="str":
        if file is None:
            raise ValueError(f"Wrong input! Please offer '{file}'")
        elif not os.path.exists(file):
            raise FileNotFoundError(f"Wrong input! '{file}' not found!")


def check_output(file_path,print_log=False):
    if not os.path.exists(file_path):
        if print_log:
            print("The file: {file_path} do not exists, please check!")
        return False
    
    if os.path.getsize(file_path) == 0:
        if print_log:
            print("The file: {file_path} is empty, please check!")
        return False

    with open(file_path, 'rb') as f:
        f.seek(-1, os.SEEK_END)
        last_byte = f.read(1)

    if last_byte != b'\n':
        if print_log:
            print("The file: {file_path} is not completed, the end charater is not \\n")
        return False

    return True


def check_dir(dir):
    if os.path.exists(dir):
        pass
    else:
        os.makedirs(dir,exist_ok=True)
    

def check_file(file,delete=False):
    if os.path.exists(file) and delete==True:
        os.remove(file)
        return False
    elif os.path.exists(file) and delete==False:
        return True
    else:
        return False
    

######### umi combine related #########
def phred_2_q(phred):
    try:
        phred=int(phred)
        q=1 - (10 ** -(phred/10))
    except ValueError as e:
        q=0
        print(f"wrong phred informat: {phred}")
        print("Error:", e)
    return q


def q_2_phred(q):
    try:
        q=float(q)
        phred=ceil(-log10(1-q) * 10)
    except ValueError as e:
        print(f"wrong q informat: {q}")
        print("Error:", e)
        # sys.exit()
        phred=0
    return phred


def trans(qual_list):
    if qual_list ==[]:
        qual_str="NA"
    else:
        qual_str=",".join(qual_list)
    return qual_str


######### cluster/ind count combine related #########
def combine_alt(series):
    """Combine the alt column"""
    combined = []
    for alt in series:
        if alt != ".":
            alt_list = alt.split(',')
            combined = list(set(combined) | set(alt_list))
    if not bool(combined):
        result = "."
    else:
        result = ",".join(combined)
    return result


def combine_UMI_count(series):
    """Combine UMI count column"""
    combined = [0, 0, 0, 0]
    for umi_count in series:
        if isinstance(umi_count, str) and umi_count!="":  # Check if umi_count is a string
            count = [int(i) for i in umi_count.split(',')]
            combined = [(a + b) for a, b in zip(combined, count)]
    # convert to string
    result = ','.join([str(i) for i in combined])
    return result


def combine_q_columns(series, epsQ=20):
    """Combine quality columns"""
    combined = {}
    for q in series:
        if isinstance(q, str) and q != "":  # Check if q_column is a string
            q_dict = str2dict(q)
            q_filter = {k:v for k,v in q_dict.items() if k>=epsQ}
            combined = dict(Counter(combined) + Counter(q_filter))
    # convert to string
    if not bool(combined):
        result = "NA"
    else:
        result = ','.join(f'{int(key)}:{value}' for key, value in combined.items())
    return result


######### ????? #########

def get_chrom_list_from_file(file):
    '''
    gnomad file list or one gnomad file
    gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz
    gnomad.genomes.v3.1.2.sites.chr2.vcf.bgz
    ...
    gnomad.genomes.v3.1.2.sites.chrX.vcf.bgz
    gnomad.genomes.v3.1.2.sites.chrY.vcf.bgz
    gnomad.genomes.v3.1.sites.chrM.vcf.bgz
    '''
    gnomad_file_list=[];chr_list=[]
    suffix=str(file).split(".")[-1]
    if suffix in ["gz","bgz","vcf"]:
        # this may be a gnomad file, not a list file
        gnomad_file_list=[file]
        chr_list=["none"]
    else:
        # this may be a gnomad list file
        with open(file ,"r")as f1:
            gnomad_file_list=[sline.strip() for sline in f1.readlines()]
        elements_all=[]
        for files in gnomad_file_list:
            elements=re.split("\.|_",str(os.path.basename(files)))
            elements_all+=elements
        chr_list=[k for k, v in dict(collections.Counter(elements_all)).items() if v == 1]

    gnomad_dict={}
    for gnomad_file, chrom in zip(gnomad_file_list, chr_list):
        gnomad_dict[str(chrom)]=gnomad_file
    return gnomad_dict


def get_chrom_list_from_list(in_list):
    '''
    the file list like:
    [file1,file2,file3]
    '''
    file_dict={}
    #print(in_list)
    if len(in_list)==1:
        chrom=re.split("\.|_",os.path.basename(in_list[0]))[0:-1][-1]
        file_dict[chrom]=in_list[0]
    else:
        elements_all=[];chr_list=[]
        for files in in_list:
            elements=re.split("\.|_",str(os.path.basename(files)))
            elements_all+=elements
        chr_list=[k for k, v in dict(collections.Counter(elements_all)).items() if v == 1]

        # print(chr_list)
        for file,chrom in zip(in_list,chr_list):
            file_dict[str(chrom)]=file
    print(file_dict)
    return file_dict

    
def intersect_pos_gnomad(dict1,dict2,out_dir):
    # dict1 is pos file dict, dict2 is annovar file dict
    out_dict={}
    # print(dict1,dict2)
    for key in dict1.keys():
        if key in dict2.keys():
            out_name=os.path.join(out_dir,"short_gnomad_"+key+".vcf")
            # print(out_name)
            # print(check_file(out_name))
            if not check_file(out_name):
                chr_POS=BedTool(dict1[key])
                chr_gnomAD=BedTool(dict2[key])
                # print(key, dict1[key],dict2[key])
                chr_gnomAD.intersect(chr_POS,u=True).saveas(out_name, trackline="#track name='short gnomad file of {chrom}' color=128,0,0".format(chrom=key)) # type: ignore
            out_dict[key]=out_name
        else:
            pass
    # print(out_dict)
    return out_dict


def intersect_pos_gnomad_only_2file(file1,file2,chrom,out_dir):
    # dict1 is pos file dict, dict2 is annovar file dict

    out_name=os.path.join(out_dir,"short_gnomad_"+chrom+".vcf")
    # print(out_name)
    # print(check_file(out_name))
    if not check_file(out_name):
        chr_POS=BedTool(file1)
        chr_gnomAD=BedTool(file2)
        # print(key, dict1[key],dict2[key])
        chr_gnomAD.intersect(chr_POS,u=True).saveas(out_name, trackline="#track name='short gnomad file of {chrom}' color=128,0,0".format(chrom=chrom)) # type: ignore

    return out_name


def split_vcf(vcf_file,split_vcf_dir,work_dir,split_type,split_line=10000):
    vcf_file = os.path.abspath(vcf_file)
    if split_type=="c": # split by chromosome

        command='cd %s;cat %s |cut -f 1 |grep -v "#" |sort -u|while read C; do vcftools  --chr "${C}" --vcf %s  --stdout --recode >split.${C}.vcf; done;' \
                'cd %s'%(split_vcf_dir, vcf_file, vcf_file,work_dir)
    elif split_type=="l": # spliy by line number
        command=f"cd {split_vcf_dir};grep -v '#' {vcf_file}|split -a 4 -d --additional-suffix=.vcf -l {split_line} ; cd {work_dir}"
    else:
        print(f"wrong input for split type {split_type}!")
    
    print(command)

    result=subprocess.run(command,shell=True)
    # print(result.returncode)
    if result.returncode!=0:
        print(f"Something wrong when run the command: {command}")
    else:
        print('finish split vcf file')

    files = os.listdir(split_vcf_dir)
    prefixes = [os.path.basename(f).split(".vcf")[0] + '\n' for f in files if f.endswith(".vcf")]
    vcf_list=open(os.path.join(split_vcf_dir,"vcf_list"),"w")
    vcf_list.writelines(prefixes)


def split_bam_by_chrom(bam_file,out_dir,work_path,run_index=False):
    # check_dir(out_dir)
    # bam_name=os.path.basename(bam_file).split(".")[0]
    bam_name=os.path.basename(bam_file)
    bam_file = os.path.abspath(bam_file)
    # bam_dir=os.path.dirname(bam_name)
    outname_prefix=bam_name+".REF_"
    split_bam_list=[]; split_bam_dict={}
    # run_index=False
    print(bam_file,out_dir)
    if os.listdir(out_dir)==[]:
        run_index=True
    # else:               
    #     for file in os.listdir(out_dir):
    #         if outname_prefix not in file:
    #             print(outname_prefix,file)
    #             run_index=True
    print(run_index)
    if run_index == True:
        #result=subprocess.run("echo yes",shell=True)
        result=subprocess.run("cd %s;rm *;ln -s %s .; samtools index %s ;bamtools split -in %s -reference;cd %s" \
                              "" %  (out_dir, bam_file, bam_name, bam_name,work_path),shell=True)

        if result.returncode!=0:
            print(f'Something wrong when running code:\n bamtools split -in {bam_file} -reference')
            # sys.exit()
        else:
            print("finish split bam file")
    else:
        print("split file already exist")

    for file in os.listdir(out_dir):
        # if outname_prefix in file and file.split(".")[-1]!="bai":
        if ".REF_" in file and file.split(".")[-1]!="bai":
            split_bam_list.append(file)
        if file.split(".")[-1]!="bai" and file+".bai" not in os.listdir(out_dir):
            result=subprocess.run("cd %s; samtools index %s; cd %s" %  (out_dir, file, work_path),shell=True)
            if result.returncode!=0:
                print(f'Something wrong when running code:\n bamtools split -in {bam_file} -reference')
                # sys.exit()
    elements_all=[];chr_list=[]
    for files in split_bam_list:
        # elements=str(os.path.basename(files).replace(bam_name+".REF_","")).replace(".bam","")
        # elements_all.append(elements) 
        prefix=os.path.basename(bam_file).replace(".bam","")+".REF_"
        elements=str(os.path.basename(files).replace(bam_name+".REF_","")).replace(".bam","").replace(prefix,"") 
        elements_all.append(elements) 
    print(elements_all)
    # chr_list=[k for k, v in dict(collections.Counter(elements_all)).items() if v == 1]

    # print(chr_list)
    for file,chrom in zip(split_bam_list,elements_all):
        split_bam_dict[str(chrom)]=os.path.join(out_dir,file)
    print(split_bam_dict)
    return split_bam_dict


def str2dict(q):
    """Convert quality string to dictionary"""
    dict = {int(i.split(':')[0]):int(i.split(':')[1]) for i in q.split(',') if q !="NA"}
    return dict
    

def str2bool(v):
    """ensure boolean input"""
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
    

def handle_posname(pos_name):
    sitem=pos_name.strip().split("_")
    chrom=str(sitem[0])
    pos=int(sitem[1])
    ref=str(sitem[2])
    alt=str(sitem[3])

    return chrom, pos, ref, alt

def handle_p_value_log10(p_value):
    try:
        in_type=type(p_value)
        if in_type==list:
            p_list=[]
            for p_val in p_value:
                if p_val==0 or p_val=="0":
                    p_val=1e-300
                elif p_val=="NA":
                    p_val=1
                p_list.append(log10(float(p_val)))
            return p_list
        elif in_type==str:
            if p_value=="0":
                p_value=1e-300
            elif p_value=="NA":
                p_value=1
            p=log10(float(p_value))
            return p
        
        elif in_type==float or in_type==np.float64:
            if p_value == 0.0:
                p_value=1e-300
            p=log10(float(p_value))
            return p
    except:
        print("p_val:",p_value)
        return "no"


def get_chr_size(fai_file):
    """
    /storage/douyanmeiLab/yangzhirui/Reference/Cellranger/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai
    """
    chr_sizes=dict()
    for line in open(fai_file,"r"):
        line=line.rstrip()
        fields=line.split('\t')
        chr_sizes[fields[0]]=chr_sizes.get(fields[0],fields[1])

    return chr_sizes


def str_to_dict(one_string):
    return_dict=collections.defaultdict(int)
    for items in one_string.strip().split(","):
        if ":" in items:
            return_dict[items.split(":")[0]]+=int(items.split(":")[1])

    return return_dict


def do_wilicox_sum_test(input_1,input_2,method="two-sided",type="dict"):
    from scipy import stats
    if type=="dict":
        input_list1=[]
        input_list2=[]
        for i in input_1.keys():
            input_list1.extend(i[0]*int(i[1]))

        for m in input_2.keys():
            input_list2.extend(m[0]*int(m[1]))    
    elif type=="list":
        input_list1=input_1
        input_list2=input_2
    else:
        print(f"wrong input in do_wilicox_sum_test: {type}")
        return "no","no"
    
    z_statistic, p_value=stats.ranksums(input_list1,input_list2,alternative=method)
    z_statistic=float(z_statistic);p_value=float(p_value)
    return z_statistic,p_value


def get_geno_from_alignment(pos,item):
    def handle_cigar(ciagr_symbol):
        '''
        ## handel cigar
        # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
        # '76M1D33M139241N11M'
        # the 1st is symbol; and the 2nd is count
        # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
        '''
        seq_length_before = 0
        pos_length_before = 0

        seq_cut_start = None; seq_cut_end = None
        pos_cut=[]
        for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
            symbol = cigars[0]
            count = cigars[1]
            if symbol in [5,6,7,8]:
                # an api for handeling mapping issues "HP=X"
                print(ciagr_symbol)  ## LOG
            elif symbol in [0, 1, 4]:
                # measure the seq length 
                seq_length_before += count
                if symbol == 0:
                    pos_length_before += count
                elif symbol == 4:
                    # whether "S" is in this read
                    if i == 1:
                        # whether the "S" is in the head or tail
                        seq_cut_start = seq_length_before
                    elif i ==len(ciagr_symbol):
                        seq_cut_end = seq_length_before
                    else:
                        print(ciagr_symbol) ## LOG
                elif symbol == 1:
                    # whether the "I" is in the cigar
                    pos_cut.append((pos_length_before,count))
            else:
                pass
        seq_cut = (seq_cut_start, seq_cut_end)
        return seq_cut, pos_cut

    def handle_seq(seq, seq_cut):
        # only support one time for cut
        cut_seq=seq[seq_cut[0]:seq_cut[1]]
        return cut_seq

    def handle_pos(pos_matrix,pos_cut):
        if len(pos_cut) == 0:
            cut_pos_matrix = pos_matrix
        elif len(pos_cut) == 1:
            times = pos_cut[0][1]
            pos = pos_cut[0][0]
            cut_pos_matrix = pos_matrix[0:pos] + [""] * times + pos_matrix[pos:]
        else:
            start = 0
            cut_pos_matrix = []
            for item in range(len(pos_cut)):
                pos = pos_cut[item][0]; times = pos_cut[item][1]
                cut_pos_matrix = cut_pos_matrix + pos_matrix[start:pos] + [""] * times
                start=pos
                #print(cut_pos_matrix)
            last_pos = pos_cut[-1][0]
            cut_pos_matrix = cut_pos_matrix + pos_matrix[last_pos:]
        return cut_pos_matrix

    def handle_quality_matrix(mutation_in_cutseq_index,seq,cut_seq):
        if len(cut_seq[mutation_in_cutseq_index:]) >= len(cut_seq[:mutation_in_cutseq_index]):
            query_str = cut_seq[mutation_in_cutseq_index:]
            raw_index = seq.index(query_str)
        else:
            query_str = cut_seq[:mutation_in_cutseq_index]
            raw_index = seq.index(query_str) + len(query_str)
        return raw_index

    Name,geno,quality='','',''
    
    try:
        CB=item.get_tag("CB").strip()
        UB=item.get_tag("UB").strip()
        Name=CB+"_"+UB
        pos_index=pos-1

        seq_cut, pos_cut = handle_cigar(item.cigar)
        cut_seq=handle_seq(item.seq, seq_cut)
        cut_pos=handle_pos(item.get_reference_positions(), pos_cut)
        raw_index = handle_quality_matrix(cut_pos.index(pos_index),item.seq,cut_seq)
        quality=item.get_forward_qualities()[raw_index]

        if pos_index in cut_pos:
            # effective_DP += 1
            geno = cut_seq[cut_pos.index(pos_index)]            
    except:
        pass

    return Name,geno,quality


def get_hard_clip_count(ciagr_symbol):
    """
    same as handle_cigar, but only used to grep hard clip
    """
    left_hardclip,right_hardclip=0,0
    for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol ==5 and i==1:
            left_hardclip=count
        elif symbol==5 and i==len(ciagr_symbol):
            right_hardclip=count

    return left_hardclip,right_hardclip


def get_indel_info(indel_info_list,pos_index_in_raw_matrix):
    """
    similar as handle_pos, input is [(10,1)]. 10 is the indel position to sequence start(exclude soft/hard slip seq)
    output:
    the plus or minus of indel will display the direction between indel and mutation. The plus means indel in the left of mutation
    """
    indel_num=0;indel_length=[];indel_distance=[]
    if len(indel_info_list) == 0:
        indel_num=0
    elif len(indel_info_list) >1 and pos_index_in_raw_matrix=="in":
        indel_num=len(indel_info_list)
        indel_distance=[]
    else:
        indel_num=len(indel_info_list)
        for item in indel_info_list:
            indel_length.append(item[1])
            if pos_index_in_raw_matrix!="in":
                indel_distance.append(pos_index_in_raw_matrix-item[0])
            elif pos_index_in_raw_matrix=="in":
                indel_distance.append(0)
    # indel_length="/".join([str(i) for i in indel_length])
    # indel_distance="/".join([str(i) for i in indel_distance])
    # if indel_num!=0:
    #     print(indel_info_list,pos_index_in_raw_matrix,indel_num,indel_length,indel_distance)
    return indel_num,indel_length,indel_distance


def judge_pos_in_indel(ins_info_list,del_info_list,get_reference_positions):
    pos_list=[]

    def get_pos_list(info_list,pos_list):
        for item in info_list:
            index_num=item[0];count_num=item[1]
            pos=get_reference_positions[index_num-1]
            for i in range(count_num):
                pos_list.append(pos+i+1)
        return pos_list

    pos_list=get_pos_list(ins_info_list,pos_list)
    pos_list=get_pos_list(del_info_list,pos_list)

    return pos_list


# def get_deletion_info(ciagr_symbol,pos_index_in_raw_matrix):
#     del_num=0;del_length=[];del_distance=[]
#     seq_length_before=0
#     for cigars, i in zip(ciagr_symbol,range(1,len(ciagr_symbol)+1)):
#         symbol = cigars[0]
#         count = cigars[1]
#         if symbol ==0:
#             # measure  the seq length 
#             seq_length_before += count
#         elif symbol==2:
#             del_num+=1
#             del_length.append(count)
#             del_distance.append(pos_index_in_raw_matrix-seq_length_before)
#         else:
#             pass
    
#     return del_num,del_length,del_distance


def combine_info_from_cigar(cigar_symbol):
    '''
    ## combine: handel cigar;  get_hard_clip_count
    # [(0, 76), (2, 1), (0, 33), (3, 139241), (0, 11)]
    # '76M1D33M139241N11M'
    # the 1st is symbol; and the 2nd is count
    # 0: Match; 1: Insertion; 2: deletion; 3: N; 4: S; 5: H; 6: P; 7: =; 8: X
    '''
    seq_length_before = 0
    pos_length_before = 0

    left_hardclip,right_hardclip=0,0

    seq_cut_start = None; seq_cut_end = None
    pos_cut=[]
    del_info=[]
    for cigars, i in zip(cigar_symbol,range(1,len(cigar_symbol)+1)):
        symbol = cigars[0]
        count = cigars[1]
        if symbol in [6,7,8]:
            # an api for handeling mapping issues "HP=X"
            print(cigar_symbol)  ## LOG
        elif symbol in [0, 1, 2, 4,5]:
            # measure the seq length 
            seq_length_before += count
            if symbol == 0:
                pos_length_before += count
            elif symbol == 4:
                # whether "S" is in this read
                if i == 1:
                    # whether the "S" is in the head or tail
                    seq_cut_start = seq_length_before
                elif i ==len(cigar_symbol):
                    seq_cut_end = seq_length_before
                else:
                    print(cigar_symbol) ## LOG
            elif symbol == 1:
                # whether the "I" is in the cigar
                pos_cut.append((pos_length_before,count))
            elif symbol == 2:
                del_info.append((pos_length_before,count))
            elif symbol==5:
                print(cigar_symbol)  ## LOG
                if i==1:
                    left_hardclip=count
                elif i==len(cigar_symbol):
                    right_hardclip=count
        else:
            pass
    seq_soft_cut = (seq_cut_start, seq_cut_end)
    seq_hard_clip = (left_hardclip,right_hardclip)

    return seq_soft_cut, pos_cut, del_info, seq_hard_clip


def list2min(value):
    """find the minimum value in the comma-separated string"""
    # check if the input is a float number
    if isinstance(value, float):
        return value
    else:
        return min(float(num) for num in value.split(','))
    
    
def list2frequent(value):
    """find the most frequent element in the comma-separated string"""
    return collections.Counter(value.split(',')).most_common(1)[0][0]


## FUNCTION: adapted from bcftools
def calc_vdb(pos_list, npos,readlen=120):
    from scipy.special import erfc
    """
    # Note well: the parameters were obtained by fitting to simulated data of
    # 100bp reads. This assumes rescaling to 100bp in bcf_call_glfgen().
    
    """
    # readlen = 100
    assert npos == readlen, "npos must be equal to readlen"

    nparam = 15
    param = np.array([
        [3, 0.079, 18], [4, 0.09, 19.8], [5, 0.1, 20.5], [6, 0.11, 21.5],
        [7, 0.125, 21.6], [8, 0.135, 22], [9, 0.14, 22.2], [10, 0.153, 22.3],
        [15, 0.19, 22.8], [20, 0.22, 23.2], [30, 0.26, 23.4], [40, 0.29, 23.5],
        [50, 0.35, 23.65], [100, 0.5, 23.7], [200, 0.7, 23.7]
    ])

    dp = 0
    mean_pos = 0
    mean_diff = 0
    for i in range(npos):
        if pos_list[i] == 0:
            continue
        dp += pos_list[i]
        mean_pos += pos_list[i] * i

    if dp < 2:
        return inf  # one or zero reads can be placed anywhere

    mean_pos /= dp
    for i in range(npos):
        if pos_list[i] == 0:
            continue
        mean_diff += pos_list[i] * abs(i - mean_pos)

    mean_diff /= dp

    ipos = int(mean_diff)  # tuned for float-to-int implicit conversion
    if dp == 2:
        return (2 * readlen - 2 * (ipos + 1) - 1) * (ipos + 1) / (readlen - 1) / (readlen * 0.5)

    if dp >= 200:
        i = nparam  # shortcut for big depths
    else:
        for i in range(nparam):
            if param[i][0] >= dp:
                break

    if i == nparam:
        # the depth is too high, go with 200x
        pscale = param[nparam - 1][1]
        pshift = param[nparam - 1][2]
    elif i > 0 and param[i][0] != dp:
        # linear interpolation of parameters
        pscale = (param[i - 1][1] + param[i][1]) * 0.5
        pshift = (param[i - 1][2] + param[i][2]) * 0.5
    else:
        pscale = param[i][1]
        pshift = param[i][2]

    return 0.5 * erfc(-(mean_diff - pshift) * pscale)


def calc_mwu_biasZ(a, b, n, left_only, do_Z):
    b_empty = all(x == 0 for x in b)
    e = l = na = nb = 0
    t = 0

    if b_empty:
        for i in reversed(range(n)):
            na += a[i]
            t += (a[i] * a[i] - 1) * a[i]  # adjustment score for ties
    else:
        for i in reversed(range(n)):
            e += a[i] * b[i]
            l += a[i] * nb  # a<b
            na += a[i]
            nb += b[i]
            p = a[i] + b[i]
            t += (p * p - 1) * p  # adjustment score for ties

    if not na or not nb:
        return float('inf')

    U = l + e * 0.5  # Mann-Whitney U score
    m = na * nb / 2.0

    # With ties adjustment
    var2 = (na * nb) / 12.0 * ((na + nb + 1) - t / ((na + nb) * (na + nb - 1)))
    if var2 <= 0:
        return 0 if do_Z else 1

    if do_Z:
        # Z-score
        return (U - m) / np.sqrt(var2)

    # U score, asymmetric for some data types
    if left_only and U > m:
        return float('inf')

    if na >= 8 or nb >= 8:
        # Normal approximation
        return np.exp(-0.5 * (U - m) ** 2 / var2)

    # Placeholder for exact calculation
    # You'll need to implement or replace mann_whitney_1947
    # return mann_whitney_1947(na, nb, U) * np.sqrt(2 * np.pi * var2)
    return "Exact calculation needed"


def logsumexp2(a, b):
    return np.log(1 + np.exp(-np.abs(a - b))) + max(a, b)


def calc_SegBias_for_one_sample(ref_dp,alt_dp):
    """
    adapted from bcftools. Here we removed the condition, which have not only one sample.
    """
    nr = alt_dp  # number of observed non-reference reads
    if not nr:
        return

    avg_dp = (ref_dp + nr)   # average depth
    M = floor(nr / avg_dp + 0.5)  # approximate number of variants sampled in the population
    M = min(max(M, 1), 1)  # Clamp M between 1 and sample number(1)
    f = M / 2.0   # allele frequency
    p = nr   # number of variant reads per sample expected if variant not real (Poisson)
    q = nr / M  # number of variant reads per sample expected if variant is real (Poisson)
    sum = 0
    log2 = log(2.0)


    oi = alt_dp  # observed number of non-ref reads
    if oi:
        tmp = logsumexp2(log(2 * (1 - f)), log(f) + oi * log2 - q)
        tmp += log(f) + oi * log(q / p) - q + p
    else:
        tmp = log(2 * f * (1 - f) * exp(-q) + f**2 * exp(-2 * q) + (1 - f)**2) + p
    sum += tmp

    return sum


def get_intergration_from_identifier_and_one_file(identifier,query_file,index_in_query="0 1 2 3"):
    """
    This function is used to get the intergration result from identifier and one special query file.
    Note: The top 4 columns of query_file must be "chrom pos ref alt"
    
    """
    current_directory = os.path.dirname(os.path.abspath(__file__))
    compare_pl_path=os.path.join(current_directory,"others/compare_files.pl")
    perl_script=compare_pl_path
    identifier_line="\t".join(identifier.strip().split("_"))

    command=f"echo -e \"{identifier_line}\" |perl {perl_script} - {query_file} {index_in_query}|sort -u"
    try:
        result=subprocess.check_output(command,text=True,shell=True)
        return result
    except:
        print(f"Something wrong when run the command: {command}")
        return ""


def closest_command(tmp_bed_file,exon_sort_gff_file):

    command=f"bedtools intersect -a {exon_sort_gff_file} -b {tmp_bed_file}"
    try:
        results=subprocess.check_output(command,text=True,shell=True)

    except:
        print(f"Something wrong when run the command: {command}")
        results=""
 
    return results


def closest_pybedtools(tmp_bed_file,exon_sort_gff_file):
    mutation_Bed=BedTool(tmp_bed_file)
    # mutation_Bed=BedTool(per_line,from_string=True)
    gff_Bed=BedTool(exon_sort_gff_file)
    mutation_closest=mutation_Bed.closest(gff_Bed)
    return mutation_closest

