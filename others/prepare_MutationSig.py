import argparse
import itertools
from Bio import SeqIO
from collections import defaultdict

def get_substitution_type(chrom_dict,chrom,pos):
    '''
    line also could be an identifier
    return_type: 
    '''
    pos=int(pos)
    # print(chrom,pos)
    up_base=chrom_dict[chrom][pos-1-1]
    down_base=chrom_dict[chrom][pos+1-1]

    return up_base,down_base


def trans_base(base):
    try:
        reverse_dict={"A":"T","T":"A","C":"G","G":"C"}
    except:
        return ""
    reverse_base=reverse_dict[base]
    return reverse_base

def get_default_labels(choice: str) -> str:
    if choice not in ["DNA", "RNA", "96", "192"]:
        raise ValueError("Choice must be 'DNA'/'96' or 'RNA'/'192")

    if choice=="DNA" or choice=="96":
        mid_list=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
        value=96
    else:
        mid_list=["C>A", "C>G", "C>T", "T>A", "T>C", "T>G","A>C","A>G","A>T","G>A","G>C","G>T"]
        value=192

    first = ["A", "T", "C", "G"]
    inner_bracket = [[x] * 16 for x in mid_list]
    inner_bracket = [item for sublist in inner_bracket for item in sublist]
    outter_bracket = [x for x in list(itertools.product(first, first))]
    result = [
        outter_bracket[f % 16][0]
        + "["
        + inner_bracket[f]
        + "]"
        + outter_bracket[f % 16][1]
        for f in range(0, value)
    ]
    return result


def fasta_to_dict(fasta_file):
    chromosome_dict = {}
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            chromosome_dict[record.id] = record.seq
    return chromosome_dict


def main():
    '''
    IN:
    chrom pos ref alt gene_strand 
    chr10   5498813 G       T       -
    chr10   5498816 A       G       -
    ## To be noticed, cause each read will have different strand direction. so we would use 
    '''

    mode=args.mode
    out_list=get_default_labels(mode)
    out_dict=defaultdict(int)
    if args.colums:
        ref_index=int(args.colums.split(",")[0])-1
        alt_index=int(args.colums.split(",")[1])-1
        strand_index=int(args.colums.split(",")[2])-1
    else:
        ref_index=2
        alt_index=3
        strand_index=4

    chromosome_dict=fasta_to_dict(args.fasta)

    with open(args.input,"r") as f1:
        for line in f1:
            sline=line.strip().split("\t")
            if sline!=[''] and sline[0][0]!="#":
                chrom=sline[0];pos=int(sline[1])
                ref=sline[ref_index][0];alts=sline[alt_index]
                for alt in alts.split(","):
                    up_base,down_base=get_substitution_type(chromosome_dict,chrom,pos)
                    strand=sline[strand_index]
                    if "N" in [up_base,down_base,ref,alt]:
                        continue
                    
                    if strand not in ["+", "-"] and mode in ["RNA","192",192]:
                        print("your input is wrong! strand is not +/-.")
                        raise TypeError
                    
                    if mode in ["DNA","96",96]:
                        if ref in "CT":
                            stat_name=up_base+"["+ref+">"+alt+"]"+down_base
                        else:
                            stat_name=trans_base(down_base)+"["+trans_base(ref)+">"+trans_base(alt)+"]"+trans_base(up_base)
                    else:
                        if args.reverse and strand in ["-"]:
                            stat_name=trans_base(down_base)+"["+trans_base(ref)+">"+trans_base(alt)+"]"+trans_base(up_base)
                        else:
                            stat_name=up_base+"["+ref+">"+alt+"]"+down_base
                out_dict[stat_name]+=1

    out_file=open(args.output,"w")
    out_file.write(f"MutationType\tCount\n")
    for key in out_list:
        out_file.write(f'{key}\t{out_dict[key]}\n')

    out_file.close()

## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--input",'-i', required=True,help="input file")
parser.add_argument("--colums",'-c', required=False,default="4,5",help="columns ref, alt alleles and strand located")
parser.add_argument("--fasta",'-f', required=True,help="fasta file")
parser.add_argument("--output",'-o',required=True, help="out file")
parser.add_argument("--reverse",action="store_true",required=False, help="Do you want to reverse by strand, and this only work on RNA mode")
parser.add_argument("--mode",'-m',required=True,choices=["DNA","RNA","96","192"],type=str, help="out file")
args = parser.parse_args()
    
if __name__ == '__main__':
    main()
