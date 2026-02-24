import pandas as pd
from scipy.stats import binomtest
import subprocess
from io import StringIO
from pybedtools import BedTool
import argparse
import os


def binomial_test(row):
    h1, h2 = map(int, row['count'].split(','))
    p_value = binomtest(h1, h1 + h2, p=0.5, alternative='two-sided').pvalue
    return p_value



def main():
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
        # "/storage/douyanmeiLab/yangzhirui/Reference/dbSNP/Homo_sapiens_assembly38.dbsnp138.vcf"
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

    # gtex_file="/storage/douyanmeiLab/yangzhirui/Reference/gtexGene.txt"
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


## parameters
parser = argparse.ArgumentParser()
parser.add_argument("--germline",required=True,help="germline file")
parser.add_argument("--count_threshold","-c", required=False,default=50,type=int, help="the minimum count value the germline sites have")
parser.add_argument("--prior_threshold","-v", required=False,default=0.0001,type=float, help="the minimum prior frequency the germline sites have")
parser.add_argument("--p_threshold","-p", required=False,default=0.01,type=float, help="the maximum p value the germline sites are, in binomial test")
parser.add_argument("--outdir","-o",required=True,help="output dir")
parser.add_argument("--outname",required=False,default="candidate_ASE_sites.txt",type=str,help="output suffix")
parser.add_argument("--gtexGene",required=True, help="gene bed file downloaded from uscs goldpath")
parser.add_argument("--dbSNP",required=False,default="", help="dbSNP vcf file")

args = parser.parse_args()

if __name__ == '__main__':
    main()
