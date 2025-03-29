import math
import numpy as np
import pandas as pd


def cluster_vaf(cluster_count, ind_geno):
    """
    Calculate the mutant allele frequency for each cluster at the candidate site

    Inputs:
        cluster_count - a dataframe with only one row containing the following information
        #chrom	site	ID	ref	alt	cluster spot_number	consensus_read_count	read_quality_A	read_quality_T	read_quality_C	read_quality_G
        ind_geno - a dataframe with only one row containing the following information
        #chrom  site    ID      germline        mutant  cluster spot_number     consensus_read_count    genotype        p_mosaic        Gi      vaf
    Outputs:
        cluster_vaf - a dataframe with only one row containing the following information
        #chrom	site	ID	germline mutant	cluster spot_number	consensus_read_count  vaf
    """
    # get info from input
    umi_count = cluster_count[7]
    germline = ind_geno[3]
    mutant = ind_geno[4]

    # calculate vaf
    nucleotide_list = ["A", "T", "C", "G"]
    counts = list(map(int, umi_count.split(',')))
    count_dict = dict(zip(nucleotide_list, counts))
    depth = sum(counts)
    vaf = count_dict[mutant] / depth if depth>0 else "NA"
    
    # output
    output_list = cluster_count[0:3] + [germline, mutant] + [str(cluster_count[5])] + cluster_count[6:8] + [vaf]

    return output_list


def calculate_percluster(cluster_count, ind_geno_file):
    results = []
    identifier = cluster_count[0] + str("_") + cluster_count[1]
    if cluster_count[0][0] == "#":
        return results
    
    try:
        # print("run vaf cluster for", identifier)
        ind_geno_list = ind_geno_file.loc[identifier]
        if isinstance(ind_geno_list, pd.DataFrame):
            for i in range(len(ind_geno_list)):
                ind_geno = ind_geno_list.iloc[i]
                output_list = cluster_vaf(cluster_count, ind_geno)
                results.append("\t".join([str(i) for i in output_list]))
        else:
            ind_geno = ind_geno_list
            output_list = cluster_vaf(cluster_count, ind_geno)
            results.append("\t".join([str(i) for i in output_list]))
    except:
        pass
    return results