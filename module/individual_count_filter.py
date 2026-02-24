import pandas as pd
from collections import Counter
from utils import combine_UMI_count, combine_alt, combine_q_columns, str2dict


def UMI_count_ind_from_cluster(cluster_filter_count_file, output_file, epsQ=20):
    """
    Generate the new UMI count file at  individual level from cluster count file
    by combining the UMI count file at the cluster level.
    So that remains the alternative alleles which pass the allele frequency (AF) binomial test 
    at the cluster level.

    Inputs:
        cluster_filter_count_file - the file contain the UMI count info at the spot level
        output_file - the path of the output file
    """
    # read in data: the spot UMI count file
    df = pd.read_csv(cluster_filter_count_file, sep="\t", header=None, 
                    names=['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_num', 'umi_count', 'qA', 'qT', 'qC', 'qG'], \
                    keep_default_na=False, comment = "#")
    df['cluster'] = df['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                              else str(x) if pd.notnull(x) else "NA")
    
    # combine the columns by the clusters
    grouped = df.groupby(['chr', 'pos', 'ID', 'ref']).agg({
        'alt': combine_alt,
        'spot_num': 'sum',
        'umi_count': combine_UMI_count,
        'qA': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qT': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qC': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qG': lambda series: combine_q_columns(series, epsQ=epsQ),
    }).reset_index()
    grouped['cluster'] = "bulk"
    columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_num', 'umi_count', 'qA', 'qT', 'qC', 'qG']
    grouped = grouped[columns]
    # write the file
    grouped.rename(columns={'chr': '#chrom'}, inplace=True)
    grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')
    
