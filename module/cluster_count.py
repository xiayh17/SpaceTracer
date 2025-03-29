import pandas as pd
from collections import Counter
from utils import combine_UMI_count, combine_alt, combine_q_columns, str2dict


def UMI_count_cluster(umi_count_file, output_file, type="ind", cluster_file=None, epsQ=20):
    """
    Generate the new UMI count file at cluster or individual level
    by combining the UMI count file at the spot level according to the clusters or combine all.
    Also delete the low quality (quality < epsQ) consensus reads.

    Inputs:
        umi_count_file - the file contain the UMI count info at the spot level
        output_file - the path of the output file
        type - "cluster" if combine the info at the cluster level
               "ind" if combine the info at the individual level
               (default="ind")
        cluster_file - the file of the clusters (default="None")
        epsQ - the threshold for the consensus read quality (default=20)
    """

    # read in data: the spot UMI count file
    count = pd.read_csv(umi_count_file, sep="\t", header=None, names=['chr', 'pos', 'ID', 'ref', 'alt', 'spot_barcode', 'umi_count', 'qA', 'qT', 'qC', 'qG'], \
                        keep_default_na=False, comment = "#")
    # add the spot number column
    count.insert(6, "spot_num", 1)
    if type == "cluster":
        # cluster info
        cluster = pd.read_csv(cluster_file, sep="\t", header=None, names=['spot_barcode', 'cluster'], na_values=[])
        cluster['cluster'] = cluster['cluster'].apply(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() 
                                                      else str(x) if pd.notnull(x) else "NA")
        # merge
        df = count.merge(cluster, on=['spot_barcode'])
    else:
        df = count
        df['cluster'] = "bulk"
    
    # combine the columns by the clusters
    grouped = df.groupby(['chr', 'pos', 'ID', 'ref', 'cluster']).agg({
        'alt': combine_alt,
        'spot_num': 'sum',
        'umi_count': combine_UMI_count,
        'qA': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qT': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qC': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qG': lambda series: combine_q_columns(series, epsQ=epsQ),
    }).reset_index()
    columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_num', 'umi_count', 'qA', 'qT', 'qC', 'qG']
    grouped = grouped[columns]
    # write the file
    grouped.rename(columns={'chr': '#chrom'}, inplace=True)
    grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')


def UMI_count_cell(umi_count_file, output_file, type="cell", cluster_file=None, epsQ=20):
    """
    Generate the new UMI count file at cluster or individual level
    by combining the UMI count file at the spot level according to the clusters or combine all.
    Also delete the low quality (quality < epsQ) consensus reads.

    Inputs:
        umi_count_file - the file contain the UMI count info at the spot level
        output_file - the path of the output file
        type - "cell" only for stereo_seq, to combine bin into cell level
               (default="ind")
        cluster_file - the file of the clusters (default="None")
        epsQ - the threshold for the consensus read quality (default=20)
    """

    # read in data: the spot UMI count file
    count = pd.read_csv(umi_count_file, sep="\t", header=None, names=['chr', 'pos', 'ID', 'ref', 'alt', 'spot_barcode', 'umi_count', 'qA', 'qT', 'qC', 'qG'], \
                        keep_default_na=False, comment = "#")
    # add the spot number column
    count.insert(6, "spot_num", 1)
  
    cell_info=pd.read_csv(cluster_file, sep="\t", header=None, names=['spot_barcode','cell'], na_values=[])
    # cell_info['spot_barcode']=cell_info["x"].astype(str)+"_"+cell_info["y"].astype(str)
    df = count.merge(cell_info, on=['spot_barcode'])
    
    # combine the columns by the clusters
    grouped = df.groupby(['chr', 'pos', 'ID', 'ref', 'cell']).agg({
        'alt': combine_alt,
        'umi_count': combine_UMI_count,
        'qA': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qT': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qC': lambda series: combine_q_columns(series, epsQ=epsQ),
        'qG': lambda series: combine_q_columns(series, epsQ=epsQ),
    }).reset_index()
    columns = ['chr', 'pos', 'ID', 'ref', 'alt', 'cell','umi_count', 'qA', 'qT', 'qC', 'qG']
    grouped = grouped[columns]
    # write the file
    grouped.rename(columns={'chr': '#chrom'}, inplace=True)
    grouped.to_csv(output_file, sep="\t", index=None, na_rep='NA')
    


# def combine_q_columns(series, epsQ=20):
#     """Combine quality columns"""
#     combined = {}
#     for q in series:
#         if isinstance(q, str) and q != "":  # Check if q_column is a string
#             q_dict = str2dict(q)
#             q_filter = {k:v for k,v in q_dict.items() if k>=epsQ}
#             combined = dict(Counter(combined) + Counter(q_filter))
#     # convert to string
#     if not bool(combined):
#         result = "NA"
#     else:
#         result = ','.join(f'{int(key)}:{value}' for key, value in combined.items())
#     return result