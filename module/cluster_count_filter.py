import pandas as pd
from scipy.stats import binomtest
from utils import str2dict


def cluster_allele_filter(cluster_count_file, output_file, alpha=0.05, epsAF=0.01):
    """
    Get the filtered consensus read counts and qualities for clusters
    here we keep the alternative alleles if they appears at least in one cluster

    Inputs:
        cluster_count_file - the umi count file at the cluster level after quality filtering
        output_file - the path for the output file
        alpha - the significance level for the binomial test when filtering AF (default=0.05)
        epsAF - the threshold for the alternative allele frequency or called the backgroud error (default=0.01)
                ignore the alternative allele with AF < epsAF by binomial test in the case of multiple alternative alleles
    """
    df = pd.read_csv(cluster_count_file, sep="\t", header=None, \
                 names=['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_num', 'umi_count', 'qA', 'qT', 'qC', 'qG'], \
                 keep_default_na=False, comment = "#")
    if not df.empty:
        # create the filtered columns
        df[['keepA', 'keepT', 'keepC', 'keepG']] = df.apply(lambda row: allele_filter(row, alpha=alpha, epsAF=epsAF), axis=1, result_type='expand')

        # combine the columns by the clusters
        def boolean_sum(series):
            return bool(series.sum())
        allele_cluster = df.groupby(['chr', 'pos']).agg({
            'keepA': boolean_sum,
            'keepT': boolean_sum,
            'keepC': boolean_sum,
            'keepG': boolean_sum
        }).reset_index()
        # rename the keep columns
        allele_cluster = allele_cluster.rename(columns={
            'keepA': 'totA',
            'keepT': 'totT',
            'keepC': 'totC',
            'keepG': 'totG'
        })
        df = df.merge(allele_cluster, on=['chr', 'pos'])

        # fill in the quality counts if keeped otherwise deleted
        df[['qA_final', 'qT_final', 'qC_final', 'qG_final']] = df.apply(quality_choose, axis=1, result_type='expand')
        df = df[['chr', 'pos', 'ID', 'ref', 'alt', 'cluster', 'spot_num', 'umi_count', 'qA_final', 'qT_final', 'qC_final', 'qG_final']]
        # rename the keep columns
        df = df.rename(columns={
            'qA_final': 'qA',
            'qT_final': 'qT',
            'qC_final': 'qC',
            'qG_final': 'qG'
        })
        # calculate the consensus read counts
        df['umi_count'] = df.apply(count_umi, axis=1, result_type='expand')
        # find the remaining alternative alleles
        df['alt'] = df.apply(check_alt, axis=1)

    # write the file
    df.rename(columns={'chr': '#chrom'}, inplace=True)
    df.to_csv(output_file, sep="\t", index=None, na_rep='NA')


def allele_filter(row, alpha=0.05, epsAF=0.01):
    """
    Delete the alternative alleles with low allele frequency by performing a binomial test for each cluster
    (AF < epsAF at alpha significance level)

    Outputs: columns show whether the allele is keeped or not
    """
    nucleotide_list = ["A", "T", "C", "G"]
    # transfrom the string to dict
    q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
    q_filter = dict(zip(nucleotide_list, q_dicts))
    
    # calculate the numbers of the nucleotides
    count_list = [sum(q_filter[x].values()) for x in nucleotide_list]
    depth = sum(count_list)
    count_dict = dict(zip(nucleotide_list, count_list))

    # binomial test for allele frequency
    count_filter = count_dict.copy()
    if depth != 0:
        for allele, count in count_dict.items():
            result = binomtest(count, n=depth, p=epsAF, alternative='greater')
            if result.pvalue >= alpha:
                count_filter[allele] = 0
                q_filter[allele] = {}

    # transform format

    allele_pass = [bool(count_filter[x]) for x in nucleotide_list]
    # true if keep the allele, false if not
    return pd.Series([allele_pass[0], allele_pass[1], allele_pass[2], allele_pass[3]])


def quality_choose(row):
    """
    keep the allele if at least one cluster keeps the allele
    otherwise delete the allele below the background error
    """
    qA = row['qA'] if row['totA'] else "NA"
    qT = row['qT'] if row['totT'] else "NA"
    qC = row['qC'] if row['totC'] else "NA"
    qG = row['qG'] if row['totG'] else "NA"
    return pd.Series([qA, qT, qC, qG])


def count_umi(row):
    """
    Count the consensus reads after filtering
    """
    nucleotide_list = ["A", "T", "C", "G"]
    # transfrom the string to dict
    q_dicts = [str2dict(row[x]) for x in ["qA", "qT", "qC", "qG"]]
    q_all = dict(zip(nucleotide_list, q_dicts))
    # calculate the numbers of the nucleotides
    count_list = [sum(q_all[x].values()) for x in nucleotide_list]

    # transform format
    count = ','.join([str(i) for i in count_list])
    # true if keep the allele, false if not
    return count


def check_alt(row):
    """
    Check the remaining alternative alleles after filtering
    """
    nucleotide_list = ["A", "T", "C", "G"]
    # get the number of consensus reads for each allele
    counts = list(map(int, row['umi_count'].split(',')))
    # reference allele
    ref = row['ref']
    # create a list of (allele, count) pairs excluding the ref allele and counts of zero
    alt_counts = [(nucleotide_list[i], counts[i]) for i in range(4) if counts[i] > 0 and nucleotide_list[i] != ref]
    # sort the list by count in descending order
    alt_counts_sorted = sorted(alt_counts, key=lambda x: x[1], reverse=True)
    # join the alleles to create the alt column
    alt = ','.join([allele for allele, _ in alt_counts_sorted]) if len(alt_counts_sorted)>0 else '.'
    return alt
