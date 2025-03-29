import os
import random
import pandas as pd
import numpy as  np
from collections import Counter
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.model_selection import train_test_split, cross_val_score, RandomizedSearchCV, GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
from sklearn.preprocessing import StandardScaler
from imblearn.over_sampling import SMOTE
from hyperopt import hp, fmin, tpe, STATUS_OK, Trials
from scipy.stats import randint
from joblib import dump, load
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from wordcloud import WordCloud
from utils import check_dir, check_file


def mutation_classification(input_file, output_dir, sample_name, model_dir="./", model_name="P6_rep2", random_seed=100, \
                            train=False, true_sites=[], artifact_sites=[], thr_altcount=0, thr_altSpotNum=None, \
                            subset=None, drop_subset=None, hard_filter=True, \
                            phase_refine=True, artifact_downsample_rate=0.5, encoder="label", save_models=True, \
                            plot=True, annotate_mosaic=True, annotate_outlier=False, n_features=20, \
                            smote=True, tune="random_search", k_neighbors=4, sampling_strategy="auto", n_jobs=None, \
                            n_estimators=100, max_depth=None, min_samples_split=2):
    """
    Apply Random Forest and Logistic Regression method to 
    classify somatic mutation from unphased sites

    Inputs:
        input_file - the direction and the name of the input file, which is the features file
        output_dir - the direction saving the outputs
        sample_name - the sample name of the output files
        model_dir - the directory of the saved trained models (default="./)
        model_name - the sample name of the saved trained model used (default="P6_rep2")
        random_seed - the random state used in the whole progress (default=100)
        train - Boolean variable shows whether training the models with input data or use exist models (default=False)
        true_sites - a list of true sites in this data set which are manually checked (default=[])
        artifact_sites - a list of artifact sites in this data set which have passed the allele count threshold, 
                         if not given then use the sites with 'haplotype>3' (default=[])
        thr_altcount - only consider the sites whose alternative allele number is larger than the threshold (default=0)
        thr_altSpotNum - the threshold of the number of spots contains the alternative alleles, only use when we want to 
                         predict for the spot-specific mutations (default=None)
        subset - the features used to classificate the true somatic mutations, use all features if 'None' (default=None)
        drop_subset - the features do not want to be used to classificate the true somatic mutations (default=None)
        hard_filter - Boolean variable of whether performing the hard filters on the predicted sites (default=True)
        phase_refine - whether doing the phase refinement (default=True)
                       True: use an additional random forest model to classificate the haplotype=3 phased sites into 
                             somatic mutation, heterozygous and artifact three groups and save the mutation sites
                       False: use the random forest model for the unphased sites with the features not related to the phasing
                              to classificate the haplotype=3 phased sites into somatic mutation and artifact groups
        artifact_downsample_rate - the downsample rate for artifacts when generating the training data set (default=0.5)
        encoder - the encode method for object columns (default="label")
                  "label": use LabelEncoder
                  "onehot": use OneHotEncoder
        save_models - Boolean variable of whether saving trained models or not (default=True)
        plot - Boolean variable of whether plotting the feature importance figure and the PCA plots (default=True)
        annotate_mosaic - Boolean variable of whether annotating the mosaic sites in the training set PCA scatter plots (default=True)
        annotate_outlier - Boolean variable of whether annotating the outlier sites in all PCA scatter plots (default=False)
        n_features - the number of the most-important features from the random forest model used in the PCA projection (default=20)
        
        # Following are the parameters used in the random forest and logistic regression model
        smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
                to treat the imbalance class values (default=True)
        tune - whether tuning the hyperparameters used in the random forest (n_estimators, max_depth, min_samples_split)
                "Bayesian_opt": Bayesian Optimization
                "random_search": Random Search
                "grid_search": Grid Search
                None: use the given parameters
                (default="random_search)
        k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default=4)
        sampling_strategy - sampling information to resample the data set (default="auto")
                'minority': resample only the minority class;
                'not minority': resample all classes but the minority class;
                'not majority': resample all classes but the majority class;
                'all': resample all classes;
                'auto': equivalent to 'not majority'.
        n_jobs - number of jobs to run in parallel (default=None) 
                None means 1 unless in a joblib.parallel_backend context. 
                -1 means using all processors.         
        n_estimators - the number of trees in the forest (default=100)
        max_depth - the maximum depth of the tree (default=None)
        min_samples_split - the minimum number of samples required to split an internal node (default=2)
    """
    # read in feature info with warning bad rows (less or more columns)
    df = pd.read_csv(input_file, sep="\t", na_values=['no', 'NA', 'NaN'], index_col=False, on_bad_lines='warn')

    # filter the alternative allele count
    df = df[df['alt_dp_consensus'] >= thr_altcount]
    if thr_altSpotNum is not None:
        df = df[df['alt_SpotNum'] <= thr_altSpotNum]
        
        
    # ====================================================================
    # Data Format
    # ====================================================================
    # transform some columns with lists to float numbers
    if 'mappabilityScore' in df.columns:
        df['mappabilityScore'] = df['mappabilityScore'].apply(list2min)
    if 'distanceExon' in df.columns:
        df['distanceExon'] = df['distanceExon'].apply(list2min)
    if 'MI_p_expression' in df.columns:
        df['MI_p_expression'] = df['MI_p_expression'].apply(list2mean)
    if 'MI_s_expression' in df.columns:
        df['MI_s_expression'] = df['MI_s_expression'].apply(list2mean)
    if 'SGB' in df.columns:
        df['SGB'] = pd.to_numeric(df['SGB'], errors='coerce')
    
    # check the identifier column name
    if '#identifier' in df.columns:
        df.rename(columns={'#identifier': 'identifier'}, inplace=True)
        
    # normalize the depth (Counts Per Million (CPM))
    df['dp_cpm'] = (df['dp'] / df['dp'].sum()) * 1e6
    df = df.drop('dp', axis=1)

    # delete the gene function related columns and other non-reasonable features
    not_related_columns = ['gene_type', 'geneStrand', 'major_read_strand', 'GCcontent', \
                           'anno', 'anno_gene', 'DNAMutationType', 'RNAMutationType', \
                           'alt_SpotNum', 'outlier_clusters', 'vaf_cluster_mean', 'vaf_cluster_std_dev', \
                           'af', 'alpha', 'beta', 'rweigh', 'true_similarity', 'major_read_strand', \
                           'refine_hFDR99', 'refine_hFDR49', 'refine_hFDR29', 'refine_hFDR9', 'refine_hFDR1', \
                           'reverse_dp', 'forward_dp', 'alt_reverse_dp', 'alt_forward_dp']
    for col in not_related_columns:
        if col in df.columns:
            df = df.drop(col, axis=1)

    # delte the columns generated for validation but not prediction
    other_validate_columns = ['alt_querypos_num', 'falt', 'fref', 'left_alt_querypos_num_remove_clip', \
                              'left_ref_querypos_num_remove_clip', 'muts_in_cluster_p', 'muts_in_cluster_s', \
                              'ref_querypos_num', 'right_alt_querypos_num_remove_clip', 'right_ref_querypos_num_remove_clip', \
                              'ase', 'chrom', 'editing', 'end', 'sig_pvalue', 'signature_filter', 'start']
    for col in other_validate_columns:
        if col in df.columns:
            df = df.drop(col, axis=1)

    # drop some features which are not distinguishable for classification
    undistinguishable_columns = ['querypos_p', 'querypos_s', 'leftpos_p', 'leftpos_s', 'seqpos_p', 'seqpos_s', \
                                 'ref_hardclip_prop', 'alt_hardclip_prop', 'ref_softclip_prop', 'alt_softclip_prop', \
                                 'spotNum', 'mut_rate_prob', 'mut_rate_likelihood', 'num_outlier_cluster', 'ref_ins_major_prop', \
                                 'dp_cpm', 'dp_consensus', 'alt_dp_consensus', 'VDB', 'RPBZ', 'hardclip_prop_p', 'hardclip_prop_odds', \
                                 'MI_p_expression', 'MI_s_expression', 'outlier_MI_p', 'outlier_MI_s', 'sb_p', 'sb_odds']
    for col in undistinguishable_columns:
        if col in df.columns:
            df = df.drop(col, axis=1)

    # all haplotype-related columns
    haplotype_columns = ['combine_nearest_phase_haplotype', 'no_combine_most_phase_haplotype', \
                         'no_combine_nearest_phase_haplotype', 'combine_most_phase_haplotype']
    # find the most frequent haplotype
    df['haplotype'] = df.apply(lambda row: merge_haplotype_columns(row, haplotype_columns), axis=1)
    # drop the original haplotype columns
    df = df.drop(haplotype_columns, axis=1)


    # delete the rows if the likelihood columns and AF columns contain NA
    df = df.dropna(subset=['mean_AFspot', 'max_AFspot', 'mosaic_likelihood'])

    # interpolate the NAs in the p-value and t-value columns as log(1)=0 and 0
    psval_cols = [col for col in df.columns if col.endswith('_p') or col.endswith('_s') or col.endswith('_odds')]
    for col in psval_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    # interpolate the NAs in the outlier-related columns but not p-values as 0
    outlier_cols = ['outlier_vaf', 'outlier_MI_p', 'outlier_MI_s', 'num_outlier_cluster']
    for col in outlier_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    # interpolate the NAs in the phasing-related columns but not p-values as 0        
    phase_rel_cols = ['combine_nearest_info_mutant_prop', 'no_combine_nearest_info_mutant_prop', \
                      'combine_most_info_mutant_prop', 'no_combine_most_info_mutant_prop', \
                      'combine_nearest_phase_distance', 'no_combine_nearest_phase_distance', \
                      'combine_most_phase_distance', 'no_combine_most_phase_distance', \
                      'combine_nearest_discordant_prop', 'no_combine_nearest_discordant_prop', \
                      'combine_most_discordant_prop', 'no_combine_most_discordant_prop']
    for col in phase_rel_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)

    # interpolate the distance related NAs with inf
    distance_cols = ['distanceExon']
    for col in distance_cols:
        if col in df.columns:
            df[col] = df[col].fillna(np.inf)
    
    # interpolate the other NAs with 0
    other_num_cols = ['alt2_proportion_consensus', 'spotNum', 'MI_expression', \
                      'ref_ins_major_prop', 'mappabilityScore','mutant_rate', 'alt_multi_map_prop', \
                      'ref_multi_map_prop', 'equal_to_previous_bases', 'cause_ploy_alt']
    for col in other_num_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0)
        
    # delete the columns with all massing values
    df = df.dropna(axis=1, how='all')
    # drop the rows if still contains NAs
    df = df.dropna()

    # replace the infinity value to finite one
    very_large_number = 1e10
    df.replace(to_replace=np.inf, value=very_large_number, inplace=True)
    df.replace(to_replace=-np.inf, value=-very_large_number, inplace=True)

    # avoid duplicate rows in the data frame
    num_duplicates = df['identifier'].duplicated().sum()
    if num_duplicates > 0:
        print(f"There exist {num_duplicates} duplicate rows, only keep the first appearance.")
        df = df.drop_duplicates(subset='identifier', keep='first')
    
    # set index
    df = df.set_index("identifier")
    # get the index for artifacts
    if artifact_sites == []:
        artifact_index = df.index[df["haplotype"] == "haplo>3"]
        print("Using the haplotype>3 sites as artifacts.")
    elif isinstance(artifact_sites, str):
        if check_file(artifact_sites):
            with open(artifact_sites, 'r') as file:
                artifact_index = [line.strip() for line in file]
            print('Training models with the artifact sites given.')
        else:
            artifact_index = df.index[df["haplotype"] == "haplo>3"]
            print(f"Can not find file {artifact_sites}, use haplotype>3 sites as artifacts.")
    elif isinstance(artifact_sites, list):
        artifact_index = artifact_sites
        print('Training models with the artifact sites given.')
    # get the index for each haplotype
    het_index = df.index[df["haplotype"] == "haplo=2"]
    unphased_index = df.index[df["haplotype"] == "unphased"]
    phased_index = df.index[df["haplotype"] != "unphased"]
    haplo3_index = df.index[df["haplotype"] == "haplo=3"]
    hyplotype_label = df['haplotype']
    # remove the hyplotype column
    df = df.drop('haplotype', axis=1)

    # for each non-numeric column, apply encoding
    if encoder == "label":
        # use label encoding
        for column in df.select_dtypes(include=['object']).columns:
            print(column)
            le = LabelEncoder()
            df[column] = le.fit_transform(df[column])
    else:
        # use one-hot encoding
        encoder = OneHotEncoder(sparse=False, handle_unknown='ignore')
        encoded = encoder.fit_transform(df.select_dtypes(include=['object']))
        encoded = pd.DataFrame(encoded, columns=encoder.get_feature_names_out())
        # drop original non-numeric columns and concatenate the new DataFrame
        df = pd.concat([df.drop(df.select_dtypes(include=['object']), axis=1), encoded], axis=1)


    # sub-dataframe for het and mutation classification
    # features for classify heterozygous sites
    het_features = ['AFind',
                    'KS_p',
                    'mut_rate_vaf',
                    'mean_AFspot',
                    'MI_p',
                    'mutant_rate',
                    'mismatches_p']
    df_het_features = df[het_features]

    # get subset features for mutation classification
    df_nophase = df.copy()
    for col in phase_rel_cols:
        if col in df_nophase.columns:
            df_nophase = df_nophase.drop(phase_rel_cols, axis=1)
    if subset is not None:
        df_nophase = df_nophase[subset]
    if drop_subset is not None:
        for col in drop_subset:
            if col in df_nophase.columns:
                df_nophase = df_nophase.drop(col, axis=1)


    # ====================================================================
    # Training Set Identify
    # ====================================================================
    if train:
        # check non-empty true mutation sites
        if true_sites == []:
            train = False
            print("No true mutation sites given, use default trained model instead.")
        elif isinstance(true_sites, str):
            if check_file(true_sites):
                with open(true_sites, "r") as file:
                    true_sites = [line.strip() for line in file]
                print('Training models with the true sites given.')
            else:
                train = False
                print(f"Can not find file {true_sites}, use default trained model instead.")
        elif isinstance(true_sites, list):
            print('Training models with the true sites given.')
    else:
        true_sites = []

    # get the training set and the remaning candidate sets
    truesites_set = set(true_sites)
    artifact_set = set(artifact_index)
    het_set = set(het_index)
    haplo3_set = set(haplo3_index)
    phased_set = set(phased_index)
    unphased_set = set(unphased_index)

    # get the phsed true sites set
    phased_true_set = truesites_set & phased_set
    phased_true_index = df.index.intersection(phased_true_set)
    true_sites = df.index.intersection(true_sites)
    # get the artifact sets without duplicate with true sets
    artifact_set = artifact_set - truesites_set
    # randomly downsample the artifact set
    random.seed(random_seed)
    sample_size = int(artifact_downsample_rate* len(artifact_set))
    downsampled_artifact_set = set(random.sample(list(artifact_set), sample_size))
    artifact_index = df.index.intersection(downsampled_artifact_set)
    # artifact_index = df.index.intersection(artifact_set)

    # get the reamining candidate sets
    # candidate_set = set(df.index) - (truesites_set | artifact_set | het_set)
    candidate_phased_set = haplo3_set - (truesites_set | artifact_set | het_set)
    candidate_unphsed_set = unphased_set - (truesites_set | artifact_set | het_set)

    # get the candidate indices
    # candidate_index = df.index.intersection(candidate_set)
    candidate_phased_index = df.index.intersection(candidate_phased_set)
    candidate_unphased_index = df.index.intersection(candidate_unphsed_set)

    # print the number of each type site
    print("=== Number of sites for each type ===")
    print("evaluated phased mosaic:", len(phased_true_index))
    print("evaluated mosaic:", len(true_sites))
    print("evaluated het:", len(het_index))
    print("evaluated artifact:", len(artifact_index))
    print("haplotype=3:", len(haplo3_index))
    print("candidate_phased:", len(candidate_phased_index))
    print("candidate_unphased:", len(candidate_unphased_index))
    print("=====================================")


    # ====================================================================
    # Training Models
    # ====================================================================
    if train:
        # ====================================================================
        # Phase Refinement Model
        # ====================================================================
        # pair each index with its label
        phased_paired_indices = [(index, "mosaic") for index in phased_true_index] + \
                                [(index, "artifact") for index in artifact_index] + \
                                [(index, "het") for index in het_index]
        # shuffle the list to mix indices of different labels
        random.seed(random_seed)
        random.shuffle(phased_paired_indices)
        # generate training data
        X_phased_df = df.loc[[idx for idx, _ in phased_paired_indices]]
        y_phased_list = [label for _, label in phased_paired_indices]
        # convert to matrix
        X_phased = X_phased_df.to_numpy()
        y_phased = np.array(y_phased_list)
        
        # train the phase refinement random forest model
        print("Start training phasing refinement random forest model...")
        rf_phased, _, _, _ = random_forest(X_phased, y_phased, test_size=0, random_state=random_seed, n_labels=3, \
                                           smote=smote, tune=tune, sampling_strategy=sampling_strategy, n_jobs=n_jobs, \
                                           n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
        # save the model
        if save_models:
            model_dir = os.path.join(output_dir, "models")
            check_dir(model_dir)
            rf_phased_file = os.path.join(model_dir, sample_name+"_rf_phased_model.joblib")
            dump(rf_phased, rf_phased_file)
            # save the feature names
            phased_feature_names = X_phased_df.columns.tolist()
            rf_phased_feature_names_file = os.path.join(model_dir, sample_name+"_rf_phased_feature_names.joblib")
            dump(phased_feature_names, rf_phased_feature_names_file)
            print("Save the phase refinement random forest model")
        # plots
        if plot:
            plot_dir = os.path.join(output_dir, "plots")
            check_dir(plot_dir)
            pca_data_dir = os.path.join(output_dir, "pca_data")
            check_dir(pca_data_dir)
            # plot the feature importance
            sorted_features_phased = feature_importamce_plot(rf_phased, X_phased_df, plot_dir, sample_name, \
                                                    title_note=" of Phasing Refinement", file_note="_phase")
            # perform PCA and plot
            pca_phased, top_features_phased, pca_phased_train_data = PCA_train(sorted_features_phased, X_phased_df, y_phased, plot_dir, sample_name, \
                                                title_note=" of Phasing Refinement", file_note="_phase", n_features=n_features, fig_size=6, \
                                                annotate_mosaic=annotate_mosaic, annotate_outlier=annotate_outlier)
            # save pca data
            pca_phased_train_file = os.path.join(pca_data_dir, sample_name+"_pca_phased_train_data.csv")
            pca_phased_train_data.to_csv(pca_phased_train_file, sep='\t', index=True, header=True)
            # save pca model
            if save_models:
                pca_phased_file = os.path.join(model_dir, sample_name+"_pca_phased_model.pkl")
                with open(pca_phased_file, 'wb') as file:
                    pickle.dump(pca_phased, file)
        #  add the labels to the data frame
        X_phased_df['evaluation'] = y_phased

        # ====================================================================
        # Logistic Regression Model For Het
        # ====================================================================
        if len(het_index) > 0:
            # pair each index with its label
            all_paired_indices = [(index, "nothet") for index in true_sites] + \
                                [(index, "nothet") for index in artifact_index] + \
                                [(index, "het") for index in het_index]
            # shuffle the list to mix indices of different labels
            random.seed(random_seed)
            random.shuffle(all_paired_indices)
            # use part features to train the het
            X_all_df = df_het_features.loc[[idx for idx, _ in all_paired_indices]]
            y_all_list = [label for _, label in all_paired_indices]
            # convert to matrix
            X_all = X_all_df.to_numpy()
            y_all = np.array(y_all_list)

            # train the logistic regression model to seperate the het sites
            print("Start training heterozygous classification logistic regression model...")
            lr_model = logistic_regression(X_all, y_all, test_size=0, random_state=random_seed, k_neighbors=3, smote=False, \
                                        sampling_strategy=sampling_strategy, n_jobs=n_jobs)
            # save the model
            if save_models:
                lr_model_file = os.path.join(model_dir, sample_name+"_lr_model.joblib")
                dump(lr_model, lr_model_file)
                print("Save the heterozygous classification logistic regression model")
            #  add the labels to the data frame
            X_all_df['evaluation'] = y_all
            no_het_lr_model = False
        else:
            no_het_lr_model = True

        # ====================================================================
        # Mutaion Classification Random Forest Model
        # ====================================================================
        # pair each index with its label
        nohet_paired_indices = [(index, "mosaic") for index in true_sites] + \
                                [(index, "artifact") for index in artifact_index]
        # shuffle the list to mix indices of different labels
        random.seed(random_seed)
        random.shuffle(nohet_paired_indices)
        # generate training data
        X_nohet_df = df_nophase.loc[[idx for idx, _ in nohet_paired_indices]]
        y_nohet_list = [label for _, label in nohet_paired_indices]
        # convert to matrix
        X_nohet = X_nohet_df.to_numpy()
        y_nohet = np.array(y_nohet_list)

        # train the mutation classification random forest model
        print("Start training somatic mutation classification random forest model...")
        rf, _, _, _ = random_forest(X_nohet, y_nohet, test_size=0, random_state=random_seed, n_labels=2, \
                                    smote=smote, tune=tune, sampling_strategy=sampling_strategy, n_jobs=n_jobs, \
                                    n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split)
        # save the model
        if save_models:
            rf_model_file = os.path.join(model_dir, sample_name+"_rf_model.joblib")
            dump(rf, rf_model_file)
            # save the feature names
            nohet_feature_names = X_nohet_df.columns.tolist()
            rf_feature_names_file = os.path.join(model_dir, sample_name+"_rf_feature_names.joblib")
            dump(nohet_feature_names, rf_feature_names_file)
            print("Save the somatic mutation classification random forest model")
        # plots
        if plot:
            # plot the feature importance
            sorted_features_nohet = feature_importamce_plot(rf, X_nohet_df, plot_dir, sample_name, 
                                    title_note=" of Mutation Classification", file_note="_mutation")
            # perform PCA and plot
            pca_nohet, top_features_nohet, pca_nohet_train_data = PCA_train(sorted_features_nohet, X_nohet_df, y_nohet, plot_dir, sample_name, \
                                            title_note=" of Mutation Classification", file_note="_mutation", n_features=n_features, \
                                            fig_size=6, annotate_mosaic=annotate_mosaic, annotate_outlier=annotate_outlier)
            # save pca data
            pca_nohet_train_file = os.path.join(pca_data_dir, sample_name+"_pca_nohet_train_data.csv")
            pca_nohet_train_data.to_csv(pca_nohet_train_file, sep='\t', index=True, header=True)
            # save pca model
            if save_models:
                pca_nohet_file = os.path.join(model_dir, sample_name+"_pca_nohet_model.pkl")
                with open(pca_nohet_file, 'wb') as file:
                    pickle.dump(pca_nohet, file)
        #  add the labels to the data frame
        X_nohet_df['evaluation'] = y_nohet          

    else:
        # get the model files
        print(f"Load the saved trained models for {model_name}")
        rf_phased_file = os.path.join(model_dir, model_name+"_rf_phased_model.joblib")
        lr_model_file = os.path.join(model_dir, model_name+"_lr_model.joblib")
        rf_model_file = os.path.join(model_dir, model_name+"_rf_model.joblib")
        rf_phased_feature_names_file = os.path.join(model_dir, model_name+"_rf_phased_feature_names.joblib")
        rf_feature_names_file = os.path.join(model_dir, model_name+"_rf_feature_names.joblib")
        # load the trained models
        rf_phased = load(rf_phased_file)
        if os.path.exists(lr_model_file):
            lr_model = load(lr_model_file)
            no_het_lr_model = False
        else:
            no_het_lr_model = True
        rf = load(rf_model_file)
        phased_feature_names = load(rf_phased_feature_names_file)
        nohet_feature_names = load(rf_feature_names_file)
        if plot:
            # get the pca model files
            pca_phased_file = os.path.join(model_dir, model_name+"_pca_phased_model.pkl")
            pca_nohet_file = os.path.join(model_dir, model_name+"_pca_nohet_model.pkl")
            # load pca models
            with open(pca_phased_file, 'rb') as file:
                pca_phased = pickle.load(file)
            with open(pca_nohet_file, 'rb') as file:
                pca_nohet = pickle.load(file)
        

    # ====================================================================
    # Predict
    # ====================================================================
    # hard filter some criteria
    if hard_filter:
        filter_df = df.copy()
        filter_df = filter_df[filter_df['alt_UMI_consistence_prop'] >= 0.8]
        filter_df = filter_df[filter_df['alt_mismatches_mean'] <= 1.5]
        filter_df = filter_df[filter_df['alt_read_number_perUMI_max'] >= 2]
        filter_df = filter_df[filter_df['cause_ploy_alt'] == 0]
        filter_df = filter_df[filter_df['alt_multi_map_prop'] <= 0.2]
        # filter_df = filter_df[filter_df['sf_test_sig'] == 1]
        filter_df = filter_df[filter_df['mutant_rate'] <= 0.5]
        filter_df = filter_df[(filter_df['median_distance_to_end_remove_clip'] >= 0.05) &
                              (filter_df['median_distance_to_end_remove_clip'] <= 0.95)]

        filter_set = set(filter_df.index)
        phased_notpassfilter_set = candidate_phased_set - filter_set
        unphased_notpassfilter_set = candidate_unphsed_set - filter_set
        # phased_notpassfilter_index = df.index.intersection(phased_notpassfilter_set)
        # unphased_notpassfilter_index = df.index.intersection(unphased_notpassfilter_set)

    # predict the phsed haplotype=3 sites
    if len(candidate_phased_index) > 0:
        pred_phase = True
        if phase_refine:
            # using phase refinement
            candidate_phased_df = df.loc[candidate_phased_index]
            candidate_phased = candidate_phased_df.to_numpy()
            candidate_phased_pred = rf_phased.predict(candidate_phased)
            # get the number of each values
            phsed_pred_counts = Counter(candidate_phased_pred)
            print("Phase refinement prediction result:")
            print(phsed_pred_counts)
        else:
            # using random forest model for non-phased sites
            candidate_phased_df = df_nophase.loc[candidate_phased_index]
            candidate_phased = candidate_phased_df.to_numpy()
            candidate_phased_pred = rf.predict(candidate_phased)
            # get the number of each values
            phsed_pred_counts = Counter(candidate_phased_pred)
            print("Phased set somatic mutation prediction result:")
            print(phsed_pred_counts)

        # get predicted true sites
        candidate_phased_df['pred'] = candidate_phased_pred
        phased_pred_true = candidate_phased_df[candidate_phased_df['pred']=="mosaic"]
        phased_pred_true_barcode = phased_pred_true.index
        phased_pred_het = candidate_phased_df[candidate_phased_df['pred']=="het"]
        phased_pred_het_barcode = phased_pred_het.index
        # save the true sites
        results_dir = os.path.join(output_dir, "results")
        check_dir(results_dir)
        phased_pred_file = os.path.join(results_dir, sample_name+"_phased_pred_truesites.txt")
        with open(phased_pred_file, 'w') as file:
            for item in phased_pred_true_barcode:
                file.write(f"{item}\n")

        # get results after hard filtering
        if hard_filter:
            phased_notpassfilter_index = df.index.intersection(phased_notpassfilter_set)
            candidate_phased_df.loc[phased_notpassfilter_index, 'pred'] = "artifact"
            phased_pred_true = candidate_phased_df[candidate_phased_df['pred']=="mosaic"]
            phased_pred_true_barcode = phased_pred_true.index
            print("After hard filter, the number of mosaic sites:", len(phased_pred_true_barcode))
            # save the true sites
            phased_pred_filter_file = os.path.join(results_dir, sample_name+"_phased_pred_truesites_filter.txt")
            with open(phased_pred_filter_file, 'w') as file:
                for item in phased_pred_true_barcode:
                    file.write(f"{item}\n")
    else:
        pred_phase = False


    # predict the het sites
    # candidate_df = df.loc[candidate_index]
    if no_het_lr_model:
        pred_nohet_barcode = candidate_unphased_index
    else:
        candidate_df = df_het_features.loc[candidate_unphased_index]
        candidate = candidate_df.to_numpy()
        candidate_pred = lr_model.predict(candidate)
        # get the number of each values
        pred_counts = Counter(candidate_pred)
        print("Heterozygous classification prediction result:")
        print(pred_counts)
        # get predicted het and nohet sites
        candidate_df['pred'] = candidate_pred
        pred_het = candidate_df[candidate_df['pred']=="het"]
        pred_het_barcode = pred_het.index
        pred_nohet = candidate_df[candidate_df['pred']=="nothet"]
        pred_nohet_barcode = pred_nohet.index
        
    # predict the somatic mutations
    candidate_nohet_df = df_nophase.loc[pred_nohet_barcode]
    candidate_nohet = candidate_nohet_df.to_numpy()
    candidate_nohet_pred = rf.predict(candidate_nohet)
    # get the number of each values
    nohet_pred_counts = Counter(candidate_nohet_pred)
    print("Somatic mutation classification prediction result:")
    print(nohet_pred_counts)
    # get predicted true sites
    candidate_nohet_df['pred'] = candidate_nohet_pred
    nohet_pred_true = candidate_nohet_df[candidate_nohet_df['pred']=="mosaic"]
    nohet_pred_true_barcode = nohet_pred_true.index
    # save predicted true sites in a file
    results_dir = os.path.join(output_dir, "results")
    check_dir(results_dir)
    pred_file = os.path.join(results_dir, sample_name+"_pred_truesites.txt")
    with open(pred_file, 'w') as file:
        for item in nohet_pred_true_barcode:
            file.write(f"{item}\n")

    # get results after hard filtering
    if hard_filter:
        nohet_notpassfilter_index = candidate_nohet_df.index.intersection(unphased_notpassfilter_set)
        candidate_nohet_df.loc[nohet_notpassfilter_index, 'pred'] = "artifact"
        nohet_pred_true = candidate_nohet_df[candidate_nohet_df['pred']=="mosaic"]
        nohet_pred_true_barcode = nohet_pred_true.index
        print("After hard filter, the number of mosaic sites:", len(nohet_pred_true_barcode))
        # save predicted true sites in a file
        pred_filter_file = os.path.join(results_dir, sample_name+"_pred_truesites_filter.txt")
        with open(pred_filter_file, 'w') as file:
            for item in nohet_pred_true_barcode:
                file.write(f"{item}\n")

    # combine the predicted sites
    if pred_phase:
        total_pred_true_barcode = phased_pred_true_barcode.union(nohet_pred_true_barcode)
        # total_pred_true_barcode = phased_pred_true_barcode + nohet_pred_true_barcode
        total_pred_filter_file = os.path.join(results_dir, sample_name+"_total_pred_truesites.txt")
        with open(total_pred_filter_file, 'w') as file:
            for item in total_pred_true_barcode:
                file.write(f"{item}\n")


    # pca scatter plot
    if plot:
        plot_dir = os.path.join(output_dir, "plots")
        check_dir(plot_dir)
        pca_data_dir = os.path.join(output_dir, "pca_data")
        check_dir(pca_data_dir)
        if train:
            # plot for phasable sites
            if pred_phase and len(candidate_phased_index) > 1:
                if phase_refine:
                    pca_phased_pred_data = PCA_pred(X_phased_df, candidate_phased_df, plot_dir, sample_name, 
                                title_note=" of Phasing Refinement", file_note="_phase", 
                                fig_size=6, pca_train=pca_phased, top_features=top_features_phased,
                                annotate_outlier=annotate_outlier)
                else:
                    pca_phased_pred_data = PCA_pred(X_nohet_df, candidate_phased_df, plot_dir, sample_name, 
                                title_note=" of Phasing Refinement", file_note="_phase", 
                                fig_size=6, pca_train=pca_nohet, top_features=top_features_nohet)
            # plot for mutation prediction
            pca_nohet_pred_data = PCA_pred(X_nohet_df, candidate_nohet_df, plot_dir, sample_name, 
                        title_note=" of Mutation Classification", file_note="_mutation", 
                        fig_size=6, pca_train=pca_nohet, top_features=top_features_nohet)
        else:
            # get feature importances
            rf_phased_importances = rf_phased.feature_importances_
            rf_importances = rf.feature_importances_
            # sort the features and their importances
            rf_phsed_sorted_importances, rf_phased_sorted_features = zip(*sorted(zip(rf_phased_importances, phased_feature_names), reverse=True))
            rf_sorted_importances, rf_sorted_features = zip(*sorted(zip(rf_importances, nohet_feature_names), reverse=True))
            # get the top features
            top_features_phased = list(rf_phased_sorted_features[:n_features])
            top_features_nohet = list(rf_sorted_features[:n_features])

            # plot PCA scatter plot
            if pred_phase and len(candidate_phased_index) > 1:
                if phase_refine:            
                    pca_phased_pred_data = PCA_pred(None, candidate_phased_df, plot_dir, sample_name, title_note=" of Phasing Refinement", file_note="_phase", fig_size=6,
                                top_features=top_features_phased)
                else:
                    pca_phased_pred_data = PCA_pred(None, candidate_phased_df, plot_dir, sample_name, title_note=" of Phasing Refinement", file_note="_phase", fig_size=6,
                                top_features=top_features_nohet)
            pca_nohet_pred_data = PCA_pred(None, candidate_nohet_df, plot_dir, sample_name, title_note=" of Mutation Classification", file_note="_mutation", fig_size=6,
                        top_features=top_features_nohet)
        
        # save pca data
        if pred_phase and len(candidate_phased_index) > 1:
            pca_phased_pred_file = os.path.join(pca_data_dir, sample_name+"_pca_phased_pred_data.csv")
            pca_phased_pred_data.to_csv(pca_phased_pred_file, sep='\t', index=True, header=True)
        pca_nohet_pred_file = os.path.join(pca_data_dir, sample_name+"_pca_nohet_pred_data.csv")
        pca_nohet_pred_data.to_csv(pca_nohet_pred_file, sep='\t', index=True, header=True)



    
# ====================================================================
# Supplementary Functions
# ====================================================================

def list2min(value):
    """find the minimum value in the comma-separated string"""
    # check if the input is a float number
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return min(float(num) if num!="no" else np.nan for num in value.split(','))
    
def list2mean(value):
    """get the mean value in the comma-separated string"""
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return np.mean([float(num) for num in value.split(',')])

def list2frequent(value):
    """find the most frequent element in the comma-separated string"""
    if isinstance(value, float):
        return value
    elif isinstance(value, int):
        return value
    else:
        return Counter(value.split(',')).most_common(1)[0][0]
    


def merge_haplotype_columns(row, haplotype_columns):
    """Merge the multiple haplotype-related columns to one haplotype columns"""
    # extract the haplotype-related values
    values = row[haplotype_columns].dropna().unique()  # remove NaN values and get unique values
    NA_name = "unphased"  #or np.nan
    
    # Rule 1: If all values are the same, then take the same value
    if len(values) == 1:
        return values[0]
    
    # Rule 2: If 'haplo=3' is one of the values, then keep 'haplo=3'
    if 'haplo=3' in values:
        return 'haplo=3'
    
    # Rule 3: Take the most frequent value among the specific columns
    if len(values) > 0:
        most_frequent = row[haplotype_columns].dropna().mode()
        return most_frequent.iloc[0] if not most_frequent.empty else NA_name
    else:
        return NA_name
    


def random_forest(X, y, test_size=0.3, random_state=100, smote=True, tune="random_search", \
                  n_labels = 2, k_neighbors=4, sampling_strategy="auto", n_jobs=None, \
                  n_estimators=100, max_depth=None, min_samples_split=2):
    """
    Generate a random forest model for the classification of each type for the candidate mutation sets

    Inputs:
        X, y - the input data and labels (in the numpy format)
        test_size - the test set used when splitting the training and test set (default = 0.3)
        random_state - the random state for the train-test splitting and random forest model (default = 100)
        smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
                to treat the imbalance class values (default = True)
        tune - whether tuning the hyperparameters used in the random forest (n_estimators, max_depth, min_samples_split)
                "Bayesian_opt": Bayesian Optimization
                "random_search": Random Search
                "grid_search": Grid Search
                (default = "random_search)
        n_labels - the number of labels (choice = {2, 3}, default = 2)
        k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default = 4)
        sampling_strategy - sampling information to resample the data set (default = "auto")
                'minority': resample only the minority class;
                'not minority': resample all classes but the minority class;
                'not majority': resample all classes but the majority class;
                'all': resample all classes;
                'auto': equivalent to 'not majority'.
        n_job - number of jobs to run in parallel (default = None) 
                None means 1 unless in a joblib.parallel_backend context. 
                -1 means using all processors.         
        n_estimators - the number of trees in the forest (default = 100)
        max_depth - the maximum depth of the tree (default = None)
        min_samples_split - the minimum number of samples required to split an internal node (default = 2)
    """
    # split dataset into training set and test set
    if test_size != 0:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)
    else:
        X_train = X
        y_train = y
    
    # ======================================================================
    # Over-sampling for the imbalance target values
    # ======================================================================
    # use SMOTE to over-sampling
    if smote:
        smote = SMOTE(random_state=random_state, sampling_strategy=sampling_strategy, k_neighbors=k_neighbors, n_jobs=n_jobs)
        X_train, y_train = smote.fit_resample(X_train, y_train)
        # get the number of each values
        print("Training set count:", dict(Counter(y_train)))


    # ======================================================================
    # Tune the hyperparameters
    # ======================================================================
    if tune == "Bayesian_opt":
        # define the search space of hyperparameters
        space = {
            'n_estimators': hp.choice('n_estimators', range(100, 500)),
            'max_depth': hp.choice('max_depth', [None, 10, 20, 30, 40]),
            'min_samples_split': hp.choice('min_samples_split', range(2, 11))
        }

        # objective function to minimize
        def objective(space):
            model = RandomForestClassifier(n_estimators=space['n_estimators'], max_depth=space['max_depth'], \
                                           min_samples_split=space['min_samples_split'])
            accuracy = cross_val_score(model, X_train, y_train, cv=5).mean()
            return {'loss': -accuracy, 'status': STATUS_OK}

        # run the algorithm
        trials = Trials()
        best = fmin(fn=objective, space=space, algo=tpe.suggest, max_evals=100, trials=trials)
        # print optimized hyperparameters
        print("Optimized hyperparameters:", best)

        # get the optimized parameter values respectively
        n_estimators_opt = range(100, 500)[best['n_estimators']]
        max_depth_opt = [None, 10, 20, 30, 40][best['max_depth']]
        min_samples_split_opt = range(2, 11)[best['min_samples_split']]
        # get the random forest model
        rf = RandomForestClassifier(n_estimators=n_estimators_opt, 
                                    max_depth=max_depth_opt, 
                                    min_samples_split=min_samples_split_opt,
                                    bootstrap=True, oob_score=True, 
                                    random_state=random_state)
    
    elif tune == "random_search":
        # define the parameter distribution
        param_dist = {
            'n_estimators': randint(100, 500),
            'max_depth': [None, 10, 20, 30, 40],
            'min_samples_split': randint(2, 11)
        }
        # initialize the classifier
        rf = RandomForestClassifier(bootstrap=True, oob_score=True, random_state=random_state)
        # initialize the Random Search model
        random_search = RandomizedSearchCV(estimator=rf, param_distributions=param_dist, n_iter=100, cv=5, verbose=0, \
                                           random_state=random_state, n_jobs=n_jobs)

        # fit the Random Search to the data
        random_search.fit(X_train, y_train)
        # get and print the optimized parameters
        best_params = random_search.best_params_
        print("Optimized hyperparameters:", best_params)
        # get the random forest model
        rf = RandomForestClassifier(**best_params, bootstrap=True, oob_score=True, random_state=random_state)
    
    elif tune == "grid_search":
        # define the parameter grid
        param_grid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [None, 10, 20, 30],
            'min_samples_split': [2, 5, 10]
        }
        # initialize the classifier
        rf = RandomForestClassifier(bootstrap=True, oob_score=True, random_state=random_state)
        # initialize the Grid Search model
        grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5, verbose=0, n_jobs=n_jobs)
        # fit the Grid Search to the data
        grid_search.fit(X_train, y_train)
        # get and print the optimized parameters
        best_params = grid_search.best_params_
        print("Optimized hyperparameters:", best_params)
        # get the random forest model
        rf = RandomForestClassifier(**best_params, bootstrap=True, oob_score=True, random_state=random_state)
    
    else:
        # create a Gaussian classifier with the default hyperparameters
        rf = RandomForestClassifier(n_estimators=n_estimators, max_depth=max_depth, min_samples_split=min_samples_split, \
                                    bootstrap=True, oob_score=True, random_state=random_state)
    

    # ======================================================================
    # Train the model
    # ======================================================================
    # train
    rf.fit(X_train, y_train)
    # print the average out-of-bag accuracy
    print("OOB Score:", rf.oob_score_)

    if test_size != 0:
        # predict for the test set
        y_pred = rf.predict(X_test)
        # predict for the test true mutation set
        true_indices = np.where(y_test == "mosaic")[0]
        X_test_true = X_test[true_indices, :]
        y_test_true = y_test[true_indices]
        y_pred_true = rf.predict(X_test_true)
        # generate confusion matrix
        if n_labels == 3:
            cm = confusion_matrix(y_test, y_pred, labels=['mosaic', 'het', 'artifact'])
        elif n_labels == 2:
            cm = confusion_matrix(y_test, y_pred, labels=['mosaic', 'artifact'])
        else:
            cm = confusion_matrix(y_test, y_pred)
        print(cm)
        # calculate the sensitivity
        test_sensitivity = cm[0, 0] / cm[0, :].sum()

        # print the test set accuracy
        test_accuracy = accuracy_score(y_test, y_pred)
        test_true_accuracy = accuracy_score(y_test_true, y_pred_true)
        print("Test Accuracy:", test_accuracy)
        print("Test True Mutation Accuracy:", test_true_accuracy)
        print(f"Test Sensitivity: {test_sensitivity}")
        
    else:
        test_accuracy, test_true_accuracy, test_sensitivity = None, None, None

    return rf, test_accuracy, test_true_accuracy, test_sensitivity



def logistic_regression(X, y, test_size=0.3, random_state=100, max_iter=1000, class_weight='balanced', \
                        smote=True, k_neighbors=4, sampling_strategy="auto", n_jobs=None):
    """
    Use Logistic regression method to seperate heterozygous sites from candidate sites

    Inputs:
        X, y - the input data and labels (in the numpy format)
        test_size - the test set used when splitting the training and test set (default = 0.3)
        random_state - the random state for the train-test splitting and random forest model (default = 100)
        max_iter - the maximum iteration number (default = 1000)
        class_weight - weights associated with classes (default = "balanced")
        smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
                to treat the imbalance class values (default = True)
        k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default = 4)
        sampling_strategy - sampling information to resample the data set (default = "auto")
                'minority': resample only the minority class;
                'not minority': resample all classes but the minority class;
                'not majority': resample all classes but the majority class;
                'all': resample all classes;
                'auto': equivalent to 'not majority'.
        n_job - number of jobs to run in parallel (default = None) 
                None means 1 unless in a joblib.parallel_backend context. 
                -1 means using all processors.
    """
    # split dataset into training set and test set
    if test_size != 0:
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)
    else:
        X_train = X
        y_train = y

    # ======================================================================
    # Over-sampling for the imbalance target values
    # ======================================================================
    # use SMOTE to over-sampling
    if smote:
        smote = SMOTE(random_state=random_state, sampling_strategy=sampling_strategy, k_neighbors=k_neighbors, n_jobs=n_jobs)
        X_train, y_train = smote.fit_resample(X_train, y_train)
        # get the number of each values
        print("Training set count:", dict(Counter(y_train)))

    # normalize the training set
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    if test_size != 0:
        X_test = scaler.transform(X_test)
    
    # ======================================================================
    # Train the model
    # ======================================================================
    # initialize the Logistic Regression model
    lr_model = LogisticRegression(random_state=random_state, max_iter=max_iter, class_weight=class_weight)
    # fit model with the training data
    lr_model.fit(X_train, y_train)

    if test_size != 0:
        # test accuracy
        y_pred = lr_model.predict(X_test)
        test_accuracy = accuracy_score(y_test, y_pred)

        # print result
        print(confusion_matrix(y_test, y_pred))
        print(classification_report(y_test, y_pred))
        print(f"Test Accuracy: {test_accuracy}")

    return lr_model



# from sklearn.model_selection import train_test_split
# from sklearn.linear_model import LogisticRegression
# from sklearn.metrics import classification_report, confusion_matrix, accuracy_score
# from imblearn.over_sampling import SMOTE
# from sklearn.preprocessing import StandardScaler

# def multinomial_logistic_regression(X, y, test_size=0.3, random_state=100, max_iter=1000, class_weight='balanced', \
#                                     smote=True, k_neighbors=4, sampling_strategy="auto", n_jobs=None):
#     """
#     Use Multinomial Logistic regression method to refine the phasing result

#     Inputs:
#         X, y - the input data and labels (in the numpy format)
#         test_size - the test set used when splitting the training and test set (default = 0.3)
#         random_state - the random state for the train-test splitting and random forest model (default = 100)
#         max_iter - the maximum iteration number (default = 1000)
#         class_weight - weights associated with classes (default = "balanced")
#         smote - whether use SMOTE (Synthetic Minority Over-sampling Technique) to over-sample the minority class 
#                 to treat the imbalance class values (default = True)
#         k_neighbors - the nearest neighbors used to define the neighborhood of samples in SMOTE (default = 4)
#         sampling_strategy - sampling information to resample the data set (default = "auto")
#                 'minority': resample only the minority class;
#                 'not minority': resample all classes but the minority class;
#                 'not majority': resample all classes but the majority class;
#                 'all': resample all classes;
#                 'auto': equivalent to 'not majority'.
#         n_job - number of jobs to run in parallel (default = None) 
#                 None means 1 unless in a joblib.parallel_backend context. 
#                 -1 means using all processors.
#     """
#     # split dataset into training set and test set
#     if test_size != 0:
#         X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, stratify=y, random_state=random_state)
#     else:
#         X_train = X
#         y_train = y


#     # ======================================================================
#     # Over-sampling for the imbalance target values
#     # ======================================================================
#     # use SMOTE to over-sampling
#     if smote:
#         smote = SMOTE(random_state=random_state, sampling_strategy=sampling_strategy, k_neighbors=k_neighbors, n_jobs=n_jobs)
#         X_train, y_train = smote.fit_resample(X_train, y_train)
#         # get the number of each values
#         print("Training set count:", dict(Counter(y_train)))

#     # normalize the training set
#     scaler = StandardScaler()
#     X_train = scaler.fit_transform(X_train)
#     if test_size != 0:
#         X_test = scaler.transform(X_test)
    
#     # ======================================================================
#     # Train the model
#     # ======================================================================
#     # initialize the Multinomial Logistic Regression model
#     mlr_model = LogisticRegression(random_state=random_state, max_iter=max_iter, class_weight=class_weight, \
#                                    multi_class='multinomial', solver='lbfgs')
#     # fit model with the training data
#     mlr_model.fit(X_train, y_train)

#     if test_size != 0:
#         # test accuracy
#         y_pred = mlr_model.predict(X_test)
#         test_accuracy = accuracy_score(y_test, y_pred)

#         # print result
#         print(confusion_matrix(y_test, y_pred))
#         print(classification_report(y_test, y_pred))
#         print(f"Test Accuracy: {test_accuracy}")

#     return mlr_model



def feature_importamce_plot(model, df, output_dir, sample_name, title_note="", file_note="", fig_size=10):
    """Plot the feature importance learnt from random forest or other classification models"""
    # get feature importances
    importances = model.feature_importances_
    feature_names = df.columns
    # sort the features and their importances
    sorted_importances, sorted_features = zip(*sorted(zip(importances, feature_names), reverse=True))
    # plot feature importance bar plot
    plt.figure(figsize=(fig_size+2, fig_size))
    plt.title('Feature Importances' + title_note)
    plt.barh(range(len(sorted_importances)), sorted_importances, align='center')
    plt.yticks(range(len(sorted_importances)), sorted_features)
    plt.xlabel('Relative Importance')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    fig_file = os.path.join(output_dir, sample_name+file_note+"_feature_importance.png")
    plt.savefig(fig_file)
    plt.close()

    # create a dict for feature importance
    feature_importance_dict = dict(zip(feature_names, importances))
    # generate the word cloud
    wordcloud = WordCloud(width=800, height=400, background_color='white').generate_from_frequencies(feature_importance_dict)
    # plot the word cloud
    plt.figure(figsize=(fig_size*2, fig_size))
    plt.imshow(wordcloud, interpolation='bilinear')
    plt.axis('off')
    plt.title('Feature Importance Word Cloud' + title_note)
    feature_wordcloud_file = os.path.join(output_dir, sample_name+file_note+"_features_wordcloud.png")
    plt.savefig(feature_wordcloud_file)
    plt.close()
    return sorted_features


def pca_circle_plot(pca, feature_names, output_dir, sample_name, title_note="", file_note="", lim=1, figsize=7):
    """Plot features in the first two principal components as a circle plot"""
    # get components
    loadings = pca.components_.T
    fig, ax = plt.subplots(figsize=(figsize, figsize))
    # color pallete
    arrow_color = '#4DBAD6'
    text_color = '#E44A33'
    # plot the features
    for i, feature in enumerate(feature_names):
        plt.arrow(0, 0, loadings[i, 0], loadings[i, 1], 
                  color=arrow_color, alpha=0.5)
        plt.text(loadings[i, 0]*1.15, loadings[i, 1]*1.15, 
                 feature, color=text_color, ha='center', va='center')
    # set plot limits
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.axhline(0, color='gray', linestyle='--')
    plt.axvline(0, color='gray', linestyle='--')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA Componants Plot' + title_note)
    fig_file = os.path.join(output_dir, sample_name+file_note+"_features_circle.png")
    plt.savefig(fig_file)
    plt.close()


def PCA_train(sorted_features, X_df, y, output_dir, sample_name, title_note="", file_note="", 
                 n_features=10, fig_size=6, annotate_mosaic=True, annotate_outlier=False, thr_outlier=99.9):
    """Perform PCA with important features"""
    # get the important feature data frame
    features_list = list(sorted_features)
    if "context" in features_list:
        features_list.remove("context")
        sorted_features = tuple(features_list)
    top_features = list(sorted_features[:n_features])
    X_df_features = X_df[top_features]

    # modify the extreme p-values to a proper number
    p_value_cols = [col for col in X_df_features.columns if col.endswith('_p')]
    X_df_features = X_df_features.copy()  # make an explicit copy of the DataFrame
    X_df_features.loc[:, p_value_cols] = X_df_features[p_value_cols].applymap(lambda x: -7 if x < -7 else x)

    # standardize the features
    X_std = StandardScaler().fit_transform(X_df_features)

    # perform PCA
    pca = PCA(n_components=2)  # we choose 2 components for visualization purposes
    X_pca = pca.fit_transform(X_std)
    # plot features
    pca_circle_plot(pca, feature_names=X_df_features.columns, output_dir=output_dir, sample_name=sample_name, 
                    title_note=title_note, file_note=file_note, lim=0.85, figsize=9)

    # create a data frame for the output
    pca_df = pd.DataFrame(data = X_pca, columns = ['PC1', 'PC2'], index=X_df_features.index)
    pca_df['evaluation'] = y
    # get the explained variance ratio
    explained_variance = pca.explained_variance_ratio_
    # calculate the distance from the origin (or centroid)
    distances = np.linalg.norm(X_pca, axis=1)
    # define a threshold (percentile) for outliers
    threshold = np.percentile(distances, thr_outlier)
    # identify outliers
    outliers = distances > threshold
    pca_df['outlier'] = outliers

    # PCA plot
    color_palette = {
        'mosaic': '#fe4a49',
        'artifact': '#009fb7',
        'het': '#FF8C00'
    }
    # if n_labels==2:
    #     color_palette = ['#fe4a49', '#009fb7']
    # elif n_labels==3:
    #     color_palette = ['#fe4a49', '#009fb7', '#ffd166']
    plt.figure(figsize=(fig_size+2, fig_size))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='evaluation', palette=color_palette, alpha=0.9, s=20)
    # add site info for true mutations
    for i, row in pca_df.iterrows():
        if row['evaluation'] == 'mosaic' and annotate_mosaic: 
            plt.annotate(row.name, 
                         (row['PC1'], row['PC2']),
                         textcoords="offset points",
                         xytext=(0,10),  # offset the label position
                         ha='center',
                         fontsize=9)
        if row['outlier'] and annotate_outlier:
            plt.annotate(row.name, 
                         (row['PC1'], row['PC2']),
                         textcoords="offset points",
                         xytext=(0,10),  # offset the label position
                         ha='center',
                         fontsize=9,
                         color="royalblue")
    plt.title('PCA Scatter Plot'+title_note)
    plt.xlabel(f'PC 1 ({explained_variance[0]*100:.2f}% variance)')
    plt.ylabel(f'PC 2 ({explained_variance[1]*100:.2f}% variance)')
    pca_scatter_plot = os.path.join(output_dir, sample_name+file_note+"_pca_scatter.png")
    plt.savefig(pca_scatter_plot)
    plt.close()

    return pca, top_features, pca_df



def PCA_pred(train_df, pred_df, output_dir, sample_name, title_note="", file_note="", 
                fig_size=6, pca_train=None, top_features=None, annotate_outlier=False, thr_outlier=99.9):
    """Perform PCA with all features and plot for both training and testing data"""
    # combine the data frames
    pred_df.rename(columns={'pred':'label'}, inplace=True)
    pred_df['type'] = 'prediction'
    if train_df is not None:
        train_df.rename(columns={'evaluation':'label'}, inplace=True)
        train_df['type'] = 'evaluation'
        combined_df = pd.concat([train_df, pred_df], ignore_index=False)
    else:
        combined_df = pred_df

    # get the labels and the types
    df = combined_df.drop(columns=['label', 'type'])
    # select the features
    if top_features is not None:
        df = df[top_features]
    # modify the extreme p-values to a proper number
    p_value_cols = [col for col in df.columns if col.endswith('_p')]
    df = df.copy()
    df.loc[:, p_value_cols] = df[p_value_cols].applymap(lambda x: -7 if x < -7 else x)

    # standardize the features
    df_std = StandardScaler().fit_transform(df)
    # perform PCA
    if pca_train is not None:
        pca = pca_train
        df_pca = pca_train.transform(df_std)
    else:
        pca = PCA(n_components=2)  # we choose 2 components for visualization purposes
        df_pca = pca.fit_transform(df_std)
    # add the PCs to the combined dataframe
    combined_df[['PC1', 'PC2']] = df_pca
    # get the explained variance ratio
    explained_variance = pca.explained_variance_ratio_
    # calculate the distance from the origin (or centroid)
    distances = np.linalg.norm(df_pca, axis=1)
    # define a threshold (percentile) for outliers
    threshold = np.percentile(distances, thr_outlier)
    # identify outliers
    outliers = distances > threshold
    combined_df['outlier'] = outliers

    # plot features
    title_note = title_note + " with Prediction"
    file_note = file_note + "_pred"
    pca_circle_plot(pca, feature_names=df.columns, output_dir=output_dir, sample_name=sample_name, 
                    title_note=title_note, file_note=file_note, lim=0.55, figsize=9)
    
    # PCA plot
    color_palette = {
        'mosaic': '#fe4a49',
        'artifact': '#009fb7',
        'het': '#FF8C00'
    }
    # marker_styles = {'evaluation': 'o', 'prediction': '^'}
    # combined_df['marker'] = combined_df['type'].map(marker_styles)
    plt.figure(figsize=(fig_size+2, fig_size))
    if train_df is not None:
        sns.scatterplot(data=combined_df, x='PC1', y='PC2', hue='label', style='type', 
                        palette=color_palette, alpha=0.9, s=20, markers=['o', '^'])
    else:
        sns.scatterplot(data=combined_df, x='PC1', y='PC2', hue='label', style='type', 
                        palette=color_palette, alpha=0.9, s=20, markers='^')
    # add site info for outliers
    if annotate_outlier:
        for i, row in combined_df.iterrows():
            if row['outlier']: 
                plt.annotate(row.name, 
                            (row['PC1'], row['PC2']),
                            textcoords="offset points",
                            xytext=(0,10),  # offset the label position
                            ha='center',
                            fontsize=9,
                            color="royalblue")
    # for marker, subset in combined_df.groupby('marker'):
    #     sns.scatterplot(data=subset, x='PC1', y='PC2', hue='label', palette=color_palette, 
    #                     alpha=0.9, s=20, markers=marker, edgecolor='w')
    plt.title('PCA Scatter Plot'+title_note)
    plt.legend(loc='best')
    plt.xlabel(f'PC 1 ({explained_variance[0]*100:.2f}% variance)')
    plt.ylabel(f'PC 2 ({explained_variance[1]*100:.2f}% variance)')
    pca_scatter_plot = os.path.join(output_dir, sample_name+file_note+"_pca_scatter.png")
    plt.savefig(pca_scatter_plot)
    plt.close()

    return combined_df