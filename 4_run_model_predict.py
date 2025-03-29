import argparse
from module.mutation_prediction import mutation_classification
from utils import check_dir, str2bool


def predict_mutation(args):
    """
    Predict somatic mutations with classification methods
    """
    check_dir(args.outdir)
    mutation_classification(args.input, args.outdir, args.outprefix, model_dir=args.model_dir, model_name=args.model_name, random_seed=args.random_seed, \
                            train=args.train, true_sites=args.true_sites, artifact_sites=args.artifact_sites, \
                            thr_altcount=args.thr_altcount, thr_altSpotNum=args.thr_altSpotNum, \
                            subset=args.subset, drop_subset=args.drop_subset, hard_filter=args.hard_filter, phase_refine=args.phase_refine, \
                            artifact_downsample_rate=args.downsample_rate, encoder=args.encoder, save_models=args.save, \
                            plot=args.plot, annotate_mosaic=args.annotate_mosaic, annotate_outlier=args.annotate_outlier, n_features=args.n_features, \
                            smote=args.smote, tune=args.tune, k_neighbors=args.k_neighbors, sampling_strategy=args.sampling_strategy, n_jobs=args.n_jobs, \
                            n_estimators=args.n_estimators, max_depth=args.max_depth, min_samples_split=args.min_samples_split)

if __name__ == '__main__':
    ## parameters
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')
    # predict somatic mutations
    parser.add_argument("--input","-i", required=True, help="the directory of features data input")
    parser.add_argument("--outdir",required=True, help="output dir")
    parser.add_argument("--outprefix", help="Custom prefix for the output files.", required=False, default="sample",type=str)
    parser.add_argument("--model_dir", required=False, default="./models_trained/tumor_skin_model", help="the directory of the trained models (default=./models_trained/tumor_skin_model)")
    parser.add_argument("--model_name", required=False, default="tumor_skin_model", help="the sample name of the trained models default=tumor_skin_model")
    parser.add_argument("--random_seed", required=False, default=100, type=int, help="the random seed (default=100)")
    parser.add_argument("--train", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for whether training the models or not (default=True)")
    parser.add_argument("--true_sites", required=False, default=[], help="the manually checked true somatic mutation list (a text file or a list, default=[])")
    parser.add_argument("--artifact_sites", required=False, default=[], help="the list of artifact sites (a text file or a list, default=[])")
    parser.add_argument("--thr_altcount", required=False, default=5, type=int, help="the threshold of the alternative alleles per site (default=5)")
    parser.add_argument("--thr_altSpotNum", required=False, default=None, type=int, help="the threshold of the number of spots contains the alternative alleles, only use when we want to predict for the spot-specific mutations (default=None)")
    parser.add_argument("--subset", required=False, default=None, help="features used to classificate the true somatic mutations, use all features if 'None' (default=None)")
    parser.add_argument("--drop_subset", required=False, default=None, help="the features do not want to be used to classificate the true somatic mutations (default=None)")
    parser.add_argument("--hard_filter", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable of whether performing the hard filters on the predicted sites (default=True)")
    parser.add_argument("--phase_refine", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for whether use the phasing refinement model (default=True)")
    parser.add_argument("--downsample_rate", required=False, default=0.5, type=float, help="the downsample rate for artifacts when generating the training data set (default=0.5)")
    parser.add_argument("--encoder", required=False, default='label', choices=['label', 'onehot'], help="the encode method for object columns ('label' or 'onehot', default='label')")
    parser.add_argument("--save", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for whether saving the models (default=True)")
    parser.add_argument("--plot","-p", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable for plotting the feature importances and the PCA figures (default=True)")
    parser.add_argument("--annotate_mosaic", required=False, default=True, type=str2bool, choices=[True, False], help="Boolean variable of whether annotating the mosaic sites in the training set PCA scatter plots (default=True)")
    parser.add_argument("--annotate_outlier", required=False, default=False, type=str2bool, choices=[True, False], help="Boolean variable of whether annotating the outlier sites in all PCA scatter plots (default=False)")
    parser.add_argument("--n_features", required=False, default=20, type=int, help="the number of the most-important features from the random forest model used in the PCA projection (default=20)")
    parser.add_argument("--smote", required=False, default=True, type=str2bool, choices=[True, False], help="whether use SMOTE to over-sample the minority class to treat the imbalance class values (default = True)")
    parser.add_argument("--tune", required=False, default='random_search', choices=['Bayesian_opt', 'random_search', 'grid_search', None], \
                                help="whether tuning the hyperparameters used in the random forest (choices=['Bayesian_opt', 'random_search', 'grid_search', None], default='random_search')")
    parser.add_argument("--k_neighbors", required=False, default=4, type=int, help="the nearest neighbors used to define the neighborhood of samples in SMOTE (default=4)")
    parser.add_argument("--sampling_strategy", required=False, default='auto', choices=['minority', 'not minority', 'not majority', 'all', 'auto'], \
                                help="sampling information to resample the data set (choices=['minority', 'not minority', 'not majority', 'all', 'auto'], default='auto')")
    parser.add_argument("--n_jobs", required=False, default=None, choices=[None, -1], help="number of jobs to run in parallel (default=None)")
    parser.add_argument("--n_estimators", required=False, default=100, type=int, help="the number of trees in the forest (default=100)")
    parser.add_argument("--max_depth", required=False, default=None, type=int, help="the maximum depth of the tree (default=None)")
    parser.add_argument("--min_samples_split", required=False, default=2, type=int, help="the minimum number of samples required to split an internal node (default=2)")    

    args = parser.parse_args()
