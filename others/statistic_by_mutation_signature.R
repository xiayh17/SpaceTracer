##Author: Zhirui Yang
##Date: 2026-02-09
##Details: This is one script to calculate hFDR.
args <- commandArgs(trailingOnly = TRUE)


if ("-h" %in% args || "--help" %in% args) {
  cat("Usage: Rscript script.R feature_all_result outpath error_profile_file\n")
  cat("Arguments:\n")
  cat("sample_id: sample id\n")
  cat("feature_all_result: Path to features dataframe file, make sure your dataframe contains these four colnames: identifier, RNAMutationType, dp, alt_dp\n")
  cat("outpath: Path to false signature file\n")
  cat("error_profile_file: Path to mutation signature profile of errors (Each col represent one Signature) \n")
  q(save = "no")
}

#===========parameters=================
sample_id<-args[1]
candidate_df_file<-args[2]
outpath<-args[3]
error_profile_file<-args[4]
reference_error_profile<-args[5]

if (!file.exists(outpath)) {
  dir.create(outpath)
  print(paste("Path", outpath, "created successfully"))
} else {
  print(paste("Path", outpath, "already exists"))
} 

#===========functions=================
library(pracma)
library(dplyr)
library(extraDistr)
library(tidyr)
library(deconstructSigs)


mutation_order <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G","G>A","G>T","G>C","A>T","A>C","A>G")
extract_mutation_type <- function(rownames) {
  mutation_types <- gsub(".*\\[([CTGA]>[ACGT])\\].*", "\\1", rownames)
  return(mutation_types)
}

cos_sim <- function(x, y) {
  res <- x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res <- as.numeric(res)
  return(res)
}

get.sig.score.3artifacts <- function(mutsigs, true.sig, artifact.sig1,artifact.sig2,artifact.sig3,prportions, eps=0.001) {
  test.spectrum <- as.numeric(mutsigs)
  
  sigs <- cbind(as.numeric(as.matrix(true.sig)), as.numeric(as.matrix(artifact.sig1)),as.numeric(as.matrix(artifact.sig2)),as.numeric(as.matrix(artifact.sig3)))
  weights <- pracma::lsqnonneg(sigs, test.spectrum)$x
  weights <- weights + eps
  
  weight1 <- proportions[proportions$Signature == "Signature1", "Proportion"]
  weight2 <- proportions[proportions$Signature == "Signature2", "Proportion"]
  weight3 <- proportions[proportions$Signature == "Signature3", "Proportion"]

  postp <- log10(artifact.sig1*weights[2]*weight1+artifact.sig2*weights[3]*weight2+artifact.sig3*weights[4]*weight3) - log10(true.sig*weights[1])
  list(postp=postp, weight.true=weights[1], weight.artifact1=weights[2], weight.artifact2=weights[3], weight.artifact3=weights[4])
}


handle_per_site_3artifact_remove_weight <- function(row, artifact.sig, proportions, eps) {
  reference_list <- rownames(artifact.sig)
  print(length(reference_list))
  true.sig <- data.frame("Count" = rep(1/192, 192))
  rownames(true.sig) <- rownames(artifact.sig)
  mut_name <- row$RNAMutationType
  index_position <- which(sapply(reference_list, function(x) x == mut_name))
  result_vector <- rep(0, length(reference_list))
  result_vector[index_position] <- 1
  art_props=normalized_refine_sig[mut_name,]
  
  result_list <- list()
  result <- get.sig.score.3artifacts(result_vector, true.sig,artifact.sig[, "Signature1"],artifact.sig[, "Signature2"],artifact.sig[, "Signature3"] ,proportions, eps)
  
  cosine_similarity1 <- cos_sim(result_vector, artifact.sig[, "Signature1"])
  cosine_similarity2 <- cos_sim(result_vector, artifact.sig[, "Signature2"])
  cosine_similarity3 <- cos_sim(result_vector, artifact.sig[, "Signature3"])
    
  result_list["true_similarity"] <- result$weight.true
  result_list["false_similarity1"] <- result$weight.artifact1
  result_list["false_similarity2"] <- result$weight.artifact2
  result_list["false_similarity3"] <- result$weight.artifact3

  weight1 <- as.numeric(proportions[proportions$Signature == "Signature1", "Proportion"][1])
  weight2 <- as.numeric(proportions[proportions$Signature == "Signature2", "Proportion"][1])
  weight3 <- as.numeric(proportions[proportions$Signature == "Signature3", "Proportion"][1])

  result_list["true_prop"] <- 1/192
  result_list["false_prop1"] <- as.numeric(art_props$Signature1[1])
  result_list["false_prop2"] <- as.numeric(art_props$Signature2[1])
  result_list["false_prop3"] <- as.numeric(art_props$Signature3[1])

  rweigh= result_list[["true_prop"]]/(weight1 * art_props$Signature1[1] + weight2 * result_list[["false_prop2"]] + weight3 * result_list[["false_prop3"]]) 
  result_list["rweigh"] <- as.numeric(rweigh)
  result_list["cosine_sim1"] <- cosine_similarity1
  result_list["cosine_sim2"] <- cosine_similarity2
  result_list["cosine_sim3"] <- cosine_similarity3
  
  return(result_list)
}

# function to get standard mutation signature matrix
get_default_labels <- function(choice) {
  if (!(choice %in% c("DNA", "RNA", "96", "192"))) {
    stop("Choice must be 'DNA'/'96' or 'RNA'/'192'")
  }
  
  if (choice == "DNA" | choice == "96") {
    mid_list <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    value <- 96
  } else {
    mid_list <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "A>C", "A>G", "A>T", "G>A", "G>C", "G>T")
    value <- 192
  }
  
  first <- c("A", "T", "C", "G")
  inner_bracket <- rep(rep(mid_list, each = 16), times = 1)
  outter_bracket <- expand.grid(first, first)
  result <- sapply(1:value, function(f) {
    paste0(outter_bracket[f %% 16 + 1, 1], "[", inner_bracket[f], "]", outter_bracket[f %% 16 + 1, 2])
  })
  
  return(result)
}

# refine input artifacts signature
complement_pairs <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
refine_counts <- function(profile){
  new_profile<-data.frame("Count"=rep(0,nrow(false_sig1 )))
  rownames(new_profile)<-rownames(profile)
  for (i in rownames(profile)) {
    char_up=substr(i,1,1)
    char_before=substr(i,3,3)
    char_after=substr(i,5,5)
    char_down=substr(i,7,7)
    complements=complement_pairs[c(char_down,char_before,char_after,char_up)]
    complement_row=paste0(complements[1],"[",complements[2],">",complements[3],"]",complements[4],sep="")
    if (complement_row %in% rownames(profile) ) {
      if (profile[i, "Count"]> profile[complement_row, "Count"]){
        new_profile[i, "Count"] <- profile[i, "Count"] - profile[complement_row, "Count"]
      } else{
        new_profile[i, "Count"] <- 0
      }
    }
  }
  return(new_profile)
}


########## handel proportion file
own_database=reference_error_profile
ref=read.table(own_database,sep="\t",row.names = 1,header = T)
sample_file=error_profile_file
ref_mutation_types <- extract_mutation_type(rownames(ref))
ref_sorted <- ref[order(match(ref_mutation_types, mutation_order)), ]

input=read.table(sample_file,sep="\t",row.names = 1,header = T)
input_mutation_types <- extract_mutation_type(rownames(input))
input_sorted <- input[order(match(input_mutation_types, mutation_order)), , drop = FALSE]
input_sorted <- input_sorted[match(rownames(ref_sorted), rownames(input_sorted)), , drop = FALSE]

input_normalized <- input_sorted
input_normalized$Count <- input_sorted$Count / sum(input_sorted$Count)

ref_normalized <- ref_sorted
ref_normalized$Signature1 <- ref_sorted$Signature1 / sum(ref_sorted$Signature1)
ref_normalized$Signature2 <- ref_sorted$Signature2 / sum(ref_sorted$Signature2)  
ref_normalized$Signature3 <- ref_sorted$Signature3 / sum(ref_sorted$Signature3)

signature <- whichSignatures(
  tumor.ref = as.data.frame(t(input_normalized)),
  sample.id = "Count", 
  signatures.ref = as.data.frame(t(ref_normalized)),
  contexts.needed = FALSE,  
  tri.counts.method = 'default'
)

weights_vector <- as.vector(t(signature$weights))

result_df <- data.frame(
  Sample = sample_id,
  Signature = c("Signature1", "Signature2", "Signature3"),
  Proportion = weights_vector
)

proportion_file<-paste0(outpath, "/signature_proportions.txt")
write.table(result_df, 
            file = proportion_file, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
            
own_database=reference_error_profile
ref=read.table(own_database,sep="\t",row.names = 1,header = T)
sample_file=error_profile_file
ref_mutation_types <- extract_mutation_type(rownames(ref))
ref_sorted <- ref[order(match(ref_mutation_types, mutation_order)), ]

input=read.table(sample_file,sep="\t",row.names = 1,header = T)
input_mutation_types <- extract_mutation_type(rownames(input))
input_sorted <- input[order(match(input_mutation_types, mutation_order)), , drop = FALSE]
input_sorted <- input_sorted[match(rownames(ref_sorted), rownames(input_sorted)), , drop = FALSE]

input_normalized <- input_sorted
input_normalized$Count <- input_sorted$Count / sum(input_sorted$Count)

ref_normalized <- ref_sorted
ref_normalized$Signature1 <- ref_sorted$Signature1 / sum(ref_sorted$Signature1)
ref_normalized$Signature2 <- ref_sorted$Signature2 / sum(ref_sorted$Signature2)  
ref_normalized$Signature3 <- ref_sorted$Signature3 / sum(ref_sorted$Signature3)

signature <- whichSignatures(
  tumor.ref = as.data.frame(t(input_normalized)),
  sample.id = "Count", 
  signatures.ref = as.data.frame(t(ref_normalized)),
  contexts.needed = FALSE,  
  tri.counts.method = 'default'
)

weights_vector <- as.vector(t(signature$weights))

result_df <- data.frame(
  Sample = sample_id,
  Signature = c("Signature1", "Signature2", "Signature3"),
  Proportion = weights_vector
)

proportion_file<-paste0(outpath, "/signature_proportions.txt")
write.table(result_df, 
            file = proportion_file, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
            
######## get hFDR
df<-read.csv(candidate_df_file,header = T,sep="\t")

candidate_df<-read.csv(candidate_df_file,header = T,sep="\t")
if (!("identifier" %in% colnames(candidate_df))) {
  colnames(candidate_df)[colnames(candidate_df) == "X.identifier"] <- "identifier"
}

candidate_df$af<-as.numeric(candidate_df$consensus_alt_allele_count)/as.numeric(candidate_df$consensus_UMI_count)

sig_eps=0.01
print(reference_error_profile)
false_sig <- read.table(reference_error_profile,header = TRUE,row.names = 1)
col_sums <- colSums(false_sig)
normalized_refine_sig <- false_sig / col_sums

normalized_refine_sig$MutationType<-rownames(false_sig)
print(normalized_refine_sig)
print(dim(normalized_refine_sig))

proportion=read.table(proportion_file,sep="\t",header=TRUE)
beta_a=0.107
beta_b=47.58
for (i in 1:nrow(candidate_df)) {
    dp<-as.numeric(candidate_df[i,"consensus_UMI_count"])
    altreads<-as.numeric(candidate_df[i,"consensus_alt_allele_count"])
    af<-as.numeric(candidate_df[i,"af"])
    proportions <- proportion[proportion$Sample == sample_id, ]

    back_info <- handle_per_site_3artifact_remove_weight(candidate_df[i,], normalized_refine_sig, proportions,sig_eps)

    for(col_name in names(back_info)) {
      # print(col_name)
      candidate_df[i, col_name] <- as.numeric(back_info[[col_name]])
    }
    # #**** NOTES: in this script, we use the real artifact distribution to evaluate alpha, and beta distribution to evaluate beta.
    alpha =1 - pbbinom(altreads-1, dp, beta_a, beta_b)
    beta=(dp-altreads+1) / (dp+1)

  candidate_df[i, "refine_hFDR"] <- alpha / (alpha + beta * candidate_df[i, "rweigh"] * 1/49)

  candidate_df[i, "alpha"] <- as.numeric(alpha)
  candidate_df[i, "beta"] <- as.numeric(beta)
}


candidate_df_no_AG <- candidate_df[!grepl("A>G", candidate_df$RNAMutationType), ]
write.table(candidate_df_no_AG,
              paste0(outpath, "/features_with_hFDR_filterAtoG.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)


suffix="_remove_AtoG"
write.table(subset(candidate_df_no_AG, refine_hFDR <= 0.8),
              paste0(outpath, "/all_RF_filter_hFDR_less_08_feature", suffix, ".txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)

              