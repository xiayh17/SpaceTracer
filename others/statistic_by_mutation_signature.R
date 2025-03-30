##Author: Zhirui Yang
##Date: 2024-07-26
##Details: This is one script to calculate hFDR.
args <- commandArgs(trailingOnly = TRUE)

if ("-h" %in% args || "--help" %in% args) {
  cat("Usage: Rscript script.R feature_all_result outpath [bkg_file lysis_profile pcr_profile]\n")
  cat("Arguments:\n")
  cat("feature_all_result: Path to features dataframe file, make sure your dataframe contains these four colnames: identifier, RNAMutationType, dp, alt_dp\n")
  cat("outpath: Path to false signature file\n")
  cat("bkg_file: Path to background frequency file \n")
  cat("error_profile_file1: Path to mutation signature ptofile of one kind of errors (eg: lysis) \n")
  cat("error_profile_file2: Path to mutation signature ptofile of other kind of errors (eg: pcr) \n")
  cat("Notes: you can find example lysis_profile/pcr_profile/bkg_file in codes\n")
  q(save = "no")
}


#===========parameters=================
candidate_df_file<-args[1]
outpath<-args[2]

if (length(args) ==5) {
  ## this is for whole input
  bkg_file<-args[3]
  error_profile_file1<-args[4]
  error_profile_file2<-args[5]
  } else {
    if (length(args) ==4){
      ## this is for those only have one artifact signature profile
      bkg_file<-args[3]
      error_profile_file1<-args[4]
      error_profile_file2<-""
    }
  }


if (!file.exists(outpath)) {
  dir.create(outpath)
  print(paste("Path", outpath, "created successfully"))
} else {
  print(paste("Path", outpath, "already exists"))
} 



#===========functions=================
library(pracma)
library(dplyr)


# calculate signature similarity
get.sig.score.2artifacts <- function(mutsigs, true.sig, artifact.sig1,artifact.sig2, eps=0.001) {
  test.spectrum <- as.numeric(mutsigs)
  
  sigs <- cbind(as.numeric(as.matrix(true.sig)), as.numeric(as.matrix(artifact.sig1)),as.numeric(as.matrix(artifact.sig2)))
  weights <- pracma::lsqnonneg(sigs, test.spectrum)$x
  weights <- weights + eps
  postp <- log10(artifact.sig1*weights[2]+artifact.sig2*weights[3]) - log10(true.sig*weights[1])
  list(postp=postp, weight.true=weights[1], weight.artifact1=weights[2], weight.artifact2=weights[3])
}

cos_sim <- function(x, y) {
  res <- x %*% y / (sqrt(x %*% x) * sqrt(y %*% y))
  # coerce matrix to numeric
  res <- as.numeric(res)
  return(res)
}

get.sig.score.1artifact <- function(mutsigs, true.sig, artifact.sig, eps=0.001) {
  test.spectrum <- as.numeric(mutsigs)
  
  sigs <- cbind(as.numeric(as.matrix(true.sig)), as.numeric(as.matrix(artifact.sig)))
  weights <- pracma::lsqnonneg(sigs, test.spectrum)$x
  weights <- weights + eps
  postp <- log10(artifact.sig*weights[2]) - log10(true.sig*weights[1])
  list(postp=postp, weight.true=weights[1], weight.artifact=weights[2])
}

# tidy input signature and return result for each row
handle_per_site_2artifacts <- function(row, artifact.sig1,artifact.sig2, eps) {
  reference_list=rownames(artifact.sig1)
  true.sig=data.frame("Count"=rep(1/192,192))
  rownames(true.sig)<-rownames(artifact.sig1)
  mut_name=row$RNAMutationType
  artifact.sig2_sorted <- artifact.sig2[reference_list, ]

  index_position  <- which(sapply(reference_list, function(x) x == mut_name))
  result_vector <- rep(0, length(reference_list))
  result_vector[index_position] <- 1
  result<-get.sig.score.2artifacts(result_vector,true.sig,artifact.sig1,artifact.sig2_sorted, eps)
  cosine_similarity1 <- cos_sim(result_vector,artifact.sig1[,1])
  cosine_similarity2 <- cos_sim(result_vector,artifact.sig2[,1])
  
  return(c(result$weight.artifact1,result$weight.artifact2, result$postp[mut_name, "Count"],cosine_similarity1,cosine_similarity2,result$weight.true))
}

handle_per_site_1artifact <- function(row, artifact.sig, eps) {
  reference_list=rownames(artifact.sig)
  true.sig=data.frame("Count"=rep(1/192,192))
  rownames(true.sig)<-rownames(artifact.sig)
  mut_name=row$RNAMutationType
  
  index_position  <- which(sapply(reference_list, function(x) x == mut_name))
  result_vector <- rep(0, length(reference_list))
  result_vector[index_position] <- 1
  result<-get.sig.score.1artifact(result_vector,true.sig,artifact.sig, eps)
  cosine_similarity <- cos_sim(result_vector,artifact.sig[,1])
  
  return(c(result$weight.artifact, result$postp[mut_name, "Count"],cosine_similarity,result$weight.true))
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


# generate signature profile
output_profile<-function(df,out_name){
  counts <- table(df$RNAMutationType)
  sig<- as.data.frame(counts)
  colnames(sig)<-c("MutationType","Count")
  
  default_list=get_default_labels("RNA")
  default_df=data.frame("MutationType" =default_list,"none"=rep(0,length((default_list))))
  merged_df<-left_join(default_df,sig,by="MutationType")
  
  merged_df$Count <- ifelse(is.na(merged_df$Count), 0, merged_df$Count)
  merged_df["none"]<-NULL
  # write.table(merged_df, paste(outpath,"/",file_name,"_",formatted_date,"_",out_name,sep=""), row.names=FALSE, col.names=TRUE,quote=FALSE,sep="\t")
  write.table(merged_df, paste(outpath,"/",file_name,out_name,sep=""), row.names=FALSE, col.names=TRUE,quote=FALSE,sep="\t")
  
}
#=============handle data=============
# date <- Sys.Date()
# formatted_date <- gsub("-", "_", date)

# get real artifacts distribution
bkg<-read.table(bkg_file,header=FALSE)
a<-hist(bkg$V1,breaks = seq(0,1,0.01))
df_hist<-data.frame("VAF"=a$mids,"prop"=a$counts/sum(a$counts))

# read features
candidate_df<-read.csv(candidate_df_file,header = T,sep="\t")
if (!("identifier" %in% colnames(candidate_df))) {
  colnames(candidate_df)[colnames(candidate_df) == "X.identifier"] <- "identifier"
}

## Please make sure your input file has these columns: "identifier","RNAMutationType","dp","alt_dp"
#candidate_df_short<-candidate_df[,c("identifier","RNAMutationType","dp","alt_dp")]
#candidate_df_short$af<-as.numeric(candidate_df_short$alt_dp)/as.numeric(candidate_df_short$dp)
candidate_df$af<-as.numeric(candidate_df$alt_dp)/as.numeric(candidate_df$dp)

sig_eps=0.01

# handle artifacts signature

false_sig1 <- read.table(error_profile_file1,header = TRUE,row.names = 1)
refine_sig1<-refine_counts(false_sig1)
col_sums <- colSums(refine_sig1)
normalized_refine_sig1 <- refine_sig1 / col_sums

refine_sig1$MutationType<-rownames(refine_sig1)

if (error_profile_file2!=""){
  false_sig2 <- read.table(error_profile_file2,header = TRUE,row.names = 1)
  refine_sig2<-refine_counts(false_sig2)
  col_sums <- colSums(refine_sig2)
  normalized_refine_sig2 <- refine_sig2 / col_sums
  
  refine_sig2$MutationType<-rownames(refine_sig2)

  mode="2artifacts"
  
} else {
  mode="1artifact"
}


na=99
nt=1

for (i in 1:nrow(candidate_df)) {
  dp<-as.numeric(candidate_df[i,"dp"])
  altreads<-as.numeric(candidate_df[i,"alt_dp"])
  af<-as.numeric(candidate_df[i,"af"])
  if (mode=="2artifacts") {
    back_info <- handle_per_site_2artifacts(candidate_df[i,], normalized_refine_sig1, normalized_refine_sig2,sig_eps)
    candidate_df[i, "true_similarity"] <- as.numeric(back_info[6])
    candidate_df[i, "false_similarity1"] <- as.numeric(back_info[1])
    candidate_df[i, "false_similarity2"] <- as.numeric(back_info[2])
    candidate_df[i, "rweigh"] <- as.numeric(10^-as.numeric(back_info[3]))
    candidate_df[i, "cosine_sim1"]<- as.numeric(back_info[4])
    candidate_df[i, "cosine_sim2"]<- as.numeric(back_info[5])
  } else {
    back_info <- handle_per_site_1artifact(candidate_df[i,], normalized_refine_sig1,sig_eps)
    candidate_df[i, "true_similarity"] <- as.numeric(back_info[4])
    candidate_df[i, "false_similarity1"] <- as.numeric(back_info[1])
    candidate_df[i, "rweigh"] <- as.numeric(10^-as.numeric(back_info[2]))
    candidate_df[i, "cosine_sim1"]<- as.numeric(back_info[3])
  }
  #**** NOTES: in this script, we use the real artifact distribution to evaluate alpha, and beta distribution to evaluate beta.
  alpha <- df_hist[which.max(df_hist$VAF   >= (af -0.005)), "prop"][1]
  beta=1-pbeta(altreads/dp,altreads+1,dp-altreads+1)
  
  hFDR <- alpha*na / (alpha*na + beta * candidate_df[i, "rweigh"] * nt)
  candidate_df[i, "alpha"] <-as.numeric(alpha)
  candidate_df[i, "beta"] <-as.numeric(beta)
  candidate_df[i, "hFDR"] <- hFDR
  
  refine_rescue.fdr99 <- alpha / (alpha + beta * candidate_df[i, "rweigh"] * 1/99)
  refine_rescue.fdr49 <- alpha / (alpha + beta * candidate_df[i, "rweigh"] * 1/49)
  refine_rescue.fdr29 <- alpha / (alpha + beta * candidate_df[i, "rweigh"] * 1/29)
  refine_rescue.fdr9 <- alpha / (alpha + beta * candidate_df[i, "rweigh"] * 1/9)
  refine_rescue.fdr1 <- alpha / (alpha + beta * candidate_df[i, "rweigh"] * 1/1)
  candidate_df[i, "refine_hFDR99"] <- refine_rescue.fdr99
  candidate_df[i, "refine_hFDR49"] <- refine_rescue.fdr49
  candidate_df[i, "refine_hFDR29"] <- refine_rescue.fdr29
  candidate_df[i, "refine_hFDR9"] <- refine_rescue.fdr9
  candidate_df[i, "refine_hFDR1"] <- refine_rescue.fdr1

}


file_name <- basename(candidate_df_file)
file_name <- sub(".txt", "", file_name)

# threshold of hFDR
threshold=0.1
candidate_false<-subset(candidate_df,subset=(hFDR > threshold))
candidate_true<-(subset(candidate_df,subset=(hFDR <= threshold)))
#write.table(candidate_true,paste(outpath,"/",file_name,"_candidate_true.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE,sep="\t")
#write.table(candidate_false,paste(outpath,"/",file_name,"_candidate_false.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE,sep="\t")
#output_profile(candidate_true,paste0(as.character(threshold),"_candidate_true_Sigcount.txt"))
#output_profile(candidate_false,paste0(as.character(threshold),"_candidate_false_Sigcount.txt"))

write.table(candidate_df,paste(outpath,"/",file_name,".add_hFDR.txt",sep=""),row.names = FALSE,col.names = TRUE,quote=FALSE,sep="\t")

png(paste(outpath,"/",file_name,"_histogram_of_hFDR.png",sep="")) 
hist(candidate_df$hFDR, main = "Histogram of hFDR score", xlab = "hFDR")
dev.off()
