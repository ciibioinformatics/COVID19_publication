rm(list = ls())

#Set working directory
setwd("~/Documents/Serochip/covid_pub")

library(readr)

#Load significant peptides from mds analysis for SARS-CoV2 vs all others
sig_pep_sars_cov2_vs_control <- read.csv("sars_cov2_vs_control_sig_pep_covid_pub.csv")

#Load filtered peptide file for Wuhan
wuhan_filtered <- read.csv("Reports/IgG/PEPTIDES_SPLIT_BY_AGENT/Wuhan_seafood_market_CoV_MN908947_CDS_proteome200129_covid_pub_EXPT_MAIN.fltred_wAtLeast1Sig_gt10000_IgG.wDescr.sorted.csv")

#Load metadata
covid_metadata <- read_csv("/Users/shreyas/Documents/Serochip/Array_200129_HCoV/MDS/covid_metadata3.csv")
#Drop info for convalescent samples  and healthy controls from USA since those are not inlcuded in this analysis
covid_metadata <- covid_metadata[!(grepl("CONVALESCENT",covid_metadata$Case_Type)),]
covid_metadata <- covid_metadata[!(grepl("Healthy_Controls_USA",covid_metadata$Subgroup_name)),]
#It is important the sequence of columns in the data and samples in metadata are in the same order
covid_metadata <- covid_metadata[order(covid_metadata$BARCODE, covid_metadata$Subarray),]
covid_metadata$sample_colnames <- colnames(wuhan_filtered[9:ncol(wuhan_filtered)])
#Shortlist peptides
wuhan_sig_pep <- wuhan_filtered[wuhan_filtered$PROBE_SEQUENCE %in% sig_pep_sars_cov2_vs_control$PROBE_SEQUENCE,]

#Make a list of individual sequences in an agent
ind_seq <- unique(wuhan_sig_pep$SEQ_NUMBER)
tmp_df <- wuhan_sig_pep
tmp_df$EPITOPE_SEQ <- NA
for(j in 1:length(ind_seq)){
  #Make data frame "tmp_df_seq" for an indicidual sequence in an agent
  tmp_df_seq <- tmp_df[which(tmp_df$SEQ_NUMBER == ind_seq[j]),1:8]
  if(nrow(tmp_df_seq) < 3){
    next
  }
  tmp_probe_no <- tmp_df_seq$PROBE_NUMBER
  y <- sort(tmp_probe_no)
  g <- cumsum(c(1, abs(y[-length(y)] - y[-1]) > 1))
  epitope_df <- data.frame(t(rbind(g, y)))
  colnames(epitope_df) <- c("index","probe_no")
  probe_freq <- data.frame(table(epitope_df$index))
  probe_freq <- probe_freq[which(probe_freq$Freq > 2),]
  if (nrow(probe_freq) ==0){
    next
  }
  #Assign epitope names
  for (k in 1:nrow(probe_freq)){
    probe_freq$epitope_name[k] <- paste("Epitope_",k,sep = "")
  }
  #identify locations above the length threshold
  freq_threshold <- as.character(probe_freq$Var1[which(probe_freq$Freq > 2)])
  #Select probe numbers that are are associated with higher frequency of g
  epitope_df <- epitope_df[epitope_df$index %in% freq_threshold,]
  epitope_df$epitope_name <- probe_freq$epitope_name[match(epitope_df$index, probe_freq$Var1)]
  epitope_df$PROBE_SEQUENCE <- tmp_df_seq$PROBE_SEQUENCE[match(epitope_df$probe_no, tmp_df_seq$PROBE_NUMBER)]
  #Make a vector of uniaue epitopes
  uni_epitopes <- unique(epitope_df$epitope_name)
  #Create empty data frame like the filtered_file in which all data frames for unique SEQ_IDs will be concatenated with the epitope sequences
  df_with_epitopes <- epitope_df
  df_with_epitopes$EPITOPE_SEQ <- NA
  df_with_epitopes <- df_with_epitopes[FALSE,]
  #Make a vector of uniaue epitopes
  uni_epitopes <- unique(epitope_df$epitope_name)
  #The following loop takes all the probes occuring in sequential manner and reassembles them into a long sequence
  for (l in 1:length(uni_epitopes)) {
    tmp_df_epitope <- epitope_df[epitope_df$epitope_name == uni_epitopes[l],]
    tmp_df_epitope$PROBE_SEQUENCE <- as.character(tmp_df_epitope$PROBE_SEQUENCE)
    rownames(tmp_df_epitope) <- c(1:nrow(tmp_df_epitope)) 
    for(m in 1:nrow(tmp_df_epitope)){
      if(m == 1){
        tmp_df_epitope$EPITOPE_SEQ[m] <- tmp_df_epitope$PROBE_SEQUENCE[m]
      } else {
        tmp_df_epitope$EPITOPE_SEQ[m] <- paste(tmp_df_epitope$EPITOPE_SEQ[m-1],substr(tmp_df_epitope$PROBE_SEQUENCE[m],12,12),sep = "") 
      }
    }
    tmp_df_epitope$EPITOPE_SEQ <- tmp_df_epitope$EPITOPE_SEQ[nrow(tmp_df_epitope)]
    df_with_epitopes <- rbind(df_with_epitopes,tmp_df_epitope)
    tmp_df_seq$EPITOPE_SEQ <- df_with_epitopes$EPITOPE_SEQ[match(tmp_df_seq$PROBE_NUMBER, df_with_epitopes$probe_no)]
  }
  tmp_df$EPITOPE_SEQ[which(tmp_df$SEQ_NUMBER == ind_seq[j])] <- tmp_df_seq$EPITOPE_SEQ[match(tmp_df$PROBE_NUMBER[which(tmp_df$SEQ_NUMBER == ind_seq[j])], tmp_df_seq$PROBE_NUMBER)]
}
#Rearrange columns
tmp_df <- tmp_df[,c(ncol(tmp_df),2:(ncol(tmp_df)-1))]
#Drop probes which are not part of the epitope
tmp_df <- na.omit(tmp_df)
#Store data in a epitope data frame
df_with_epitopes <- tmp_df
write.csv(df_with_epitopes, "wuhan_epitopes_from_sig_pep_covid_pub.csv")

#Making positive vs negative call
#Convert data to binary. For a single epitope within a sample, if there are three consecutive peptides that are above the threshold,
#or if there are 3 out of 4 peptides (allowing a gap of one), a positive call will be made
tmp_df <- df_with_epitopes[,9:ncol(df_with_epitopes)]
tmp_df[tmp_df > 9999] <- 1
tmp_df[tmp_df != 1] <- 0
epitopes_bin <- cbind(df_with_epitopes[,1:8],tmp_df)

#make empty dataframe for positive negative calls
call_df <- df_with_epitopes[FALSE,]


epitopes <- unique(df_with_epitopes$EPITOPE_SEQ)
for(i in 1:length(epitopes)){
  tmp_df <- epitopes_bin[which(epitopes_bin$EPITOPE_SEQ == epitopes[i]),]
  call_df[i,1:8] <- tmp_df[1,1:8]
  tmp_df <- tmp_df[,9:ncol(tmp_df)]
  if (nrow(tmp_df) < 3){
    call_df[i,9:ncol(call_df)] <- 0
  } else if (nrow(tmp_df) <= 4) {
    call_df[i,(which(colSums(tmp_df) > 2))+8] <- 1
    call_df[i,(which(colSums(tmp_df) < 3))+8] <- 0
  } else {
    list_vec <- list()
    for(j in 1:ncol(tmp_df)){
      sum_vec <- vector()
      tmp_vec <- tmp_df[,j]
      for(k in 1:length(tmp_vec)){
        if(k+3 <= length(tmp_vec)){
          sum_vec <- c(sum_vec,sum(tmp_vec[k:(k+3)]))
        }
      }
      list_vec[[j]] <- sum_vec
      if(length(which(sum_vec > 2)) > 0){
        call_df[i,j+8] <- 1
      }else{
        call_df[i,j+8] <- 0
      }
    }
  }
}

tmp <- call_df[,9:ncol(call_df)]
tmp <- tmp[,order(names(tmp))]
call_df <- cbind(call_df[,1:8],tmp)
call_df <- call_df[which(rowSums(tmp) != 0),]

# nm_epitopes <- read_csv("nm_epitopes.csv")
# nm_epitopes <- nm_epitopes$nm_epitope
# 
# tmp_df <- call_df
# tmp_df$nm_epitope <- NA
# 
# for(i in 1:nrow(tmp_df)){
#   if (length(nm_epitopes[which(grepl(tmp_df$PROBE_SEQUENCE[i], nm_epitopes))]) > 0){
#     tmp_df$nm_epitope[i] <- nm_epitopes[which(grepl(tmp_df$PROBE_SEQUENCE[i], nm_epitopes))]
#   }
# }

# tmp_df <- tmp_df[,c(ncol(tmp_df),1:(ncol(tmp_df)-1))]
# tmp <- na.omit(tmp_df)
# call_df <- tmp_df
#Select relevant columns
call_df <- call_df[,c(1,6:ncol(call_df))]
write.csv(call_df,"epitope_call_covid_pub.csv")
