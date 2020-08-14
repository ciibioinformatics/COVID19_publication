rm(list = ls())

#Set working directory
setwd("~/Documents/Serochip/covid_pub")

args <- vector()
args[1] <- "IgG"
args[2] <- "10000"
args[3] <- paste(args[1],"_files_aggregated.txt",sep = "")
#args[3] <- "Herpes_IgG_Epitopes_test.csv"
#library(data.table)
library(dplyr)
library(readr)
library(stringr)
#library(tidyverse)
#library("Biostrings")

#Following code is to analysis the data from new Coronavirus HCoV chip design NTX127
#Aggregate files for IgG were generated from the raw data using the preexisting code from the pipeline using programe aggregate_signals.pl

#Make column names informative based on the metadata
#Load metadata
covid_metadata <- read_csv("covid_metadata3.csv")
#Drop info for convalescent samples since those are not inlcuded in this analysis
covid_metadata <- covid_metadata[!(grepl("CONVALESCENT",covid_metadata$Case_Type)),]
#It is important the sequence of columns in the data and samples in metadata are in the same order
covid_metadata <- covid_metadata[order(covid_metadata$BARCODE, covid_metadata$Subarray),]
#Make new column names
covid_metadata$new_colnames <- paste(covid_metadata$Subgroup_name, covid_metadata$BARCODE, covid_metadata$Subarray, "IgG",sep = "_")

#Load aggregate files
aggregate_file <- read_delim(args[3], delim = "\t")
data_cols <- aggregate_file[,18:ncol(aggregate_file)]
data_cols <- data_cols[,order(colnames(data_cols))]
aggregate_file <- data.frame(aggregate_file[,1:17],data_cols)
#Remove X from column name
colnames(aggregate_file)[18:ncol(aggregate_file)] <- gsub("^X", "",  colnames(aggregate_file)[18:ncol(aggregate_file)])
#Assign new column names
colnames(aggregate_file)[18:ncol(aggregate_file)] <- covid_metadata$new_colnames

#Drop healthy samples from USA
aggregate_file <- aggregate_file %>% select(-contains("USA"))
write.csv(aggregate_file, "IgG_aggregate_file_covid_pub.csv", row.names = F)
#Store random probes in a separate variable and then drop them
aggregate_random <- aggregate_file[which(grepl("random_probe",aggregate_file$PROBE_ID)==T),]
#Filter aggregate file 
aggregate_filtered <- aggregate_file[apply(aggregate_file[,18:(ncol(aggregate_file))]>=as.numeric(args[2]), 1, any),]
#84148 peptide probes remaining after filtration (THIS INCLUDES RANDOM PROBES)
#Drop random probes
aggregate_filtered <- aggregate_filtered[-which(grepl("random_probe",aggregate_filtered$PROBE_ID)==T),]
#Separate aggregate_filtered into redundant and non-redundant
aggregate_redund <- aggregate_filtered[aggregate_filtered$SEQ_ID == "REDUNDANT",]
aggregate_nonredund <- aggregate_filtered[aggregate_filtered$SEQ_ID != "REDUNDANT",]
#Save memory by deleting data frames
#rm(aggregate_file)
#Add all matching probes from aggregate files to correspondence key
#Load Correspondence key
correspondence_key <- read_delim("/Users/shreyas/Documents/Serochip/Array_200129_HCoV/Design_files/NTX-127_all_11_files_correspondence_key.txt", delim = "\t")
#Change column name for peptide sequence in correspondence key so that it is easier to merge
colnames(correspondence_key)[which(colnames(correspondence_key) == "PEPTIDE_SEQUENCE")] <- "PROBE_SEQUENCE"
#Select rows from correspondence key that are present in filtered aggregate_redund file
correspondence_key <- correspondence_key[(correspondence_key$PROBE_SEQUENCE %in% aggregate_redund$PROBE_SEQUENCE),]
#Add all matching probes from aggregate_redund to correspondence key
#Use left_join from dplyr package
merged_redund <- left_join(correspondence_key,aggregate_redund, by=c("PROBE_SEQUENCE","SEQ_ID","POSITION"))
#rm(sndf_file)
#Assign ORIGINAL_PROBE_ID to PROBE_ID in the merged file since this field will be used to make epitopes and then with the header map
merged_redund$PROBE_ID <- merged_redund$ORIGINAL_PROBE_ID
#Select relevant columns from merged_redund
merged_redund <- merged_redund[,c(1,6,8,17,21:ncol(merged_redund))]
#Select same columns from aggregate_nonredund
aggregate_nonredund <- aggregate_nonredund[,c(5,6,2,13,18:ncol(aggregate_nonredund))]
#Combine data frames. this will serve as the main analysis data frame
analysis_df <- rbind(aggregate_nonredund, merged_redund)
#Delete aggregate_nonredund
rm(aggregate_nonredund)
rm(merged_redund)
#Drop correspondence key
rm(correspondence_key)

#Select rows that have at least one column above the threshold
#filtered_file <- analysis_df[apply(analysis_df[,5:(ncol(analysis_df))]>=as.numeric(args[2]), 1, any),]
filtered_file <- analysis_df
#Drop analysis_df
#rm(analysis_df)

#Add colummn for probe_number
filtered_file$PROBE_NUMBER <- str_split_fixed(filtered_file$PROBE_ID,";",2)[,2]
#Add column for sequence source file
filtered_file$SEQ_NUMBER <- str_split_fixed(filtered_file$PROBE_ID,";",2)[,1]

#Load Header Map
header_map <- read_delim("/Users/shreyas/Documents/Serochip/Array_200129_HCoV/Design_files/Corona_chip_header_map_no_seq.txt", delim = "\t")
#Change column names of header map
colnames(header_map) <- c("SEQ_SOURCE","SEQ_NUMBER","SEQ_INFO")
#Add sequence details from header map to filtered file
filtered_file <- left_join(filtered_file, header_map, by = "SEQ_NUMBER")
#Drop data frame for header maps
rm(header_map)
#Sort data frame
filtered_file$PROBE_NUMBER <- as.numeric(filtered_file$PROBE_NUMBER)
filtered_file <- filtered_file[order(filtered_file$SEQ_NUMBER,filtered_file$PROBE_NUMBER),]

#Rearrange columns
filtered_file <- filtered_file[,c(1:4,(ncol(filtered_file)-3):ncol(filtered_file),5:(ncol(filtered_file)-4))]
#Save filtered file
write.csv(filtered_file,"Reports/HCoV_200129_Mar26_covid_pub_IgG_EXPT_MAIN.fltred_wAtLeast1Sig_gt10000_IgG.wDescr.sorted.csv", row.names = F)

#Separate filtered_file by individual agent
agent_source <- unique(filtered_file$SEQ_SOURCE)
agent_source <- c(agent_source[10],agent_source[1:9])
for(i in 1:length(agent_source)){
  #i <- 9
  #Make data frame "tmp_df" that is for a an individual agent
  print(paste("Starting analysis for ",agent_source[i]))
  tmp_df <- filtered_file[which(filtered_file$SEQ_SOURCE == agent_source[i]),]
  tmp_df <- tmp_df[order(tmp_df$SEQ_NUMBER, tmp_df$PROBE_NUMBER),]
  #Save this data frame in a file that will be provided in reports
  outfile_name <- paste("Reports/IgG/PEPTIDES_SPLIT_BY_AGENT/",agent_source[i],"200129_covid_pub_EXPT_MAIN.fltred_wAtLeast1Sig_gt10000_IgG.wDescr.sorted.csv",sep = "")
  write.csv(tmp_df, outfile_name, row.names = F)
  #Make a list of individual sequences in an agent
  ind_seq <- unique(tmp_df$SEQ_NUMBER)
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
  #Drop column for probe number
  tmp_df$PROBE_NUMBER <- NULL
  ####REARRANGE COLUMNS####
  test <- tmp_df[,8:ncol(tmp_df)]
  test <- test[,order(colnames(test))]
  tmp_df <- cbind(tmp_df[,1:7],test)
  
  #Assign apitope file name
  epitope_outfile <- paste("Reports/IgG/PEPTIDES_SPLIT_BY_AGENT/epitopes/",agent_source[i],"_200129_covid_pub_IgG_gt_10000_EXPT_MAIN.min3consecutive.nonredund.csv",sep = "")
  print(paste("Epitope detection completed for ",agent_source[i]))
  print(paste("Writing results to file for ",agent_source[i]))
  write.csv(tmp_df, epitope_outfile, row.names = F)
  print(paste("Results saved for ",agent_source[i]))
  #assign(agent_source[i], tmp_df)
}

#Extract wuhan epitopes from VIPR data
vipr_epi <- read_csv("HCoV_Date20200326_ExpIDMar26_allsamples/Reports_3/PEPTIDES_SPLIT_BY_AGENT/epitopes/VIPRdb_CoronaProteinFasta_no_duplicates_200129_Mar26_allsamples_EXPT_MAIN.min3consecutive.nonredund.csv")
wuhan_vipr <- vipr_epi[which(grepl("wuhan",vipr_epi$SEQ_INFO,ignore.case = T)),]
rm(vipr_epi) #save memory
wuhan_vipr$organism <- str_split_fixed(wuhan_vipr$SEQ_INFO,"[|]",6)[,4]
wuhan_vipr$strain <- str_split_fixed(wuhan_vipr$SEQ_INFO,"[|]",6)[,5]
wuhan_vipr$protein <- str_split_fixed(wuhan_vipr$SEQ_INFO,"[|]",6)[,6]
#Rearrange columns
wuhan_vipr <- wuhan_vipr[,c(1:7,140:142,8:139)]
wuhan_vipr$SEQ_INFO <- str_split_fixed(wuhan_vipr$SEQ_INFO,"[|]",6)[,2]
wuhan_vipr_epitope_outfile <- paste("HCoV_Date20200326_ExpIDMar26_allsamples/Reports/PEPTIDES_SPLIT_BY_AGENT/epitopes/","Wuhan_VIPR","_200129_Mar26_allsamples_EXPT_MAIN.min3consecutive.nonredund.csv")
write_csv(wuhan_vipr,wuhan_vipr_epitope_outfile)

wuhan_epi <- read_csv("HCoV_Date20200326_ExpIDMar26_allsamples/Reports_3/PEPTIDES_SPLIT_BY_AGENT/epitopes/Wuhan_seafood_market_CoV_MN908947_CDS_proteome_200129_Mar26_allsamples_EXPT_MAIN.min3consecutive.nonredund.csv")

