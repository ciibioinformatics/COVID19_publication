rm(list = ls())

#Set working directory
setwd("~/Documents/Serochip/covid_pub")

library(readr)
library(tidyr)
library(stringr)
library(dplyr)

epitope_call <- read.csv("epitope_call_IgG_covid_pub.csv") 
epitope_call$X <- NULL
epitope_call <- epitope_call[,c(1,5:ncol(epitope_call))]
#A12 of 218278 belongs to SARS 2003 subgroup as well
colnames(epitope_call)[which(colnames(epitope_call) == "Other_CoV_218278_A12_IgG")] <- "SARS_2003_218278_A12_IgG"
#Select the epitopes that are not in Nischay's list
call_df <- epitope_call
rownames(call_df) <- call_df$EPITOPE_SEQ

call_df$EPITOPE_SEQ <- NULL
call_df <- t(call_df)
rownames(call_df) <- gsub("Healthy_Controls","Healthy-Controls",rownames(call_df))
rownames(call_df) <- gsub("Mild_COVID","Mild-COVID",rownames(call_df))
rownames(call_df) <- gsub("Severe_COVID","Severe-COVID",rownames(call_df))
rownames(call_df) <- gsub("Severe_COVID","Severe-COVID",rownames(call_df))
#rownames(call_df) <- gsub("Asymptomatic_First","Asymptomatic-First",rownames(call_df))
#rownames(call_df) <- gsub("Asymptomatic_Second","Asymptomatic-Second",rownames(call_df))

call_df <- as.data.frame(call_df)
call_df$Barcode <- str_split_fixed(rownames(call_df),"_",5)[,3]
call_df$Subarray <- str_split_fixed(rownames(call_df),"_",5)[,4]
call_df$category1 <- str_split_fixed(rownames(call_df),"_",4)[,1]
call_df$category2 <- str_split_fixed(rownames(call_df),"_",4)[,2]

call_df <- call_df[,c((ncol(call_df)-3):ncol(call_df),1:(ncol(call_df)-4))]

rownames(call_df) <- 1:nrow(call_df)

#make category single variable
call_df$category <- paste(call_df$category1, call_df$category2, sep = "_")
call_df <- call_df[,c(ncol(call_df),1:(ncol(call_df)-1))]
#call_df$category[which(grepl("Asymptomatic",call_df$category))] <- "Asymtomatic"

category_size <- call_df %>% group_by(category) %>% tally()
antibody_freq <- call_df[,c(1,6:ncol(call_df))]
antibody_freq <- antibody_freq %>% group_by(category) %>% summarise_all(sum)
antibody_freq <- as.data.frame(antibody_freq)
rownames(antibody_freq) <- antibody_freq$category
antibody_freq$category <- NULL
antibody_freq <- as.data.frame(t(antibody_freq))

antibody_report <- antibody_freq
for(i in 1:nrow(category_size)){
  antibody_report[,i] <- paste(antibody_freq[,i],"/",category_size$n[i]," (",ceiling((antibody_report[,i]/category_size$n[i])*100),"%)",sep = "")
}

write.csv(antibody_report, "IgG_antibody_frequency_covid_pub_updated_Jul20.csv") #This is table S5 in paper

antibody_report$epitope_seq <- rownames(antibody_report)
shortlisted_epitopes <- read.csv("shortlisted_epitopes.csv")
shortlisted_epitopes <- shortlisted_epitopes[103:131,]
shortlisted_epitopes <- left_join(shortlisted_epitopes, antibody_report, by = "epitope_seq")
write.csv(shortlisted_epitopes,"Table1_covid_pub.csv")
