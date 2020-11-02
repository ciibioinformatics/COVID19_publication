rm(list = ls())

#Set working directory
setwd("~/Documents/Serochip/covid_pub")

library(reshape2)
library(ggplot2)
library(dplyr)

#Load wuhan epitopes from significant peptides
wuhan_epitopes <- read.csv("wuhan_epitopes_from_sig_pep_IgG_covid_pub.csv")
wuhan_epitopes$X <- NULL
colnames(wuhan_epitopes)[which(colnames(wuhan_epitopes) == "Other_CoV_218278_A12_IgG")] <- "SARS_2003_218278_A12_IgG"
plot_data <- wuhan_epitopes[,c(1,9:ncol(wuhan_epitopes))]
plot_data[,2:ncol(plot_data)] <- log(plot_data[,2:ncol(plot_data)])

#Load file containing shortlisted epitopes
epitopes_shortlist <- read.csv("shortlisted_epitope_frequency.csv")
epitopes_shortlist <- epitopes_shortlist[,1:2]
colnames(epitopes_shortlist) <- c("EPITOPE_NAME","EPITOPE_SEQ")
epitopes_shortlist$EPITOPE_NAME <- factor(epitopes_shortlist$EPITOPE_NAME, levels = unique(epitopes_shortlist$EPITOPE_NAME))
#Add signal values for shortlisted epitopes
epitopes_shortlist <- left_join(epitopes_shortlist, plot_data, by = "EPITOPE_SEQ")
#Melt the data to long format
plot_data_melt <- melt(epitopes_shortlist)
sample_names <- as.character(plot_data_melt$variable)
sample_names[which(grepl("Asymptomatic",sample_names))] <- "Asymptomatic COVID"
sample_names[which(grepl("Healthy_Controls_China",sample_names))] <- "Healthy controls"
#sample_names[which(grepl("Healthy_Controls_USA",sample_names))] <- "Healthy_Controls_USA"
sample_names[which(grepl("Mild_COVID",sample_names))] <- "Mild-COVID"
sample_names[which(grepl("Other_CoV",sample_names))] <- "Other HCoV controls"
sample_names[which(grepl("SARS_2003",sample_names))] <- "SARS-CoV-1"
sample_names[which(grepl("Severe_COVID",sample_names))] <- "Severe-COVID"
plot_data_melt$case_type <- sample_names
plot_data_melt$variable <- NULL
#The box plots have to be arranged in a specific order
plot_data_melt$case_type <- factor(plot_data_melt$case_type, levels = c('Healthy controls','Other HCoV controls','SARS-CoV-1','Asymptomatic COVID','Mild-COVID','Severe-COVID'),ordered = TRUE)
#Make color blind pallette
cbPalette <- c("#CC79A7", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00")
#Boxplot
plot_data_melt %>%
  ggplot(aes(x=case_type,y=value, fill=case_type)) +
  geom_boxplot() + #geom_jitter(position=position_jitter(0.2), cex = 0.1, alpha = 0.3) +
  scale_fill_manual(values=cbPalette) +
  facet_wrap(~EPITOPE_NAME) +
  xlab("Case Type") +
  ylab("Reactivity") +
  guides(fill=guide_legend(title = "Case Type")) +
  theme(axis.text.x = element_blank())
# Change stripchart colors by groups
plot_data_melt %>%
  ggplot(aes(x=case_type,y=value, color=case_type)) +
  geom_jitter(position=position_jitter(0.2)) +
  scale_color_manual(values=cbPalette) +
  facet_wrap(~EPITOPE_NAME) +
  theme(axis.text.x = element_blank())
