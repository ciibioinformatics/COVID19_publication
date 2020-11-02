rm(list = ls())

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
#BiocManager::install("DESeq2")
#BiocManager::install("SummarizedExperiment")

library(edgeR)
library(readr)
library(ggplot2)
library(ggthemes)
library(plyr)
library(MASS)

#Set working directory
setwd("~/Documents/Serochip/covid_pub")

args <- vector()
args[1] <- "IgG_aggregate_file_covid_pub.csv"
args[2] <- "10000"

#Load aggregate files
aggregate_file <- read.csv(args[1])
#Filter aggregate file
aggregate_filtered <- aggregate_file[apply(aggregate_file[,18:(ncol(aggregate_file))]>=as.numeric(args[2]), 1, any),]

#For differential analysis between groups, functions from package "limma" have been used
#To be able to load values from aggregate file into DGEList object, make dataset containing probe sequences and signal intensities
#from filtered aggregate file
analysis_df <- aggregate_filtered[,c(13,18:ncol(aggregate_filtered))]
#Assign probe ids as row names since they are unique identifier
rownames(analysis_df) <- analysis_df$PROBE_ID
analysis_df$PROBE_ID <- NULL
analysis_df <- ceiling(analysis_df)

#Load metadata
covid_metadata <- read_csv("/Users/shreyas/Documents/Serochip/Array_200129_HCoV/MDS/covid_metadata3.csv")
#Drop info for convalescent samples  and healthy controls from USA since those are not inlcuded in this analysis
covid_metadata <- covid_metadata[!(grepl("CONVALESCENT",covid_metadata$Case_Type)),]
covid_metadata <- covid_metadata[!(grepl("Healthy_Controls_USA",covid_metadata$Subgroup_name)),]
#It is important the sequence of columns in the data and samples in metadata are in the same order
covid_metadata <- covid_metadata[order(covid_metadata$BARCODE, covid_metadata$Subarray),]
#For comparison between SARS-COV2 and control, make new variable agent_2
covid_metadata$agent_2 <- covid_metadata$AGENT
covid_metadata$agent_2[which(grepl("SARS-COV-2",covid_metadata$agent_2))] <- "SARS_COV_2"
covid_metadata$agent_2[which(!grepl("SARS_COV_2",covid_metadata$agent_2))] <- "Control"

covid_metadata$new_colnames <- paste(colnames(aggregate_file[,18:ncol(aggregate_file)]),covid_metadata$agent_2,sep = "_")
colnames(analysis_df) <- covid_metadata$new_colnames
covid_metadata <- covid_metadata[,c(ncol(covid_metadata),1:(ncol(covid_metadata)-1))]
#covid_metadata$group <- interaction(covid_metadata$AGENT, covid_metadata$Case_Type)

colnames(analysis_df) <- covid_metadata$new_colnames

#Make DGEList object
#covid_metadata$Stage <- factor(covid_metadata$Stage)
y <- DGEList(analysis_df, group = as.character(covid_metadata$agent_2))
y <- calcNormFactors(y)
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
y <- estimateDisp(y, design)
fit1 <- glmQLFit(y, design)

#MAke combinarions for contrast
con_SARS_vs_Control <- "SARS_COV_2 - Control"

#Generate linear model
#contrast_pair <- as.expression(noquote(contrast_list[1]))
con <- makeContrasts(con_SARS_vs_Control, levels=design)
qlf <- glmQLFTest(fit1, contrast=con)
#Shortlist significant peptides based on adjusted p-value
sig_pep <- qlf$table
sig_pep$p_adjusted <- p.adjust(sig_pep$PValue, method = "BH")
sig_pep <- sig_pep[which(sig_pep$p_adjusted <= 0.05),]
#Print summary 
print(summary(decideTests(qlf)))
#Save decideTests results in a data frame
test <- data.frame(decideTests(qlf))
test$probe_name <- rownames(test)
rownames(test) <- 1:nrow(test)
length(which(test$X.1.Control.1.SARS_COV_2 == 1))
test <- test[which(test$X.1.Control.1.SARS_COV_2 != 0),]
#Make variable for storing significant probes
probe_var <- test$probe_name
sig_pep_df <- aggregate_filtered[which(aggregate_filtered$PROBE_ID %in% probe_var),]
write.csv(sig_pep_df, "sars_cov2_vs_control_sig_pep_covid_pub.csv", row.names = F)

#Function to generate pattern to be used with grepl
gen_pattern <- function(contrast_pair){
  contrast_pair <- as.character(contrast_pair)
  #Store elements of pair in a variable
  e1 <- paste(strsplit(contrast_pair," - ")[[1]][1],"$",sep = "")
  e2 <- paste(strsplit(contrast_pair," - ")[[1]][2],"$",sep = "")
  grep_pattern <- paste(e1,e2,sep = "|")
  return(grep_pattern)
}

name_combo <- gen_pattern(con_SARS_vs_Control)

#Function to generate log normalized dataframe with significant peptides
#Takes list of significant peptides as input
log_sigdf <- function(probe_var,name_combo){
  sig_df <- analysis_df[which(rownames(analysis_df) %in% probe_var),]
  #Extract columns of interest
  sig_df <- sig_df[,grepl(name_combo,colnames(sig_df),perl = T)]
  #sig_df <- analysis_df
  sig_log_df <- data.frame(log(sig_df))
  return(sig_log_df)
}

sig_log_df <- log_sigdf(probe_var, name_combo)

df_mds <- dist(t(sig_log_df))
fit2 <- isoMDS(df_mds) 
x <- fit2$points[,1]
y <- fit2$points[,2]
min_lim = min(min(x),min(y))
max_lim = max(max(x),max(y))
plot_df = data.frame(x =x, y =y)
plot_df$id = names(x)
color_var <- as.factor(covid_metadata$agent_2[grepl(name_combo,covid_metadata$agent_2,perl = T)])
shape_var <- as.factor(covid_metadata$agent_2[grepl(name_combo,covid_metadata$agent_2,perl = T)])
plot_df$color_var <- color_var
plot_df$shape_var <- shape_var
plot_df$hull = as.factor(paste0(color_var,shape_var))
# Points lying on convex hull
find_hull <- function(df) df[chull(df$x, df$y), ]
hulls <- ddply(plot_df, "hull", find_hull)

p <- ggplot(plot_df, aes(x= x, y=y,colour=color_var,shape=color_var,label = id)) 
p
png <- p + 
  #geom_point(alpha=5,size = 1) +
  geom_point(aes(shape=shape_var, color=shape_var)) +
  scale_shape_manual(values=(16:(length(levels(plot_df$shape_var))+16))) +
  scale_color_manual(values=c('#D55E00','#009E73')) +
  #scale_shape_manual(values=1) +
  #xlim(min_lim, max_lim)+ylim (min_lim, max_lim) +
  xlim(min_lim, max_lim)+ylim (-100, 150) +
  geom_polygon(data = hulls,aes(group = hull, fill=color_var), alpha = 0, size = 0.2) +
  theme_stata() +
  #theme(legend.text=element_text(size=6),legend.position = "bottom",legend.title = element_blank()) +
  theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1), legend.title = element_blank(), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  xlab("Distance: Dimension 1") +
  ylab("Distance: Dimension 2") 
print(png)

