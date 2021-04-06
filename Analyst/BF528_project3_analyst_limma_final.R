#--------------------------------------
#Title: "BF528 Project 3 Analyst Part 5"
#Author: "Janvee Patel"
#Date: 03/21/2021
#--------------------------------------

#import libraries
library(limma)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

#-----------------------------------------------------------------------------------------------------------------------------------
#Part 5
#differential expression of samples against the controls using the RMA normalized expression matrix

#Chose group_1 microarray file
#Contains information about samples with array_id, chemical, and vehicle columns
samples <- read.csv('/project/bf528/project_3/groups/group_1_mic_info.csv',as.is=TRUE)

#read in the full RMA normalized expression matrix of all experiments
rma <- read.table('/project/bf528/project_3/samples/liver-normalization-rma.txt', sep='\t', as.is=TRUE, header=TRUE, row.names=1)

#subset the full RMA normalized expression matrix to just those array_ids in this comparison (within the samples table)
rma_subset <- rma[paste0('X',samples$array_id)]

#Run limma and get CSV file for the 3 treatment vs control comparisons
gen_limma_results <- function(treatment, vehicle) {
  #subset to get only a particular treatment sample and controls with same vehicle value
  sample_subset <- samples[samples$chemical == treatment | (samples$chemical == "Control" & samples$vehicle == vehicle),]
  sample_subset_ready <- paste0('X', sample_subset$array_id)
  
  #subset rma_subset to include the array_ids for that treatment sample and appropriate controls
  rma_subset_ready <- subset(rma_subset, select = c(sample_subset_ready))
  
  #construct a design matrix modeling treatment vs control for use by limma
  design <- model.matrix(~factor(sample_subset$chemical, levels=c('Control', treatment)))
  colnames(design) <- c('Intercept', treatment)
  
  #run limma
  fit <- lmFit(rma_subset_ready, design)
  fit <- eBayes(fit)
  t <- topTable(fit, coef=2, n=nrow(rma_subset), adjust='BH')
  
  #write results to csv file
  write.csv(t[order(t$adj.P.Val), ], paste0(treatment, "_limma_results.csv"))
}

#3-METHYLCHOLANTHRENE
gen_limma_results("3-METHYLCHOLANTHRENE", "CMC_.5_%")

#CLOTRIMAZOLE
gen_limma_results("CLOTRIMAZOLE", "CORN_OIL_100_%")

#CHLOROFORM
gen_limma_results("CHLOROFORM", "CORN_OIL_100_%")


#report number of genes significant at p-adjust < 0.05 for each analysis
#3-METHYLCHOLANTHRENE
DEG_3ME_m <- read.csv('3-METHYLCHOLANTHRENE_limma_results.csv',as.is=TRUE)
DEG_3ME_m <- DEG_3ME_m %>% filter(adj.P.Val < 0.05) %>% rename(PROBEID = X)

#CLOTRIMAZOLE
DEG_CLO_m <- read.csv('CLOTRIMAZOLE_limma_results.csv',as.is=TRUE)
DEG_CLO_m <- DEG_CLO_m %>% filter(adj.P.Val < 0.05) %>% rename(PROBEID = X)

#CHLOROFORM
DEG_CHL_m <- read.csv('CHLOROFORM_limma_results.csv',as.is=TRUE)
DEG_CHL_m <- DEG_CHL_m %>% filter(adj.P.Val < 0.05) %>% rename(PROBEID = X)


#report top 10 differentially expressed genes by p-value for each analysis
#converted probeid to gene symbol using the refseq_affy_map.csv file and the DAVID database
DEG_3ME_m %>% arrange(P.Value) %>% top_n(10) %>% select(PROBEID)
DEG_CLO_m %>% arrange(P.Value) %>% top_n(10) %>% select(PROBEID)
DEG_CHL_m %>% arrange(P.Value) %>% top_n(10) %>% select(PROBEID)


#generate histograms of fold changes for all 3 treatments
p1 <- ggplot(data=DEG_3ME_m, mapping=aes(x=logFC)) + geom_histogram(bins=50, fill="transparent", color="darkgray") + xlab("Log Fold Change") + ylab("Frequency") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line()) + labs(title="3-Methylcholanthrene", tag="A") + theme(plot.title = element_text(hjust = 0.5))
p2 <- ggplot(data=DEG_CLO_m, mapping=aes(x=logFC)) + geom_histogram(bins=70, fill="transparent", color="darkgray") + xlab("Log Fold Change") + ylab("Frequency") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line()) + labs(title="Clotrimazole", tag="B") + theme(plot.title = element_text(hjust = 0.5))
p3 <- ggplot(data=DEG_CHL_m, mapping=aes(x=logFC)) + geom_histogram(bins=80, fill="transparent", color="darkgray") + xlab("Log Fold Change") + ylab("Frequency") + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line()) + labs(title="Chloroform", tag="C") + theme(plot.title = element_text(hjust = 0.5))
overall <- grid.arrange(p1, p2, p3, ncol=3)
annotate_figure(overall, top=text_grob("Distributions of Log Fold Change for Significant Differentially Expressed Genes\n", size=16), fig.lab.size = 10)


#scatter plots of fold changes vs nominal p-value for all 3 treatments shown as volcano plots
#there were no 0 fold change values found within the DEG files
gen_volcanoplt <- function(treatment, title, tag) {
  #split into not significant
  not_sig <- treatment[treatment$adj.P.Val > 0.05, ]
  nsx = not_sig$logFC
  nsy = -log10(not_sig$P.Value)
  
  #split into significant genes
  sig <- treatment[treatment$adj.P.Val < 0.05, ]
  
  #split into positive fold change
  positive_fc <- sig[sig$logFC > 0, ]
  px = positive_fc$logFC
  py = -log10(positive_fc$P.Value)
  
  #split into negative fold change
  negative_fc <- sig[sig$logFC < 0, ]
  nx = negative_fc$logFC
  ny = -log10(negative_fc$P.Value)
  
  #plot the volcano plot
  clrs <- c("Not Significant"="grey", "Increased"="red", "Decreased"="blue")
  ggplot() + geom_point(data=not_sig, aes(x=nsx, y=nsy, color="Not Significant"), alpha=0.3) + geom_point(data=positive_fc, aes(x=px, y=py, color="Increased"), alpha=0.3) + geom_point(data=negative_fc, aes(x=nx, y=ny, color="Decreased"), alpha=0.3) + xlab("Log Fold Change") + ylab("-log10(P-value)") + theme_bw() + labs(title=title, tag=tag, color="Differentially Expressed") + scale_color_manual(values=clrs) + theme(plot.title = element_text(hjust = 0.5))
}

#the input tables are the full limma results
#3-METHYLCHOLANTHRENE
p1 <- gen_volcanoplt(DEG_3ME_m, "3-Methylcholanthrene", "A")

#CLOTRIMAZOLE
p2 <- gen_volcanoplt(DEG_CLO_m, "Clotrimazole", "B")

#CHLOROFORM
p3 <- gen_volcanoplt(DEG_CHL_m, "Chloroform", "C")

p4 <- grid.arrange(p1, p2, p3, layout_matrix = matrix(c(1,1,2,2,4,3,3,4), nrow = 2, byrow=TRUE))
annotate_figure(p4, top=text_grob("Volcano Plots of Significant Differentially Expressed Genes\n", size=16), fig.lab.size = 10)
