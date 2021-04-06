#--------------------------------------------
#Title: "BF528 Project 3 Analyst Concordance"
#Author: "Janvee Patel"
#Date: 03/21/2021
#--------------------------------------------

#import libraries
library(limma)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

#-----------------------------------------------------------------------------------------------------------------------------------
#Part 6
#Concordance method from online methods in paper

#Map Affymetrix probe ids to RefSeq identifiers
#read in the mapping file for probe ids and refseq ids
idmap <- read.csv("/project/bf528/project_3/refseq_affy_map.csv")[,c(1,2,3)]

#remove rows with NA values in the probeid and symbol columns combined
idmap <- na.omit(idmap)

#identify duplicates within the id map
#duplct <- idmap %>% group_by(SYMBOL) %>% dplyr::filter(n() > 1)

#-----------------------------------------------------------------------------------------------------------------------------------
#Tried 2 different ways of filtering by Adjusted P-value < 0.05 and P-value < 0.05
#Method 1 - Filtered by Adjusted P-value < 0.05 (taken) without a log FC cutoff
#Did not use log FC cutoff as there were quite less total number of DEG for 3-Methylcholanthrene

#read in the limma results from the Microarray analysis for the 3 treatment analysis
#3-METHYLCHLOANTHRENE
DEG_3ME_m1 <- read.csv('3-METHYLCHOLANTHRENE_limma_results.csv',as.is=TRUE) %>% rename(PROBEID = X) %>% filter(adj.P.Val < 0.05)
#CLOTRIMAZOLE
DEG_CLO_m1 <- read.csv('CLOTRIMAZOLE_limma_results.csv',as.is=TRUE) %>% rename(PROBEID = X) %>% filter(adj.P.Val < 0.05)
#CHLOROFORM
DEG_CHL_m1 <- read.csv('CHLOROFORM_limma_results.csv',as.is=TRUE) %>% rename(PROBEID = X) %>% filter(adj.P.Val < 0.05)

#read in the DESeq results from the RNA-Seq analysis for the 3 treatment analyses
#3-METHYLCHOLANTHRENE
DEG_3ME_r1 <- read.csv("/projectnb/bf528/users/frizzled/project_3/Programmer/results/full_AhR_DESeq2.csv")  %>% rename(REFSEQ = X) %>% filter(padj < 0.05)
#CLOTRIMAZOLE
DEG_CLO_r1 <- read.csv("/projectnb/bf528/users/frizzled/project_3/Programmer/results/full_CAR.PXR_DESeq2.csv")  %>% rename(REFSEQ = X) %>% filter(padj < 0.05)
#CHLOROFORM
DEG_CHL_r1 <- read.csv("/projectnb/bf528/users/frizzled/project_3/Programmer/results/full_Cytotoxic_DESeq2.csv")  %>% rename(REFSEQ = X) %>% filter(padj < 0.05)

#-----------------------------------------------------------------------------------------------------------------------------------

#Concordance method from the online methods section of the Wang et al. paper
#2 x intersect(DEGsmicroarray, DEGsrnaseq) / (DEGsmicroarray + DEGsrnaseq) with agreement in direction of fold change

#N = number of items in whole sets
#n1 and n2 = number of items in individual sets
#n0 = number of items in observed intersection
#nx = background corrected intersection

#Tried 2 different ways of mapping the probe ids to refseq ids
#Both methods gave the same number of total genes
#Did not use FC > 1.5 cutoff for the adjusted P-value < 0.05 as the number of total genes drastically decreased, and there weren't enough genes for concordance

#-----------------------------------------------------------------------------------------------------------------------------------
#Overall

#3-METHYLCHOLANTHRENE
#Method 1 - used to validate number of DEG determined in the second method
mapped3ME <- idmap[which(idmap$PROBEID %in% DEG_3ME_m1$PROBEID & idmap$REFSEQ %in% DEG_3ME_r1$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, DEG_3ME_m1, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, DEG_3ME_r1, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
#researched the number of coding genes in the rat genome (Rnor_6.0) which I found to be approximately 22,250 (used an approximation, not an exact value)
N = 22250  #number of approximate coding genes in rat genome
n1 = nrow(DEG_3ME_m1) #independent set of significant DEG from microarray
n2 = nrow(DEG_3ME_r1)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
c1 = (2*(nx))/(n1+n2)*100 #concordance calculation


#CLOTRIMAZOLE
#Method 1 - used to validate number of DEG determined in the second method
mappedCLO <- idmap[which(idmap$PROBEID %in% DEG_CLO_m1$PROBEID & idmap$REFSEQ %in% DEG_CLO_r1$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, DEG_CLO_m1, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, DEG_CLO_r1, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change
#4 mappings did not agree in direction of fold change - these were omitted for concordance analysis

#concordance calculation
N = 22250  #number of approximate coding genes in the rat genome
n1 = nrow(DEG_CLO_m1) #independent set of significant DEG from microarray
n2 = nrow(DEG_CLO_r1)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
c2 = (2*(nx))/(n1+n2)*100  #concordance calculation


#CHLOROFORM
#Method 1 - used to validate number of DEG determined in the second method
mappedCHL <- idmap[which(idmap$PROBEID %in% DEG_CHL_m1$PROBEID & idmap$REFSEQ %in% DEG_CHL_r1$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, DEG_CHL_m1, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, DEG_CHL_r1, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change
#17 mappings did not agree in direction of fold change - these were omitted for concordance analysis

#Concordance calculation
N = 22250 #number of approximate genes in rat genome
n1 = nrow(DEG_CHL_m1) #independent set of significant DEG from microarray
n2 = nrow(DEG_CHL_r1)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
c3 = ((2*(nx))/(n1+n2))*100  #concordance calculation

#-----------------------------------------------------------------------------------------------------------------------------------
#generate plots of concordance vs number of DEG from RNA-Seq Analysis
r <- data.frame(DEG=c(nrow(DEG_3ME_r1), nrow(DEG_CLO_r1), nrow(DEG_CHL_r1)), Concordance=c(c1, c2, c3))
Treatments = c("3-Methylcholanthrene", "Clotrimazole", "Chloroform")
clrs <- colorRampPalette(brewer.pal(3, "Dark2"))(3)
g1 <- ggplot(r, aes(DEG, Concordance, color=Treatments)) + geom_point(size=3) + geom_smooth(method="lm", se=FALSE, color="darkgrey", linetype="dashed", formula=y~x) + stat_regline_equation(label.y = 25, aes(label = ..rr.label..), color="black", size=5) + xlab("Number of DEG from RNA-Seq") + ylab("Concordance of DEG") + annotate("text", x = 313, y = 10, label = "3ME") + annotate("text", x = 930, y = 30, label = "CLO") + annotate("text", x = 1728, y = 29, label = "CHL") + theme_bw() +labs(title="RNA-Seq") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=clrs)

#generate plots of concordance vs number of DEG from RNA-Seq Analysis
m <- data.frame(DEG=c(nrow(DEG_3ME_m1), nrow(DEG_CLO_m1), nrow(DEG_CHL_m1)), Concordance=c(c1, c2, c3))
Treatments = c("3-Methylcholanthrene", "Clotrimazole", "Chloroform")
clrs <- colorRampPalette(brewer.pal(3, "Dark2"))(3)
g2 <- ggplot(m, aes(DEG, Concordance, color=Treatments)) + geom_point(size=3) + geom_smooth(method="lm", se=FALSE, color="darkgrey", linetype="dashed", formula=y~x) + stat_regline_equation(label.y = 25, aes(label = ..rr.label..), color="black", size=5) + xlab("Number of DEG from Microarray") + ylab("Concordance of DEG") + annotate("text", x = 58, y = 10, label = "3ME") + annotate("text", x = 2692, y = 30, label = "CLO") + annotate("text", x = 9458, y = 29, label = "CHL") + theme_bw() +labs(title="Microarray") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values=clrs)

overall <- ggarrange(g1, g2, ncol=2, common.legend=TRUE, legend="top")
annotate_figure(overall, top=text_grob("Overall Concordance Versus Number of Significant Differentially Expressed Genes\n", size=16), fig.lab.size = 10)


#-----------------------------------------------------------------------------------------------------------------------------------
#Split individuals sets into above-median and below-median groups

#Above-Median
#split the individual sets by the median of the AveExpr/baseMean columns - this will be n1 and n2
#then combine the above and below groups separately and find how many overlap - this will be n0
#above groups contain the median if they were an odd number set

#3-METHYLCHOLANTHRENE
#Microarray DEG
#compute median
medm <- median(DEG_3ME_m1$AveExpr)
mabove <- subset(DEG_3ME_m1, DEG_3ME_m1$AveExpr >= medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_3ME_r1$baseMean)
rabove <- subset(DEG_3ME_r1, DEG_3ME_r1$baseMean >= medr)  #n2

#Method 1 - used to validate number of DEG determined in method 2
mapped3ME <- idmap[which(idmap$PROBEID %in% mabove$PROBEID & idmap$REFSEQ %in% rabove$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mabove, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rabove, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
#researched the number of genes in the rat genome which was determined to be approximately 25,000 (used an approximation, not an exact value)
N=22250  #number of coding genes in rat genome
n1 = nrow(mabove) #independent set of above-median significant DEG from microarray
n2 = nrow(rabove)  #independent set of above-median significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
a1 = (2*(nx))/(n1+n2)*100 #concordance calculation


#CLOTRIMAZOLE
#Microarray DEG
#compute median
medm <- median(DEG_CLO_m1$AveExpr)
mabove <- subset(DEG_CLO_m1, DEG_CLO_m1$AveExpr >= medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CLO_r1$baseMean)
rabove <- subset(DEG_CLO_r1, DEG_CLO_r1$baseMean >= medr)  #n2

#Method 1
mappedCLO <- idmap[which(idmap$PROBEID %in% mabove$PROBEID & idmap$REFSEQ %in% rabove$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mabove, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rabove, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
N = 22250  #number of coding genes in the rat genome
n1 = nrow(mabove) #independent set of significant DEG from microarray
n2 = nrow(rabove)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
a2 = ((2*(nx))/(n1+n2))*100  #concordance calculation

#CHLOROFORM
#Microarray DEG
#compute median
medm <- median(DEG_CHL_m1$AveExpr)
mabove <- subset(DEG_CHL_m1, DEG_CHL_m1$AveExpr >= medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CHL_r1$baseMean)
rabove <- subset(DEG_CHL_r1, DEG_CHL_r1$baseMean >= medr)  #n2

#Method 1
mappedCHL <- idmap[which(idmap$PROBEID %in% mabove$PROBEID & idmap$REFSEQ %in% rabove$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mabove, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rabove, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mabove) #independent set of significant DEG from microarray
n2 = nrow(rabove)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
a3 = ((2*(nx))/(n1+n2))*100  #concordance calculation


#-----------------------------------------------------------------------------------------------------------------------------------
#Below-Median

#3-METHYLCHOLANTHRENE
#Microarray DEG
#compute median
medm <- median(DEG_3ME_m1$AveExpr)
mbelow <- subset(DEG_3ME_m1, DEG_3ME_m1$AveExpr < medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_3ME_r1$baseMean)
rbelow <- subset(DEG_3ME_r1, DEG_3ME_r1$baseMean < medr)  #n2

#Method 1
mapped3ME <- idmap[which(idmap$PROBEID %in% mbelow$PROBEID & idmap$REFSEQ %in% rbelow$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mbelow, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rbelow, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
#researched the number of genes in the rat genome which was determined to be approximately 25,000 (used an approximation, not an exact value)
N=22250  #number of coding genes in rat genome
n1 = nrow(mbelow) #independent set of above-median significant DEG from microarray
n2 = nrow(rbelow)  #independent set of above-median significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
b1 = (2*(nx))/(n1+n2)*100 #concordance calculation


#CLOTRIMAZOLE
#Microarray DEG
#compute median
medm <- median(DEG_CLO_m1$AveExpr)
mbelow <- subset(DEG_CLO_m1, DEG_CLO_m1$AveExpr < medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CLO_r1$baseMean)
rbelow <- subset(DEG_CLO_r1, DEG_CLO_r1$baseMean < medr)  #n2

#Method 1 (not taken)
mappedCLO <- idmap[which(idmap$PROBEID %in% mbelow$PROBEID & idmap$REFSEQ %in% rbelow$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mbelow, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rbelow, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mbelow) #independent set of significant DEG from microarray
n2 = nrow(rbelow)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
b2 = ((2*(nx))/(n1+n2))*100  #concordance calculation

#CHLOROFORM
#Microarray DEG
#compute median
medm <- median(DEG_CHL_m1$AveExpr)
mbelow <- subset(DEG_CHL_m1, DEG_CHL_m1$AveExpr < medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CHL_r1$baseMean)
rbelow <- subset(DEG_CHL_r1, DEG_CHL_r1$baseMean < medr)  #n2

#Method 1 (not taken)
mappedCHL <- idmap[which(idmap$PROBEID %in% mbelow$PROBEID & idmap$REFSEQ %in% rbelow$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mbelow, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rbelow, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mbelow) #independent set of significant DEG from microarray
n2 = nrow(rbelow)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
b3 = ((2*(nx))/(n1+n2))*100  #concordance calculation


#------------------------------------------------------------------------------------------------------------------------------------
#Plot combined bar plot with overall, above-median, and below-median concordance
treatment <- c(rep("3-Methylcholanthrene", 3), rep("Clotrimazole", 3), rep("Chloroform", 3)) #treatments
Groups <- rep(c("Overall", "Above-Median", "Below-Median"), 3)
concordance_vals <- c(c1, a1, b1, c2, a2, b2, c3, a3, b3)  #contain concordance values
df <- data.frame(treatment, groups, concordance_vals)
Groups <- factor(Groups, levels=c("Overall", "Above-Median", "Below-Median"))
df$treatment <- factor(df$treatment, levels=c("3-Methylcholanthrene", "Clotrimazole", "Chloroform"))

#plot barplot
clrs <- colorRampPalette(brewer.pal(3, "Dark2"))(3)
ggplot(df, aes(fill=Groups, x=treatment, y=concordance_vals)) + geom_bar(position="dodge", stat="identity") + scale_fill_brewer(palette = "Dark2") + xlab("Treatments") + ylab("Concordance") + labs(title="Concordance of DEG in Overall, Above-Median, and Below-Median Sets") + theme_bw() + theme(plot.title = element_text(size = 16)) + theme(plot.title = element_text(hjust = 0.5))




#------------------------------------------------------------------------------------------------------------------------------------
#Additional Code for Method 2
#Computing overall concordance using p-value < 0.05 and FC > 1.5
#Data from this method is included in the supplemental section

#read in the limma results from the Microarray analysis for the 3 treatment analysis
#3-METHYLCHLOANTHRENE
DEG_3ME_m2 <- read.csv('3-METHYLCHOLANTHRENE_limma_results.csv',as.is=TRUE) %>% rename(PROBEID = X) %>% filter(P.Value < 0.05) %>% filter(abs(logFC) > 0.58)
#CLOTRIMAZOLE
DEG_CLO_m2 <- read.csv('CLOTRIMAZOLE_limma_results.csv',as.is=TRUE) %>% rename(PROBEID = X) %>% filter(P.Value < 0.05) %>% filter(abs(logFC) > 0.58)
#CHLOROFORM
DEG_CHL_m2 <- read.csv('CHLOROFORM_limma_results.csv',as.is=TRUE) %>% rename(PROBEID = X) %>% filter(P.Value < 0.05) %>% filter(abs(logFC) > 0.58)

#read in the DESeq results from the RNA-Seq analysis for the 3 treatment analyses
#3-METHYLCHOLANTHRENE
DEG_3ME_r2 <- read.csv("/projectnb/bf528/users/frizzled/project_3/Programmer/results/full_AhR_DESeq2.csv")  %>% rename(REFSEQ = X) %>% filter(pvalue < 0.05) %>% filter(abs(log2FoldChange) > 0.58)
#CLOTRIMAZOLE
DEG_CLO_r2 <- read.csv("/projectnb/bf528/users/frizzled/project_3/Programmer/results/full_CAR.PXR_DESeq2.csv")  %>% rename(REFSEQ = X) %>% filter(pvalue < 0.05) %>% filter(abs(log2FoldChange) > 0.58)
#CHLOROFORM
DEG_CHL_r2 <- read.csv("/projectnb/bf528/users/frizzled/project_3/Programmer/results/full_Cytotoxic_DESeq2.csv")  %>% rename(REFSEQ = X) %>% filter(pvalue < 0.05) %>% filter(abs(log2FoldChange) > 0.58)



#Overall
#3-METHYLCHOLANTHRENE
#Method 1 - used to validate number of DEG determined in the second method
mapped3ME <- idmap[which(idmap$PROBEID %in% DEG_3ME_m2$PROBEID & idmap$REFSEQ %in% DEG_3ME_r2$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, DEG_3ME_m2, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, DEG_3ME_r2, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
#researched the number of genes in the rat genome which I found to be approximately 25,000 (used an approximation, not an exact value)
N = 22250  #number of approximate genes in rat genome
n1 = nrow(DEG_3ME_m2) #independent set of significant DEG from microarray
n2 = nrow(DEG_3ME_r2)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
c1 = (2*(nx))/(n1+n2)*100 #concordance calculation


#CLOTRIMAZOLE
#Method 1 - used to validate number of DEG determined in the second method
mappedCLO <- idmap[which(idmap$PROBEID %in% DEG_CLO_m2$PROBEID & idmap$REFSEQ %in% DEG_CLO_r2$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, DEG_CLO_m2, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, DEG_CLO_r2, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change
#1 mapping did not agree in direction of fold change - this was omitted for concordance analysis

#concordance calculation
N = 22250  #number of approximate genes in the rat genome
n1 = nrow(DEG_CLO_m2) #independent set of significant DEG from microarray
n2 = nrow(DEG_CLO_r2)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
c2 = (2*(nx))/(n1+n2)*100  #concordance calculation


#CHLOROFORM
#Method 1 - used to validate number of DEG determined in the second method
mappedCHL <- idmap[which(idmap$PROBEID %in% DEG_CHL_m2$PROBEID & idmap$REFSEQ %in% DEG_CHL_r2$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, DEG_CHL_m2, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, DEG_CHL_r2, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change
#7 mappings did not agree in direction of fold change - these were omitted for concordance analysis

#Concordance calculation
N = 22250 #number of approximate genes in rat genome
n1 = nrow(DEG_CHL_m2) #independent set of significant DEG from microarray
n2 = nrow(DEG_CHL_r2)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
c3 = ((2*(nx))/(n1+n2))*100  #concordance calculation



#Above-Median
#3-METHYLCHOLANTHRENE
#Microarray DEG
#compute median
medm <- median(DEG_3ME_m2$AveExpr)
mabove <- subset(DEG_3ME_m2, DEG_3ME_m2$AveExpr >= medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_3ME_r2$baseMean)
rabove <- subset(DEG_3ME_r2, DEG_3ME_r2$baseMean >= medr)  #n2

#Method 1 (not taken)
mapped3ME <- idmap[which(idmap$PROBEID %in% mabove$PROBEID & idmap$REFSEQ %in% rabove[,1]), ]

#Method 2 (taken)
comf <- merge(idmap, mabove, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rabove, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
#researched the number of genes in the rat genome which was determined to be approximately 25,000 (used an approximation, not an exact value)
N=22250  #number of genes in rat genome
n1 = nrow(mabove) #independent set of above-median significant DEG from microarray
n2 = nrow(rabove)  #independent set of above-median significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
a1 = (2*(nx))/(n1+n2)*100 #concordance calculation


#CLOTRIMAZOLE
#Microarray DEG
#compute median
medm <- median(DEG_CLO_m2$AveExpr)
mabove <- subset(DEG_CLO_m2, DEG_CLO_m2$AveExpr >= medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CLO_r2$baseMean)
rabove <- subset(DEG_CLO_r2, DEG_CLO_r2$baseMean >= medr)  #n2

#Method 1 (not taken)
mappedCLO <- idmap[which(idmap$PROBEID %in% mabove$PROBEID & idmap$REFSEQ %in% rabove$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mabove, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rabove, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change
#2 mappings did not agree in direction of fold change - these were omitted from concordance analysis

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mabove) #independent set of significant DEG from microarray
n2 = nrow(rabove)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
a2 = ((2*(nx))/(n1+n2))*100  #concordance calculation

#CHLOROFORM
#Microarray DEG
#compute median
medm <- median(DEG_CHL_m2$AveExpr)
mabove <- subset(DEG_CHL_m2, DEG_CHL_m2$AveExpr >= medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CHL_r2$baseMean)
rabove <- subset(DEG_CHL_r2, DEG_CHL_r2$baseMean >= medr)  #n2

#Method 1 - used to validate number of DEG determined in second method
mappedCHL <- idmap[which(idmap$PROBEID %in% mabove$PROBEID & idmap$REFSEQ %in% rabove$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mabove, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rabove, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change
#5 mappings that did not agree in direction of fold change - these were omitted from concordance analysis

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mabove) #independent set of significant DEG from microarray
n2 = nrow(rabove)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
a3 = ((2*(nx))/(n1+n2))*100  #concordance calculation



#Below-Median
#3-METHYLCHOLANTHRENE
#Microarray DEG
#compute median
medm <- median(DEG_3ME_m2$AveExpr)
mbelow <- subset(DEG_3ME_m2, DEG_3ME_m2$AveExpr < medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_3ME_r2$baseMean)
rbelow <- subset(DEG_3ME_r2, DEG_3ME_r2$baseMean < medr)  #n2

#Method 1 (not taken)
mapped3ME <- idmap[which(idmap$PROBEID %in% mbelow$PROBEID & idmap$REFSEQ %in% rbelow$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mbelow, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rbelow, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
#researched the number of genes in the rat genome which was determined to be approximately 25,000 (used an approximation, not an exact value)
N=22250  #number of genes in rat genome
n1 = nrow(mbelow) #independent set of above-median significant DEG from microarray
n2 = nrow(rbelow)  #independent set of above-median significant DEG from RNA-Seq
n0 = nrow(comf)  #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
b1 = (2*(nx))/(n1+n2)*100 #concordance calculation


#CLOTRIMAZOLE
#Microarray DEG
#compute median
medm <- median(DEG_CLO_m2$AveExpr)
mbelow <- subset(DEG_CLO_m2, DEG_CLO_m2$AveExpr < medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CLO_r2$baseMean)
rbelow <- subset(DEG_CLO_r2, DEG_CLO_r2$baseMean < medr)  #n2

#Method 1 (not taken)
mappedCLO <- idmap[which(idmap$PROBEID %in% mbelow$PROBEID & idmap$REFSEQ %in% rbelow$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mbelow, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rbelow, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mbelow) #independent set of significant DEG from microarray
n2 = nrow(rbelow)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
b2 = ((2*(nx))/(n1+n2))*100  #concordance calculation

#CHLOROFORM
#Microarray DEG
#compute median
medm <- median(DEG_CHL_m2$AveExpr)
mbelow <- subset(DEG_CHL_m2, DEG_CHL_m2$AveExpr < medm)  #n1

#RNA-Seq DEG
#compute median
medr <- median(DEG_CHL_r2$baseMean)
rbelow <- subset(DEG_CHL_r2, DEG_CHL_r2$baseMean < medr)  #n2

#Method 1 (not taken)
mappedCHL <- idmap[which(idmap$PROBEID %in% mbelow$PROBEID & idmap$REFSEQ %in% rbelow$REFSEQ), ]

#Method 2 (taken)
comf <- merge(idmap, mbelow, by="PROBEID")  #merge the idmap with the DEG from microarray by probeid
comf <- merge(comf, rbelow, by="REFSEQ")  #merge the above with the DEG from rna-seq by refseq
comf <- comf[, c(3,2,4,5,7,8,1,10,11,14,15)]
comf$direction_m <- ifelse(comf$logFC >= 0, "+", "-")  #added column for the direction of fold change from the microarray analysis
comf$direction_r <- ifelse(comf$log2FoldChange >= 0, "+", "-")  #added column for the direction of fold change from the rna-seq analysis
comf <- comf[comf$direction_m==comf$direction_r, ]  #only include rows that match in direction of fold change

#concordance calculation
N = 22250  #number of genes in the rat genome
n1 = nrow(mbelow) #independent set of significant DEG from microarray
n2 = nrow(rbelow)  #independent set of significant DEG from RNA-Seq
n0 = nrow(comf) #number of observed intersection
nx <- ((n0*N) - (n1*n2))/(n0+N-n1-n2)  #background corrected intersection
b3 = ((2*(nx))/(n1+n2))*100  #concordance calculation

