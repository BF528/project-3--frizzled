#Importing packages.
library(dplyr)

#Reading data.
AhR <- read.csv("AhR.csv",header=TRUE,row.names=1)

#Applying filters.
AhR_filtered <- filter(AhR,AhR$adj.P.Val<0.05 & abs(AhR$logFC)>log2(1.5))

CAR_PXR <- read.csv("CAR_PXR.csv",header=TRUE,row.names=1)

CAR_PXR_filtered <- filter(CAR_PXR,CAR_PXR$adj.P.Val<0.05 & abs(CAR_PXR$logFC)>log2(1.5))

Cytotoxic <- read.csv("Cytotoxic.csv",header=TRUE,row.names = 1)

Cytotoxic_filtered <- filter(Cytotoxic,Cytotoxic$adj.P.Val<0.05 & abs(Cytotoxic$logFC)>log2(1.5))

#Writing csv files.
write.table(rownames(AhR_filtered),file="AhR_significant",row.names = FALSE,col.names=FALSE,quote=FALSE)
write.table(rownames(CAR_PXR_filtered),file="CAR_PXR_significant",row.names = FALSE,col.names=FALSE,quote=FALSE)
write.table(rownames(Cytotoxic_filtered)[1:3000],file="Cytotoxic_significant",row.names = FALSE,col.names=FALSE,quote=FALSE)


AhR_deseq <- read.csv("significant_AhR_final.csv",header=TRUE)

AhR_filtered_deseq <- filter(AhR_deseq,AhR_deseq$AhR_padj<0.05 & abs(AhR_deseq$AhR_log2FoldChange)>log2(1.5))

CAR_PXR_deseq <- read.csv("significant_CAR.PXR_final.csv",header=TRUE)

CAR_PXR_filtered_deseq <- filter(CAR_PXR_deseq,CAR_PXR_deseq$CAR.PXR_padj<0.05 & abs(CAR_PXR_deseq$CAR.PXR_log2FoldChange)>log2(1.5))

Cytotoxic_deseq <- read.csv("significant_Cytotoxic_final.csv",header=TRUE)

Cytotoxic_filtered_deseq <- filter(Cytotoxic_deseq,Cytotoxic_deseq$Cytotoxic_padj<0.05 & abs(Cytotoxic_deseq$Cytotoxic_log2FoldChange)>log2(1.5))

write.table(AhR_filtered_deseq[,1],file="AhR_significant_deseq",row.names = FALSE,col.names=FALSE,quote=FALSE)
write.table(CAR_PXR_filtered_deseq[,1],file="CAR_PXR_significant_deseq",row.names = FALSE,col.names=FALSE,quote=FALSE)
write.table(Cytotoxic_filtered_deseq[,1],file="Cytotoxic_significant_deseq",row.names = FALSE,col.names=FALSE,quote=FALSE)
