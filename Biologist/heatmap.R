#Importing libraries.
library(gplots,quietly = TRUE)
library(dplyr,quietly = TRUE)
library(RColorBrewer,quietly = TRUE)


#Reading Data.
df <- read.csv("deseq_normalized_counts_final.csv",header=TRUE,row.names = 1)

#Calculating coefficient of variation for every gene.
coeff_vec <- sd(as.double(df[1,]))/mean(as.double(df[1,]))

for (x in 2:(nrow(df)))
{
  coeff <- sd(as.double(df[x,]))/mean(as.double(df[x,]))
  coeff_vec <- c(coeff_vec,coeff)
}

df$cv <- coeff_vec
#0.3567 works best

#Filtering genes with low CV.
df_filtered <- filter(df,df$cv>0.3567)

df_filtered <- select(df_filtered,-cv)

#Annotating samples according to groups.
ahr <- c("SRR1177997","SRR1177999","SRR1178002")

car <- c("SRR1178020","SRR1178036","SRR1178046")

cytotoxic <- c("SRR1177987","SRR1177988","SRR1177989")

summary(rowMeans(df_filtered[,ahr], na.rm=TRUE))
summary(rowMeans(df_filtered[,car], na.rm=TRUE))
summary(rowMeans(df_filtered[,cytotoxic], na.rm=TRUE))

#Assigning colors for MoA.
colors <- ""

for(x in colnames(df_filtered)){
  if(x %in% ahr){
    colors <- c(colors,"darkred")
  }
  else if(x %in% car){
    colors <- c(colors,"steelblue")
  }
  else if(x %in% cytotoxic){
    colors <- c(colors,"darkgreen")
  }
  else{
    colors <- c(colors,"hotpink")
  }
}

colors <- colors[2:16]

#Plotting heatmap to visualize hierarchical clustering of samples according to MoAs and gene expression.
theme = brewer.pal(11,"Spectral")

heatmap.2(as.matrix(df_filtered), xlab='Sample', ylab='Gene', 
          main='RNA Seq Gene Expression',
          trace='none', density.info = 'none',ColSideColors = colors,col = theme,
          colCol = colors,key.xlab='Gene Expression Level', scale='row',dendrogram = 'column',labRow = FALSE,margins = c(8,2))

legend("bottomleft",legend = c("AhR","CAR.PXR","Cytotoxic","Control"),xpd=TRUE,cex=1.20,fill=c("darkred","steelblue","darkgreen","hotpink"),title="Toxgroup MoA")

summary(df)

