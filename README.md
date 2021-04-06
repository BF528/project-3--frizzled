# Project Description

Concordance of microarray and RNA-Seq differential gene expression

This study by Wang et. al. aims to determine the concordance between microarrays and RNA-seq differential gene expression by comparing the results from both techniques. The study primarily tests the differences between the two methods by comparing the results from both methods in differential gene expression and creating predictive classifiers. The group devised a comprehensive study design to generate Illumina RNA-seq and Affymetrix microarray data from the same group of liver samples of rats under varying degrees of perturbation by 27 chemicals that represented multiple modes of action (MOA). The cross-platform concordance in terms of differentially expressed genes (DEGs) or enriched pathways was found to be highly correlated with treatment effect size, gene-expression abundance and the biological complexity of the MOA.

In this project, we reproduced the results from Figure 2a and Figure 3b+c, as well as compare the pathway enrichment results reported in the paper. The goals for this project:
 - Align short reads to the rat genome using STAR and quantify expression using a read counting strategy
 - Perform differential expression analysis of RNA-Seq data with DESeq2
 - Perform differential expression analysis of pre-normalized microarray expression data using limma
 - Map between Affymetrix and refSeq identifier systems

Original paper:
 - Wang, C., Gong, B., Bushel, P. R., Thierry-Mieg, J., Thierry-Mieg, D., Xu, J., Fang, H., Hong, H., Shen, J., Su, Z., Meehan, J., Li, X., Yang, L., Li, H., Łabaj, P. P., Kreil, D. P., Megherbi, D., Gaj, S., Caiment, F., van Delft, J., … Tong, W. (2014). The concordance between RNA-seq and microarray data depends on chemical treatment and transcript abundance. Nature biotechnology, 32(9), 926–932. https://doi.org/10.1038/nbt.3001

# Contributors

Zhuorui Sun - Data Curator @sunzhuorui

Camilla Belamarich - Programmer @cmbelama

Janvee Patel - Analyst @Janvee-Patel

Yashrajsinh Jadeja - Biologist @Yashrajsinh-Jadeja

# Repository Contents

Provide a brief description of each script/code file in this repo, what it does, and how to execute it

Data Curator:

  - STAR.qsub : This file used on command line to run STAR on fastq files. Each time the input should be two paired fastq files with same sample name. (For example we have XXX_1.fastq.gz and XXX_2.fastq.gz, input ```qsub STAR.qsub XXX``` to run)
  - Multiqc.qsub: This file used on command line to run MultiQC based on the FastQC and STAR results we generated before.
  - Multiqc_report_1.html: The report by Multiqc.
  
Programmer:

  - programmer_workflow.txt: This file should be run first. It contains the script used on the command line to run feature counts on each of the samples (including controls). It also includes how multiqc was performed. Lastly, it shows how the data was subset to correctly make the counts matrix csv file.
  - featurecounts_script.Rmd: This file contains the additional subsetting script used to organize the counts matrix in R. Additionally, the script to plot the gene counts boxplot is included. This is an R markdown file.
 - DE_seq_script.Rmd: This file contains the full script used to perform the DESeq2 analysis. Additionally, this script contains the steps used to filter and extract significantly differentailly expressed genes for each treatment group and make a normalized counts file for all sample. Lastly, this script was used to generate histogram and scatter plots. 
  - It's important to metion that I labelled variables with the chemical treatment groups as their mode of action, but in the written report used the chemical name to reference them.

Analyst:
  - BF528_project3_analyst_limma_final.R: This file should be run first. This contains the script used to perform Limma analysis for the microarray data. The total and top 10 significant differentially expressed genes were determined for each treatment group. 
  - BF528_project3_analyst_final.R: This file should be run after all the above Data Curator, Programmer, and Analyst files have been executed. This contains a script to perform concordance analysis of DEG between the RNA-Seq and Microarray platforms for each of the 3 treatments: 3-Methylcholanthrene, Clotrimazole, and Chloroform. In addition, above and below-median subsets have been determined using the Average Expression of Microarray and baseMean of RNA-Seq, and concordance was computed on these subsets as well. Note that this file contains additional code at the end which used a different filtering method for determining DEG; the data corresponding to this section is shown in Supplemental Table 1 within the Project 3 Final Write-Up.  

Biologist:
