# Project Description

A brief description of what this repository is for and what it contains

# Contributors

Zhuorui Sun - Data Curator @sunzhuorui

Camilla Belamarich - Programmer @cmbelama

Janvee Patel - Analyst @Janvee-Patel

Yashrajsinh Jadeja - Biologist @Yashrajsinh-Jadeja

# Repository Contents

Provide a brief description of each script/code file in this repo, what it does, and how to execute it

Programmer:

  - programmer_workflow.txt: This file should be run first. It contains the script used on the command line to run feature counts on each of the samples (including controls). It also includes how multiqc was performed. Lastly, it shows how the data was subset to correctly make the counts matrix csv file.
  - featurecounts_script.Rmd: This file contains the additional subsetting script used to organize the counts matrix in R. Additionally, the script to plot the gene counts boxplot is included. This is an R markdown file.
 - DE_seq_script.Rmd: This file contains the full script used to perform the DESeq2 analysis. Additionally, this script contains the steps used to filter and extract significantly differentailly expressed genes for each treatment group and make a normalized counts file for all sample. Lastly, this script was used to generate histogram and scatter plots. 
  - It's important to metion that I labelled variables with the chemical treatment groups as their mode of action, but in the written report used the chemical name to reference them.
