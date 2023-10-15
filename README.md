# Paper

- Title: Host Genetic Regulation of Human Gut Microbial Structural Variation
- Journal: Nature
- Year: 2023

# Summary

This repository contains the code for the manuscript "Host Genetic Regulation of Human Gut Microbial Structural Variation"

# Contents
The code is organised in two folders:

- microbiomeAnalysis contains scripts for the gut microbial SV analysis and downstream microbiome-related analyses
- GWASAnalysis contains scripts for host genetics-related analyses, including association analyses and heritability

## Description of each folder
Here we describe the general workflow, for details please see the comments and descriptions provided in each file.

### microbiomeAnalysis folder

- **s01.cleanData.Rmd**: Clean the raw input data.

- **s02.SV_summary.Rmd**: Summary statistics for gut microbial SV profiles.(Main Fig. 1b and c; Extended Data Fig. 1a-d)

- **s03.Fprau_SV.Rmd**: Analysis of SVs of F. prausnitzii, including the calculation of the populational genetic structure of F. prausnitzii, the correlation between SVs, and associations between SVs and top principal components of the populational genetic structure of F. prausnitzii. (Extended Data Fig. 5)

- **s04.GalNAc_SV.Rmd**: The gene organization and phylogenetic analysis of F. prausnitzii strains used in growth experiments. (Main Fig. 3c; Supplementary Figure 6a-c)

- **s05.GalNAc_gene_search.Rmd**: Summary of homologous of GalNAc utilization genes in genomes of the species previously reported associated with ABO blood type and FUT2 genotype. (Main Fig. 4a-c)

- **s06.GalNAc_gene_assoc.Rmd**: Association analysis of GalNAc gene abundance with gut microbiome diversity/richness and host phenotypes. (Main Fig. 5c-e; Extended Data Fig. 10a and b)

- **s07.Experiment.Rmd**: Visualization of growth curves and qPCR results. (Main Fig. 3d-e; Extended Data Fig. 6; Extended Data Fig. 7)

### GWASAnalysis folder

- **s01.prepare_files_for_fastGWA.sh**: Prepare and filter SV profiles to use in GWAS, and submit GWAS meta-analysis jobs for each SV

- **s02.process_gwas_results.sh**: combine and format the GWAS results

- **s03.run_heritability_analysis.sh**: run family-based heritability estimation in DMP cohort

- **gwas_scripts_misc/ **: this folder contains helper scripts used in s01 and s02 heritability/ : this folder contains helper scripts used in heritability estimation plotting/ : this folder contains scripts used for making plots (manhattan, SV per ABO blood group barplots, etc)

NOTE: Most scripts point to paths on our cluster, so if you want to replicate any of the scripts and need input file format descriptions, please contact us:

- Microbiome-related analysis: Dr. Daoming Wang, wangdaoming94@outlook.com

- Human genetic analysis: Dr. Daria Zhernakova, dashazhernakova@gmail.com



