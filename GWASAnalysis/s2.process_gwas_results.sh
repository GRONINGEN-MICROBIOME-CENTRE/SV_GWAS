#!/usr/bin/env bash

#
# Process the GWAS results: calculate FDR, annotate, reformat for upload
#

d=/data/umcg-tifn/SV/SV_GWAS/
resdir=${d}/results_all_summary_stats/${svtype}/meta_combined_fdr/
script_dir=/home/umcg-dzhernakova/scripts/umcg_scripts/SV_GWAS/clean/
genotype_dir=/data/umcg-tifn/SV/SV_GWAS/genotypes/

#
# dSV
#
svtype="dSV"
#TODO: make file dSV_species.txt and vSV_species.txt

# Combine top results and calculate permutation-based FDR using eqtl mapping pipeline
${script_dir}/gwas_scripts_misc/calculate_FDR.sh $d ${svtype}

# Annotate the results
# TODO: fix how the per cohort z-scores are added
${script_dir}/gwas_scripts_misc/annotate_res.sh ${resdir}/eQTLs.txt.gz ${resdir}/${svtype}_meta_all_5e-08 

# Combine all summary stats per species, reformat for upload
${script_dir}/gwas_scripts_misc/combine_reformat_summary_stats.sh ${svtype}

#
# vSV
#
svtype="vSV"
#TODO: make file dSV_species.txt and vSV_species.txt

# Combine top results and calculate permutation-based FDR using eqtl mapping pipeline
${script_dir}/gwas_scripts_misc/calculate_FDR.sh $d ${svtype}

# Annotate the results
# TODO: fix how the per cohort z-scores are added
${script_dir}/gwas_scripts_misc/annotate_res.sh ${resdir}/eQTLs.txt.gz ${resdir}/${svtype}_meta_all_5e-08 

# Combine all summary stats per species, reformat for upload
${script_dir}/gwas_scripts_misc/combine_reformat_summary_stats.sh ${svtype}


