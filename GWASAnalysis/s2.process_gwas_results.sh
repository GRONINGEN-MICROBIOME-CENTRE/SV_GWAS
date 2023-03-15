#!/usr/bin/env bash

#
# Process the GWAS results: calculate FDR, annotate, reformat for upload
#

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
resdir=${d}/results/${svtype}/meta_combined_fdr/
script_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SV_GWAS/GWASAnalysis/
genotype_dir=${d}/genotypes/

#
# dSV
#
svtype="dSV"
cut -f1 -d ":" ${d}/data/${svtype}_per_cohort.txt | tail -n+2 | sort | uniq > ${d}/data/${svtype}_species.txt

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
cut -f1 -d ":" ${d}/data/${svtype}_per_cohort.txt | tail -n+2 | sort | uniq > ${d}/data/${svtype}_species.txt

# Combine top results and calculate permutation-based FDR using eqtl mapping pipeline
${script_dir}/gwas_scripts_misc/calculate_FDR.sh $d ${svtype}

# Annotate the results
# TODO: fix how the per cohort z-scores are added
${script_dir}/gwas_scripts_misc/annotate_res.sh ${resdir}/eQTLs.txt.gz ${resdir}/${svtype}_meta_all_5e-08 

# Combine all summary stats per species, reformat for upload
${script_dir}/gwas_scripts_misc/combine_reformat_summary_stats.sh ${svtype}


