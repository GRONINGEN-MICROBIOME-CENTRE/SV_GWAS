d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/

cd ${d}/data_fastGWA

svtype=dSV


#!/usr/bin/env bash

#
# Format and filter the SV tables, prepare all necessary files and scripts, submit GWAS jobs
#

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/
genotype_dir=${d}/genotypes/

#
# 1. Prepare annotation and SV name conversion files
#

# dSVs
bash ${script_dir}/create_sv_name_conversion_tables.sh \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_deletionStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvInfo/20230116_full_Informative_species_information.tsv \
dSV

# vSVs
bash ${script_dir}/create_sv_name_conversion_tables.sh \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_variableStructuralVariation_12388samples.tsv  \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvInfo/20230116_full_Informative_species_information.tsv \
vSV


#
# 2. Format the covariate file: bind species abundances and phenotypes (age, read number)
#

cd ${d}/data_fastGWA/


cohorts=("300OB" "500FG" "LLD" "DAG3")
for cohort in ${cohorts[@]}
do
    
    # Format covariates in plink 1.9, make sure they include all samples with relatives
    sed 's:"::g' abund/${cohort}_abundances.tsv | \
    sed '1s:.:#IID\t&:' | \
    python ${script_dir}/rename_header_based_on_file.py stdin ${svtype}_name_conversion_table.txt 1 3 | \
    python ${script_dir}/add_columns_from_file.py -i stdin  -f pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 | \
    grep -w -v "NA" | \
    awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print }' \
    > ${cohort}.covariates.txt
    
    Rscript ${script_dir}/scale_read_count.R  ${cohort}.covariates.txt ${cohort}.covariates.scaled.txt
    mv ${cohort}.covariates.scaled.txt ${cohort}.covariates.txt
    tail -n+2 ${cohort}.covariates.txt > ${cohort}.covariates.noheader.txt


done

module load R/3.6.1-foss-2018a
cd ${d}

# dSVs:
Rscript ${script_dir}/fastGWA/filter_SV_table.R \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_deletionStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/cleaned_file_list.tsv \
dSV

# vSVs:
Rscript ${script_dir}/fastGWA/filter_SV_table.R \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_variableStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/cleaned_file_list.tsv \
vSV




# make GRMs