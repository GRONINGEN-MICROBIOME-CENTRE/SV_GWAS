#!/usr/bin/env bash

#
# Format and filter the SV tables, prepare all necessary files and scripts, submit GWAS jobs
#

d=/data/umcg-tifn/SV/SV_GWAS/
script_dir=/home/umcg-dzhernakova/scripts/umcg_scripts/SV_GWAS/clean/
genotype_dir=/data/umcg-tifn/SV/SV_GWAS/genotypes/

cd ${d}/data/
cohorts=("300OB" "500FG" "LLD" "DAG3")
cohorts_nodag3=("300OB" "500FG" "LLD")


#
# 1. Prepare annotation and SV name conversion files
#

# dSVs
bash ${script_dir}/gwas_scripts_misc/create_sv_name_conversion_tables.sh \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/20211212_full_deletionStructuralVariation_12388samples.tsv \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_info/20211212_full_Informative_species_information.tsv \
dSV

# vSVs
bash ${script_dir}/gwas_scripts_misc/create_sv_name_conversion_tables.sh \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/20211212_full_variableStructuralVariation_12388samples.tsv  \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_info/20211212_full_Informative_species_information.tsv \
vSV


#
# 2. Format the covariate file: bind species abundances, genotype PCs and phenotypes (age, read number)
#
for cohort in ${cohorts_nodag3[@]}
do
    genotypeDir=${genotype_dir}/${cohort}/

    cd ${d}/data
    
    sed 's:"::g' abund/${cohort}_abundances.tsv | \
    sed '1s:.:#IID\t&:' | \
    python ${script_dir}/gwas_scripts_misc/rename_header_based_on_file.py stdin dSV_name_conversion_table.txt 1 3 | \
    python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f ${genotypeDir}/${cohort}.PC1-2.txt -f_m 1 -f_cols 2,3 | \
    python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 | \
    grep -w -v "NA" \
    > ${cohort}.covariates.txt  

    num_l=`wc -l ${cohort}.covariates.txt | awk '{print $1}'` 
    echo "${cohort} has ${num_l} samples with covariates"
done

# Covariate file for DAG3
cohort="DAG3"
genotypeDir=${genotype_dir}/${cohort}/

cd ${d}/data

sed 's:"::g' abund/${cohort}_abundances.tsv | \
sed '1s:.:#IID\t&:' | \
python ${script_dir}/gwas_scripts_misc/rename_header_based_on_file.py stdin dSV_name_conversion_table.txt 1 3 | \
python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f ${genotypeDir}/${cohort}.PC1-5.txt -f_m 1 -f_cols 2,3,4,5,6 | \
python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 \
> ${cohort}.covariates.txt

num_l=`wc -l ${cohort}.covariates.txt | awk '{print $1}'` 
echo "${cohort} has ${num_l} samples with covariates"

#
# 3. filter and rename SVs
#
module load R/3.6.1-foss-2018a

# dSVs:
Rscript ${script_dir}/gwas_scripts_misc/filter_SV_table.R \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/20211212_full_deletionStructuralVariation_12388samples.tsv \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/cleaned_file_list.tsv \
dSV

# vSVs:
Rscript ${script_dir}/gwas_scripts_misc/filter_SV_table.R \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/20211212_full_variableStructuralVariation_12388samples.tsv \
/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/cleaned_file_list.tsv \
vSV


#
# 4. Prepare the GWAS scripts and submit them
#

#dSV
svtype="dSV"
mkdir ${d}/scripts_${svtype}/
cd ${d}/scripts_${svtype}/
cut -f1 ${d}/data/${svtype}_per_cohort.txt | tail -n+2 > all_bacs_${svtype}.txt
rm *split*
split -l$((`wc -l < all_bacs_${svtype}.txt`/40)) all_bacs_${svtype}.txt split. -da 2
for f in split.*
do
    sed "s:__BACLIST__:${f}:g" ${script_dir}//gwas_scripts_misc/run_GWAS_per_dSV_TEMPLATE.sh > run_${f}.sh
    sbatch run_${f}.sh
done

#vSV
svtype="vSV"
mkdir ${d}/scripts_${svtype}/
cd ${d}/scripts_${svtype}/
cut -f1 ${d}/data/${svtype}_per_cohort.txt | tail -n+2 > all_bacs_${svtype}.txt
rm *split*
split -l$((`wc -l < all_bacs_${svtype}.txt`/40)) all_bacs_${svtype}.txt split. -da 2
for f in split.*
do
    sed "s:__BACLIST__:${f}:g" ${script_dir}//gwas_scripts_misc/run_GWAS_per_vSV_TEMPLATE.sh > run_${f}.sh
    sbatch run_${f}.sh
done