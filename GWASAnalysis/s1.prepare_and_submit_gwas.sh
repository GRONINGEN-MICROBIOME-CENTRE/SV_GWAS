#!/usr/bin/env bash

#
# Format and filter the SV tables, prepare all necessary files and scripts, submit GWAS jobs
#

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
script_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SV_GWAS/GWASAnalysis/
genotype_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/

cd ${d}/data/
cohorts=("300OB" "500FG" "LLD" "DAG3")
cohorts_nodag3=("300OB" "500FG" "LLD")


#
# 1. Prepare annotation and SV name conversion files
#

# dSVs
bash ${script_dir}/gwas_scripts_misc/create_sv_name_conversion_tables.sh \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_deletionStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvInfo/20230116_full_Informative_species_information.tsv \
dSV

# vSVs
bash ${script_dir}/gwas_scripts_misc/create_sv_name_conversion_tables.sh \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_variableStructuralVariation_12388samples.tsv  \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvInfo/20230116_full_Informative_species_information.tsv \
vSV


#
# 2. Format the covariate file: bind species abundances, genotype PCs and phenotypes (age, read number)
#

cd ${d}/data/

for cohort in ${cohorts[@]}
do
    genotypeDir=${genotype_dir}/${cohort}/
    sed -i "s:\tIID:\t#IID:g" ${genotypeDir}/${cohort}.PC1-2.txt
    
    # 2.1 Gather age and read count for each cohort
    ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i sample_ids/${cohort}_samples_with_SVs.txt -f pheno/${cohort}_age.txt -f_cols 1 | cut -f1,4 | ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin -f /groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/readPairNum/${cohort}.read_pair_number.tsv -f_cols 1 | grep -v "NA" | sed 1i"#IID\tage\tread_number" > pheno/${cohort}_pheno.txt

    # 2.2 Check the number of samples with SVs, phenoypes and genotypes
    echo "Counting samples for ${cohort}."
    echo "N samples with SVs --- of them N with genotypes --- of them N with phenotypes --- of them N with both phenotypes and genotypes"
    
    python3 ${script_dir}/gwas_scripts_misc/check_sample_overlap.py sample_ids/${cohort}_samples_with_SVs.txt 0  ${genotype_dir}/${cohort}/${cohort}_filtered.fam 1 pheno/${cohort}_pheno.txt 0 "id\twith_geno\twith_pheno" > sample_ids/${cohort}_id_overlap.txt

    # 2.3 Create a covariate file: species abundances, age, read count and genotype PCs
    
    if [ $cohort == "DAG3" ]
    then
    
        sed 's:"::g' abund/${cohort}_abundances.tsv | \
        sed '1s:.:#IID\t&:' | \
        python ${script_dir}/gwas_scripts_misc/rename_header_based_on_file.py stdin dSV_name_conversion_table.txt 1 3 | \
        python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f ${genotypeDir}/${cohort}.PC1-5.txt -f_m 1 -f_cols 2,3,4,5,6 | \
        python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 | \
        grep -w -v "NA" \
        > ${cohort}.covariates.txt
            
    else
        sed 's:"::g' abund/${cohort}_abundances.tsv | \
        sed '1s:.:#IID\t&:' | \
        python ${script_dir}/gwas_scripts_misc/rename_header_based_on_file.py stdin dSV_name_conversion_table.txt 1 3 | \
        python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f ${genotypeDir}/${cohort}.PC1-2.txt -f_m 1 -f_cols 2,3 | \
        python ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin  -f pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 | \
        grep -w -v "NA" \
        > ${cohort}.covariates.txt  
    fi
    num_l=`tail -n+2 ${cohort}.covariates.txt | wc -l | awk '{print $1}'`
    echo "${cohort} has ${num_l} samples with covariates"
done


#
# 3. filter and rename SVs
#
module load R/3.6.1-foss-2018a
cd ${d}

# dSVs:
Rscript ${script_dir}/gwas_scripts_misc/filter_SV_table.R \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_deletionStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/cleaned_file_list.tsv \
dSV

# vSVs:
Rscript ${script_dir}/gwas_scripts_misc/filter_SV_table.R \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_variableStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/cleaned_file_list.tsv \
vSV


#
# 4. Prepare the GWAS scripts and submit them
#

#dSV
svtype="dSV"
mkdir ${d}/scripts/scripts_${svtype}/
cd ${d}/scripts/scripts_${svtype}/
mkdir logs

cut -f1 ${d}/data/${svtype}_per_cohort.txt | tail -n+2 > all_bacs_${svtype}.txt

while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/run_${sv}.out \
    -e logs/run_${sv}.err \
    -J dsv_${sv} \
    ${script_dir}//gwas_scripts_misc/run_GWAS_per_dSV_TEMPLATE.sh "${sv}"
done < all_bacs_${svtype}.txt


# Check which SVs are missing
result_dir=${d}/results/${svtype}/meta/
result_dir="/groups/umcg-fu/tmp01/projects/SV_GWAS/GWAS_tmp//results/dSV/meta/"
while read line
do
    sv=$line
    sv_resdir=${result_dir}/${sv}/
    if [ ! -f "${sv_resdir}/${sv}.meta_res.annot.tbl.gz"  ] ||  [ ! -f "${sv_resdir}/${sv}.meta_res.eQTLs.txt.gz" ]
    then
        echo -e "$sv\t0"
    fi
    for i in `seq 1 10`
    do
       if [ ! -f "${sv_resdir}/${sv}.meta_res.eQTLs.perm${i}.txt.gz"  ]
        then
            echo -e "$sv\t${i}"
        fi 
    done
    
done < all_bacs_${svtype}.txt


#vSV
svtype="vSV"
mkdir ${d}/scripts/scripts_${svtype}/
cd ${d}/scripts/scripts_${svtype}/
mkdir logs

cut -f1 ${d}/data/${svtype}_per_cohort.txt | tail -n+2 > all_bacs_${svtype}.txt
while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/run_${sv}.out \
    -e logs/run_${sv}.err \
    -J dsv_${sv} \
    ${script_dir}//gwas_scripts_misc/run_GWAS_per_dSV_TEMPLATE.sh "${sv}"
done < all_bacs_${svtype}.txt

# Check which SVs are missing
result_dir=${d}/results/${svtype}/meta/
while read line
do
    sv=$line
    sv_resdir=${result_dir}/${sv}/
    all_finished=1
    if [ ! -f "${sv_resdir}/${sv}.meta_res.annot.tbl.gz"  ] ||  [ ! -f "${sv_resdir}/${sv}.meta_res.eQTLs.txt.gz" ]
    then
        echo -e "$sv\tmain results missing"
    fi
    for i in `seq 1 10`
    do
       if [ ! -f "${sv_resdir}/${sv}.meta_res.eQTLs.perm${i}.txt.gz"  ]
        then
            echo -e "$sv\t${i}th permutation results missing"
        fi 
    done
    
done < all_bacs_${svtype}.txt

