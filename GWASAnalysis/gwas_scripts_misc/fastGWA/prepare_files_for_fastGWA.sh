d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/

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
    
    # 2.1 Format covariates in plink 1.9, make sure they include all samples with relatives
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

cd data_fastGWA
for cohort in ${cohorts[@]}
do
    echo "Counting samples for ${cohort}."
    echo "N samples with SVs --- of them N with genotypes --- of them N with phenotypes --- of them N with both phenotypes and genotypes"
    
    python3 ${script_dir}/check_sample_overlap.py ${cohort}.dSV.filtered.txt 1  ${genotype_dir}/${cohort}/with_relatives/${cohort}_filtered_withrel.fam 1 pheno/${cohort}_pheno.txt 0 "id\twith_geno\twith_pheno" > sample_ids/${cohort}_id_overlap.txt


done

for f in  *.filtered.txt
do
 tail -n+2 $f > ${f%txt}filtered.noheader.txt
done
# make GRMs
gcta=${d}/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
for cohort in ${cohorts[@]}
do
    geno_file=${d}genotypes/${cohort}/with_relatives/${cohort}_filtered_withrel
    grm=${d}genotypes/${cohort}/with_relatives/GCTA/GRM_${cohort}
    mkdir ${d}genotypes/${cohort}/with_relatives/GCTA/

    gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/Programs/gcta_1.93.2beta/gcta64

    $gcta --bfile ${geno_file} --extract ${d}/genotypes/genotyped_SNPs.snplist --maf 0.05 --make-grm --out $grm --thread-num 2

    $gcta --grm $grm --make-bK-sparse 0.05 --out ${grm}_sparse 

done