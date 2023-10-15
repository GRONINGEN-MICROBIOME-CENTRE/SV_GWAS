#!/usr/bin/env bash

#
# Format and filter the SV tables, prepare all necessary files and scripts, submit GWAS jobs
#

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/
script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/
genotype_dir=${d}/genotypes/

#
# 1. Prepare annotation and SV name conversion files
#

# dSVs
bash ${script_dir}/utils/create_sv_name_conversion_tables.sh \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_deletionStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvInfo/20230116_full_Informative_species_information.tsv \
dSV

# vSVs
bash ${script_dir}/utils/create_sv_name_conversion_tables.sh \
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
    python ${script_dir}/utils/rename_header_based_on_file.py stdin ${svtype}_name_conversion_table.txt 1 3 | \
    python3 ${script_dir}/utils/add_columns_from_file_v2.py -i stdin  -f pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 | \
    grep -w -v "NA" | \
    awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print }' \
    > ${cohort}.covariates.txt
    
    # 2.2 scale read counts
    Rscript ${script_dir}/utils/scale_read_count.R  ${cohort}.covariates.txt ${cohort}.covariates.scaled.txt
    mv ${cohort}.covariates.scaled.txt ${cohort}.covariates.txt
    
    # make the covariate file without header for fastGWA
    #tail -n+2 ${cohort}.covariates.txt > ${cohort}.covariates.noheader.txt

    cohorts=(LLD DAG3 300OB 500FG)
    
    # CLR transform the abundances
    mkdir tmp
    cut -f2-111 ${c}.covariates.txt > tmp/tmp.${c}.covariates.txt
    Rscript clr_transform.R tmp/tmp.${c}.covariates.txt tmp/tmp.${c}.covariates.clr.txt

    cut -f1 ${c}.covariates.txt | \
        paste - tmp/tmp.${c}.covariates.clr.txt <(cut -f112- ${c}.covariates.txt) | 
        sed "s:FID\t\t:FID\tIID\t:1" \
        > ${c}.covariates.txt
    tail -n+2 ${c}.covariates.txt > ${c}.covariates.noheader.txt

  
done

module load R/3.6.1-foss-2018a
cd ${d}

#
# 3. Filter the SV tables (call rate, frequency). INT the vSV table
#

# dSVs:
Rscript ${script_dir}/filter_SV_table.R \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_deletionStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/cleaned_file_list.tsv \
dSV

# vSVs:
Rscript ${script_dir}/filter_SV_table.R \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/20230116_full_variableStructuralVariation_12388samples.tsv \
/groups/umcg-fu/tmp01/projects/SV_GWAS/data/SvUnfiltered/SvProfile/cleaned_file_list.tsv \
vSV

# make the SV tables files without header for fastGWA
for f in  *.filtered.txt
do
 tail -n+2 $f > ${f%txt}filtered.noheader.txt
done

#
# 4. Count the number of overlapping samples
#

cd data_fastGWA
for cohort in ${cohorts[@]}
do
    echo "Counting samples for ${cohort}."
    echo "N samples with SVs --- of them N with genotypes --- of them N with phenotypes --- of them N with both phenotypes and genotypes"
    
    python3 ${script_dir}/utils/check_sample_overlap.py ${cohort}.dSV.filtered.txt 1  ${genotype_dir}/${cohort}/with_relatives/${cohort}_filtered_withrel.fam 1 pheno/${cohort}_pheno.txt 0 "id\twith_geno\twith_pheno" > sample_ids/${cohort}_id_overlap.txt
done


# make GRMs
gcta=${d}/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
for cohort in ${cohorts[@]}
do
    geno_file=${d}genotypes/${cohort}/with_relatives/${cohort}_filtered_withrel
    grm=${d}genotypes/${cohort}/with_relatives/GCTA/GRM_${cohort}
    mkdir ${d}genotypes/${cohort}/with_relatives/GCTA/
    
    $gcta --bfile ${geno_file} --extract ${d}/genotypes/genotyped_SNPs.snplist --maf 0.05 --make-grm --out $grm --thread-num 2

    $gcta --grm $grm --make-bK-sparse 0.05 --out ${grm}_sparse 
    
    $gcta --bfile ${geno_file} --extract ${d}/genotypes/genotyped_SNPs.snplist --maf 0.05 --make-grm-gz --out ${grm}.text --thread-num 2
done


#
# 5. Submit fastGWA jobs
#

svtype="dSV"
cd ${d}/scripts/scripts_fastGWA_${svtype}

cut -f1 ${d}/data_fastGWA/${svtype}_per_cohort.txt | tail -n+2 > all_bacs.${svtype}.txt
while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/run_${sv}.out \
    -e logs/run_${sv}.err \
    -J dsv_${sv} \
    -t 00:30:00 \
    ${script_dir}/run_fastGWA.sh $sv $svtype
done < all_bacs.${svtype}.txt

svtype="vSV"
cd ${d}/scripts/scripts_fastGWA_${svtype}

cut -f1 ${d}/data_fastGWA/${svtype}_per_cohort.txt | tail -n+2 > all_bacs.${svtype}.txt
while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/run_${sv}.out \
    -e logs/run_${sv}.err \
    -J dsv_${sv} \
    -t 00:20:00 \
    ${script_dir}/run_fastGWA.sh $sv $svtype
done < all_bacs.${svtype}.txt


#
# 6. Check the GWAS results files
#

result_dir="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/results_fastGWA/${svtype}/meta/"
while read line
do
    sv=$line
    sv_resdir=${result_dir}/${sv}/
    if [ ! -f "${sv_resdir}/${sv}.meta_res.annot.tbl.gz"  ] ||  [ ! -f "${sv_resdir}/${sv}.meta_res.annot.5e-8.tbl.gz" ] 
    then
        echo -e "$sv"
    elif [ `ls -l ${sv_resdir}/${sv}.meta_res.annot.tbl.gz | awk '{print $5}'` -lt 10000 ]
    then
       	echo -e "$sv"
    fi     
done < all_bacs.${svtype}.txt