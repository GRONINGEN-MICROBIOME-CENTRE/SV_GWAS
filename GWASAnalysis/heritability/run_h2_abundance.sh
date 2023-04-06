#!/bin/bash
#SBATCH --job-name=h2__BL__
#SBATCH --output=logs/h2_abund___BL__.out
#SBATCH --error=logs/h2_abund___BL__.err
#SBATCH --mem=10gb
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
gcta=${d}/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1

pheno_dir=${d}/data_fastGWA/
mgrm=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/GCTA/mgrm.txt
grm=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/GCTA/GRM_DAG3_norel
gender_file=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/DAG3_gender.txt

sp=$1
res_dir=${d}/results_fastGWA/heritability_GCTA_abundance/${sp}/
mkdir -p $res_dir

#get the column number for the species in the table
col=`head -1 ${pheno_dir}/abund/sv.abun.dag3.renamed.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sp} | cut -d ":" -f1`

# write the selected quantitative covariates to a tmp file (species abundance, age, read number)
awk  'BEGIN {FS=OFS="\t"}; {print $1, $2, $(NF-1), $(NF) }' ${pheno_dir}/DAG3.covariates.noheader.txt  \
>  ${res_dir}/tmp.qcovar.txt

# Run reml on family data
$gcta --reml \
    --mgrm ${mgrm} \
    --pheno ${pheno_dir}/abund/sv.abun.dag3.renamed.txt \
    --mpheno $col \
    --qcovar ${res_dir}/tmp.qcovar.txt \
    --covar ${gender_file} \
    --out ${res_dir}/${sp}_bKsK 


$gcta --reml \
    --grm ${grm} \
    --pheno ${pheno_dir}/abund/sv.abun.dag3.renamed.txt \
    --mpheno $col \
    --qcovar ${res_dir}/tmp.qcovar.txt \
    --covar ${gender_file} \
    --out ${res_dir}/${sp}_norel 
    
rm ${res_dir}/tmp.qcovar.txt
