#!/bin/bash
#SBATCH --job-name=h2__BL__
#SBATCH --output=logs/h2___BL__.out
#SBATCH --error=logs/h2___BL__.err
#SBATCH --mem=10gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1

gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/Programs/gcta_1.93.2beta/gcta64
d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
pheno_dir=${d}/data/plink19_format/
mgrm=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/GCTA/mgrm.txt
gender_file=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/DAG3_gender.txt

sv=$1
svtype=$2
echo "sv=$sv"
echo "svtype=$svtype"

res_dir=${d}/results/${svtype}/heritability/${sv}/
mkdir -p $res_dir

#get the column number for the SV in the SV table
col=`head -1 ${pheno_dir}/DAG3.${svtype}.filtered.plink19.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sv} | cut -d ":" -f1`

# get the species column numbers from the species abundance file 
sp=`grep -w "$sv" ${d}/data/${svtype}_name_conversion_table.txt | cut -f4 | uniq`
col_cov=`head -1 ${pheno_dir}/DAG3.covariates.plink19.txt | sed "s:\t:\n:g"  | grep -w -n ${sp} | cut -d ":" -f1`

# write the selected quantitative covariates to a tmp file (species abundance, age, read number)
awk -v c=$col_cov 'BEGIN {FS=OFS="\t"}; {print $1, $2, $c, $(NF-1), $(NF) }' ${pheno_dir}/DAG3.covariates.plink19.noheader.txt  \
>  ${res_dir}/tmp.qcovar.txt

# Run reml on family data
$gcta --reml \
    --mgrm ${mgrm} \
    --pheno ${pheno_dir}/DAG3.${svtype}.filtered.plink19.noheader.txt \
    --mpheno $col \
    --qcovar ${res_dir}/tmp.qcovar.txt \
    --covar ${gender_file} \
    --out ${res_dir}/${sv}_bKsK \
    --reml-no-constrain

rm ${res_dir}/tmp.qcovar.txt
