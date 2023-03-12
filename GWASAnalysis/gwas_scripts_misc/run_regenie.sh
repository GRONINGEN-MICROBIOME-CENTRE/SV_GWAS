#!/bin/bash
#SBATCH --job-name=regenie
#SBATCH --output=logs/regenie_fp.out
#SBATCH --error=logs/regenie_fp.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ml regenie

d=/groups/umcg-fu/tmp01/projects/SV_GWAS/GWAS_tmp/
sp="F.prausnitzii"
cohort="DAG3"

#col_cov=`head -1 ${d}/data/${cohort}.covariates.plink19.txt | sed "s:\t:\n:g"  | grep -w -n ${sp} | cut -d ":" -f1`

# write the selected quantitative covariates to a tmp file (species abundance, age, read number)
#awk -v c=$col_cov 'BEGIN {FS=OFS="\t"}; {print $1, $2, $c, $(NF-1), $(NF) }' ${d}/data/${cohort}.covariates.plink19.txt \
#>  ${d}/data/per_species/${cohort}/${cohort}.${sp}.covariates.txt

# don't forget to add sex to covars
# use all genotypes with relatives

regenie \
  --step 1 \
  --bed ${d}/genotypes/${cohort}/${cohort}_filtered \
  --extract ${d}/genotypes/${cohort}/genotyped_SNPs/${cohort}_genotyped_SNPs.snplist \
  --covarFile ${d}/data/per_species/${cohort}/${cohort}.${sp}.covariates.txt \
  --phenoFile ${d}/data/per_species/${cohort}/${cohort}.dSV.${sp}.txt \
  --bsize 100 \
  --bt --lowmem --threads 4 \
  --lowmem-prefix tmp_rg \
  --out fit_bin_out 
  
regenie \
  --step 2 \
  --bed ${d}/genotypes/${cohort}/${cohort}_filtered \
  --covarFile ${d}/data/per_species/${cohort}/${cohort}.${sp}.covariates.txt \
  --phenoFile ${d}/data/per_species/${cohort}/${cohort}.dSV.${sp}.txt \
  --bsize 200 --threads 4 \
  --bt \
  --firth --approx \
  --pred fit_bin_out_pred.list \
  --out test_bin_out_firth