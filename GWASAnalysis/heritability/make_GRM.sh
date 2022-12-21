#!/bin/bash
#SBATCH --job-name=GRM
#SBATCH --output=GRM.out
#SBATCH --error=GRM.err
#SBATCH --mem=35gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=2

geno_file=/groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/nonimputed/DAG3_genotyped_SNPs
gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/Programs/gcta_1.93.2beta/gcta64

$gcta --bfile ${geno_file} --autosome --maf 0.05 --make-grm --out GRM_DAG3 --thread-num 2
$gcta --grm GRM_DAG3 --make-bK 0.05 --out GRM_DAG3_bK