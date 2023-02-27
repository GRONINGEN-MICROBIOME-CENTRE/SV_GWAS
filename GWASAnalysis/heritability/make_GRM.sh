#!/bin/bash
#SBATCH --job-name=GRM
#SBATCH --output=GRM.out
#SBATCH --error=GRM.err
#SBATCH --mem=35gb
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=2

#plink2 --bfile /groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/nonimputed/DAG3_genotyped_SNPs --keep /groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/noLLD/DAG3_filtered_with_relatives_passQC_noLLD.fam --make-bed --out DAG3_genotyped_SNPs_noLLD

geno_file=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/genotyped_SNPs/DAG3_genotyped_SNPs_noLLD
gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-sandreusanchez/Immuno_markers/Genetics/Heritability/Programs/gcta_1.93.2beta/gcta64

$gcta --bfile ${geno_file} --autosome --maf 0.05 --make-grm --out GRM_DAG3 --thread-num 2
$gcta --grm GRM_DAG3 --make-bK 0.05 --out GRM_DAG3_bK

cat "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/GCTA/GRM_DAG3" > mgrm.txt
cat "/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/DAG3/GCTA/GRM_DAG3_bK" >> mgrm.txt" 
