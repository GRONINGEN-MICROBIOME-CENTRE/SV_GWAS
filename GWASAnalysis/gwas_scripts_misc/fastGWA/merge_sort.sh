#!/bin/bash
#SBATCH --job-name=fdr
#SBATCH --output=logs/fdr.out
#SBATCH --error=logs/fdr.err
#SBATCH --time=60:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=85gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp=90gb


svtype=$2

meta_comb_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/results_fastGWA/${svtype}/meta_combined/

mkdir ${meta_comb_dir}
cd ${meta_comb_dir}

gunzip -f ../meta/*/*.meta_res.annot.5e-8.tbl.gz

echo -e "SV\tSNP\tEffect_allele\tOther_allele\tN\tZ\tPvalue\tDirection\tCohorts\tN_per_cohort" \
> ${svtype}.fastGWA.5e-08.txt

sort -m -k7,7g -T $TMPDIR -S 70G --buffer-size=1000 ../meta/*/*.meta_res.annot.5e-8.tbl >> ${svtype}.fastGWA.5e-08.txt

gzip ../meta/*/*.meta_res.annot.5e-8.tbl