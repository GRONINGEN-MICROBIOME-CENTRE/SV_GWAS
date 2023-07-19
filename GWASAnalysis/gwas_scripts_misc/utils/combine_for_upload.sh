#!/bin/bash
#SBATCH --job-name=combine
#SBATCH --output=combine.out
#SBATCH --error=combine.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=80gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

svtype=dSV
d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/results_fastGWA/${svtype}/meta_combined/

zcat ${d}/summary_stats/A.hadrus.${svtype}.fastGWA.meta-analysis.annot.txt.gz | head -1 > ${d}/all_species.${svtype}.fastGWA.meta-analysis.p0.05.txt

for f in ${d}/summary_stats/*.${svtype}.fastGWA.meta-analysis.annot.txt.gz
do
    echo $f
    zcat $f | tail -n+2 | awk 'BEGIN {FS=OFS="\t"}; {if ($9 < 0.05) print $0}' >> ${d}/all_species.${svtype}.fastGWA.meta-analysis.p0.05.txt

done

bzip2 -9 ${d}/all_species.${svtype}.fastGWA.meta-analysis.p0.05.txt

svtype=vSV
d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/results_fastGWA/${svtype}/meta_combined/

zcat ${d}/summary_stats/A.hadrus.${svtype}.fastGWA.meta-analysis.annot.txt.gz | head -1 > ${d}/all_species.${svtype}.fastGWA.meta-analysis.p0.05.txt

for f in ${d}/summary_stats/*.${svtype}.fastGWA.meta-analysis.annot.txt.gz
do
    echo $f
    zcat $f | tail -n+2 | awk 'BEGIN {FS=OFS="\t"}; {if ($9 < 0.05) print $0}' >> ${d}/all_species.${svtype}.fastGWA.meta-analysis.p0.05.txt

done

bzip2 -9 ${d}/all_species.${svtype}.fastGWA.meta-analysis.p0.05.txt