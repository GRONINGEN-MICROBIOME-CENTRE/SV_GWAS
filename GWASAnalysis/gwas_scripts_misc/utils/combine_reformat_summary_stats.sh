#!/bin/bash
#SBATCH --job-name=combine_SV
#SBATCH --output=logs/combine_SV.out
#SBATCH --error=logs/combine_SV.err
#SBATCH --time=150:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=110gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp=50gb

#
# combine full summary stats of individual SVs into files per species  
#

sp=$1
svtype=$2
echo "SV type=${svtype}, species=$sp"
d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/results_fastGWA/
script_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/utils/
# merge sorted files
meta_comb_dir=${d}/${svtype}/meta_combined/${sp}/

mkdir -p ${meta_comb_dir}

echo -e "SV_id\tSNP_position_hg19\teffect_allele\tother_allele\tmeta_beta\tmeta_SE\tmeta_pvalue\teffect_direction_per_cohort\theterogeneity_pvalue\tmeta_samplesize\tcohorts_with_SV\tsamplesize_per_cohort" \
> ${TMPDIR}/${sp}.${svtype}.fastGWA.meta-analysis.txt

gunzip ${d}/${svtype}/meta/${sp}\:*/*meta_res.annot.tbl.gz
sort -m -k7,7g -T $TMPDIR -S 100G --buffer-size=1000  ${d}/${svtype}/meta/${sp}\:*/*meta_res.annot.tbl |
cut -f1-8,12- \
>> ${TMPDIR}/${sp}.${svtype}.fastGWA.meta-analysis.txt

gzip -c ${TMPDIR}/${sp}.${svtype}.fastGWA.meta-analysis.txt > ${meta_comb_dir}/${sp}.${svtype}.fastGWA.meta-analysis.txt.gz
gzip  ${d}/${svtype}/meta/${sp}\:*/*meta_res.annot.tbl

python3 ${script_dir}/reformat_summary_stats.py $sp ${svtype} /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/

rm ${meta_comb_dir}/${sp}.${svtype}.fastGWA.meta-analysis.txt.gz
gzip ${meta_comb_dir}/${sp}.${svtype}.fastGWA.meta-analysis.annot.txt
    
