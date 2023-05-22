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

svtype=$1
echo "SV type=${svtype}"
d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/results_fastGWA/
while read line
do
    sp=$line
    # merge sorted files
    meta_comb_dir=${d}/${svtype}/meta_combined/${sp}/
   
    mkdir -p ${meta_comb_dir}
    
    echo -e "SV_id\tSNP_position_hg19\teffect_allele\tother_allele\tmeta_beta\tmeta_SE\tmeta_pvalue\teffect_direction_per_cohort\theterogeneity_pvalue\tmeta_samplesize\tcohorts_with_SV\tsamplesize_per_cohort" \
    > ${TMPDIR}/${sp}.${svtype}.fastGWA.meta-analysis.txt
    
    gunzip ${d}/${svtype}/meta/${sp}\:*/*meta_res.annot.tbl.gz
    sort -m -k7,7g -T $TMPDIR -S 30G --buffer-size=1000  ${d}/${svtype}/meta/${sp}\:*/*meta_res.annot.tbl |
    cut -f1-8,12- \
    >> ${TMPDIR}/${sp}.${svtype}.fastGWA.meta-analysis.txt

    gzip -c ${TMPDIR}/${sp}.${svtype}.fastGWA.meta-analysis.txt > ${meta_comb_dir}/${sp}.${svtype}.fastGWA.meta-analysis.txt.gz
    gzip  ${d}/${svtype}/meta/${sp}\:*/*meta_res.annot.tbl
    
    python3 reformat_summary_stats.py $sp dSV

    # TODO: Check what files to remove!

    # retCode=$?
    # if [ $retCode -ne 0 ]; then
    # echo "Error!"
    # else
    # echo "Finished succesfully"
    # rm ${f}/*meta-analysis.txt.gz
    # fi

    tail -n+2 ${f}/${sp}.dSV.meta-analysis.annot.txt | awk 'BEGIN {FS=OFS="\t"}; {if ($9 < 0.01) print $1, $2,$3,$4,$9}' >> ${for_plot}
    gzip ${f}/${sp}.dSV.meta-analysis.annot.txt
    
done < dSV_species.txt