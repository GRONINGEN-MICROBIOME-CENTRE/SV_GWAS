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

while read line
do
    sp=$line
    # merge sorted files
    meta_comb_dir=/data/umcg-tifn/SV/SV_GWAS/results_all_summary_stats/dSV/meta_combined/${sp}/
    #gunzip  /data/umcg-tifn/SV/SV_GWAS/results_all_summary_stats/dSV/meta/${sp}\:*/*meta_res.annot.tbl.gz

    mkdir ${meta_comb_dir}
    echo -e "SV_id\tSNP\tchromosome\tposition\teffect_allele\tother_allele\tmeta_samplesize\tmeta_zscore\tmeta_pvalue\tmeta_EAF\tcohorts_with_SV\teffect_direction_per_cohort\tsamplesize_per_cohort" \    
    > ${TMPDIR}/${sp}.dSV.meta-analysis.txt

    # CHECK!
    sort -m -k9,9g -T $TMPDIR -S 100G --buffer-size=1000  /data/umcg-tifn/SV/SV_GWAS/results_all_summary_stats/dSV/meta/${sp}\:*/*meta_res.annot.fmt.txt \
    >> ${TMPDIR}/${sp}.dSV.meta-analysis.txt

    gzip -c ${TMPDIR}/${sp}.dSV.meta-analysis.txt > ${meta_comb_dir}/${sp}.dSV.meta-analysis.txt.gz
    gzip  /data/umcg-tifn/SV/SV_GWAS/results_all_summary_stats/dSV/meta/${sp}\:*/*meta_res.annot.fmt.txt
    
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