d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir="${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/"
gcta=${d}/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
pheno_dir=${d}/data_fastGWA/

cohort="DAG3"
res_dir=${d}/results_fastGWA/shortbred/
mkdir -p $res_dir
geno_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_filtered_withrel
grm=${d}/genotypes/${cohort}/with_relatives/GCTA/GRM_${cohort}_sparse
gender_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_gender.txt 


# write the selected quantitative covariates to a tmp file (species abundance, age, read number)
awk 'BEGIN {FS=OFS="\t"}; {print $1, $2, $(NF-1), $(NF) }' ${pheno_dir}/${cohort}.covariates.noheader.txt  \
>  ${res_dir}/tmp.qcovar.txt

genes=(HTF-238_02532 HTF-238_02533 HTF-238_02539 HTF-238_02547 HTF-238_02549 HTF-238_02550 HTF-238_02552 HTF-238_02553 HTF-238_02554 HTF-238_02556 Int-Tn_4 afr_2-iolU agaC agaS immR_6 gatY-kbaY gatZ-kbaZ lacC agaV malX agaF agaD nagA ndoA_2 pgmB ptsH hpaA-rhaR_2)
for gene in "${genes[@]}"
do
#gene=agaC
col=`head -1 ${pheno_dir}/shortbred_final.res.log.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${gene} | cut -d ":" -f1`


$gcta \
      --bfile $geno_file \
      --grm-sparse ${grm} \
      --fastGWA-mlm \
      --pheno ${pheno_dir}/shortbred_final.res.log.noheader.txt \
      --mpheno $col \
      --qcovar  ${res_dir}/tmp.qcovar.txt \
      --covar ${gender_file} \
      --out ${res_dir}/${gene}
      
sort -k10,10g ${res_dir}/${gene}.fastGWA | gzip -c > ${res_dir}/${gene}.fastGWA.gz
rm ${res_dir}/${gene}.fastGWA

#done


zcat ${res_dir}/${gene}.fastGWA.gz | \
awk 'BEGIN {FS=OFS="\t"}; {if ($10 < 5e-08) print}' |  \
python3 ${script_dir}/utils/add_columns_from_file_v2.py \
-i stdin -i_m 1 \
-f /groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/300TZFG_replication/results_fastGWA/shortbred/${gene}.fastGWA.gz -f_m 1 \
-f_cols 4,7,9 \
> ${res_dir}/${gene}.fastGWA.replication.txt

done