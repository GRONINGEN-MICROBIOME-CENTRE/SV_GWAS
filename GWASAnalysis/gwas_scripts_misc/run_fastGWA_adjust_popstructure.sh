d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/"
script_dir="${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/"
gcta=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v3/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
pheno_dir=${d}/data_fastGWA/
cohort="DAG3"

geno_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_filtered_withrel
grm=${d}/genotypes/${cohort}/with_relatives/GCTA/GRM_${cohort}_sparse
gender_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_gender.txt 
covar_file=${pheno_dir}/${cohort}.covariates.Fprau_population.noheader.txt

svs=("F.prausnitzii:102" "F.prausnitzii:101" "F.prausnitzii:9" "F.prausnitzii:69")
svtype=dSV

gcta_mode="--fastGWA-mlm-binary"



for sv in ${svs[@]}
do
    echo -e "\n\nRunning the analysis for ${sv}\n\n"

    res_dir=${d}/results_fastGWA/${svtype}/adj_population_PCs/${sv}/
    mkdir -p $res_dir

    #get the column number for the SV in the SV table
    col=`head -1 ${pheno_dir}/${cohort}.${svtype}.filtered.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sv} | cut -d ":" -f1`


    $gcta \
      --bfile $geno_file \
      --grm-sparse ${grm} \
      ${gcta_mode} \
      --pheno ${pheno_dir}/${cohort}.${svtype}.filtered.noheader.txt \
      --mpheno $col \
      --qcovar ${covar_file} \
      --covar ${gender_file} \
      --out ${res_dir}/${sv}_adj_population_PCs.chr9 \
      --chr 9
    
   
    sort -k13,13g ${res_dir}/${sv}_adj_population_PCs.chr9.fastGWA | gzip -c > ${res_dir}/${sv}_adj_population_PCs.chr9.fastGWA.gz
    rm ${res_dir}/${sv}_adj_population_PCs.chr9.fastGWA
    
done


sv="F.prausnitzii:33"
svtype=vSV

gcta_mode="--fastGWA-mlm"


res_dir=${d}/results_fastGWA/${svtype}/adj_population_PCs/${sv}/
mkdir -p $res_dir

#get the column number for the SV in the SV table
col=`head -1 ${pheno_dir}/${cohort}.${svtype}.filtered.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sv} | cut -d ":" -f1`


$gcta \
  --bfile $geno_file \
  --grm-sparse ${grm} \
  ${gcta_mode} \
  --pheno ${pheno_dir}/${cohort}.${svtype}.filtered.noheader.txt \
  --mpheno $col \
  --qcovar ${covar_file} \
  --covar ${gender_file} \
  --out ${res_dir}/${sv}_adj_population_PCs.chr9 \
  --chr 9


sort -k13,13g ${res_dir}/${sv}_adj_population_PCs.chr9.fastGWA | gzip -c > ${res_dir}/${sv}_adj_population_PCs.chr9.fastGWA.gz
rm ${res_dir}/${sv}_adj_population_PCs.chr9.fastGWA



