script_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SV_GWAS/GWASAnalysis/
d=/groups/umcg-fu/tmp01/projects/SV_GWAS/GWAS_tmp/
svtype="dSV"


cohorts=("300OB" "500FG" "LLD" "DAG3")
for cohort in ${cohorts[@]}
do
  
  python3 ${script_dir}/gwas_scripts_misc/prepare_phenotypes_for_regenie.py \
    ${d}/data/${cohort}.${svtype}.filtered.txt \
    ${d}/data/${svtype}_per_cohort.txt \
    ${cohort} \
    ${svtype} \
    ${d}/data/per_species/

  awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "FID\tIID"; else $1 = "0\t" $1; print }' ${d}/data/${cohort}.covariates.txt > ${d}/data/plink19_format/${cohort}.covariates.plink19.txt

done


for cohort in ${cohorts[@]}
do
  species_list=${d}/data/per_species/${cohort}/${cohort}.${svtype}.species_list
  rm $species_list
  for f in ${d}/data/per_species/${cohort}/${cohort}.${svtype}.*txt
    do tmp=${f##*SV.}
    sp=${tmp%".txt"}
    echo $sp >> $species_list
  done
  
  while read line
  do
    sp=$line
    awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "FID\tIID"; else $1 = "0\t" $1; print }' \
      ${d}/data/per_species/${cohort}/${cohort}.${svtype}.${sp}.txt \
      > ${d}/data/per_species/${cohort}/${cohort}.${svtype}.${sp}.plink19.txt
    
    col_cov=`head -1 ${d}/data/plink19_format/${cohort}.covariates.plink19.txt | sed "s:\t:\n:g"  | grep -w -n ${sp} | cut -d ":" -f1`
    awk -v c=$col_cov 'BEGIN {FS="\t"; OFS=" "}; {print $1, $2, $c, $(NF-1), $(NF) }' \
      ${d}/data/plink19_format/${cohort}.covariates.plink19.txt | \
    python2 ${script_dir}/gwas_scripts_misc/add_columns_from_file.py \
    -i stdin -i_m 1 \
    -d " " \
    -f ${d}/genotypes/${cohort}/${cohort}_filtered.fam -f_m 1 \
    -f_cols 4 | \
    sed '1s:NA$:sex:;s: :\t:g' \
    > ${d}/data/per_species/${cohort}/${cohort}.${sp}.covariates.txt
    
  done < $species_list
done 





#
sp="F.prausnitzii"
cohort="300OB"

  sbatch \
    -o logs/run_${sp}_${cohort}.out \
    -e logs/run_${sp}_${cohort}.err \
    -J ${sp}_${cohort} \
    run_regenie.sh "${sp}" ${cohort}

