ml GCCcore/11.3.0

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
pheno_dir=${d}/data/plink19_format/
score=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SCORE/build/SCORE

sv=$1
svtype=$2
echo "sv=$sv"
echo "svtype=$svtype"

#svtype="dSV"
#sv="F.prausnitzii:102"
res_dir=${d}/results/${svtype}/heritability_score/${sv}/
mkdir -p $res_dir

#get the column number for the SV in the SV table
col=`head -1 ${pheno_dir}/DAG3.${svtype}.filtered.plink19.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sv} | cut -d ":" -f1`
tmp_pheno=${res_dir}/tmp.pheno.txt
awk -v c=$col 'BEGIN {FS=OFS="\t"}; {print $1, $2, $c}' ${pheno_dir}/DAG3.${svtype}.filtered.plink19.txt | grep -v -w "NA" | sed 's:#FID:FID:g' > ${tmp_pheno}

# get the species column numbers from the species abundance file 
sp=`grep -w "$sv" ${d}/data/${svtype}_name_conversion_table.txt | cut -f4 | uniq`
covar_file=/groups/umcg-fu/tmp01/projects/SV_GWAS/GWAS_tmp/data/per_species/DAG3/DAG3.${sp}.covariates.txt

${score} \
  -g /groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/LL_genotypes/UGLI_mgs_v1/nonimputed/DAG3_genotyped_SNPs.filtered.norel \
  -p ${tmp_pheno} \
  -c ${covar_file} \
  -o ${res_dir}/${sv}.h2_score.txt

rm $tmp_pheno