#!/usr/bin/env bash
#
# Run family-based heritability analysis of SVs using GCTA
#

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
#cur_script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cur_script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/heritability/


#
# Prepare SV and covariate files in plink v1.9 format
#
cd ${d}/data/
mkdir plink19_format

svtype=dSV
awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print  }' DAG3.${svtype}.filtered.txt > plink19_format/DAG3.${svtype}.filtered.plink19.txt
tail -n+2 plink19_format/DAG3.${svtype}.filtered.plink19.txt > plink19_format/DAG3.${svtype}.filtered.plink19.noheader.txt

svtype=vSV
awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print  }' DAG3.${svtype}.filtered.txt > plink19_format/DAG3.${svtype}.filtered.plink19.txt
tail -n+2 plink19_format/DAG3.${svtype}.filtered.plink19.txt > plink19_format/DAG3.${svtype}.filtered.plink19.noheader.txt

sed 's:"::g' abund/DAG3_abundances.tsv | \
sed '1s:.:#IID\t&:' | \
python ${cur_script_dir}/../gwas_scripts_misc/rename_header_based_on_file.py stdin dSV_name_conversion_table.txt 1 3 | \
python ${cur_script_dir}/../gwas_scripts_misc/add_columns_from_file.py -i stdin  -f pheno/DAG3_pheno.txt -f_m 0 -f_cols 1,2 | \
grep -w -v "NA" | \
awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print }' \
> plink19_format/DAG3.covariates.plink19.txt

tail -n+2 plink19_format/DAG3.covariates.plink19.txt > plink19_format/DAG3.covariates.plink19.noheader.txt

#
# vSVs
#
svtype=vSV
mkdir -p ${d}/scripts/scripts_${svtype}/heritability/logs/
cd ${d}/scripts/scripts_${svtype}/heritability/

grep "DAG3" ${d}/data/${svtype}_per_cohort.txt | cut -f1 | tail -n+2 > ${svtype}_per_cohort.DAG3.txt

#submit scripts
while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/h2_${sv}.out \
    -e logs/h2_${sv}.err \
    -J h2_${sv} \
    ${cur_script_dir}/run_h2_SV_TEMPLATE.sh "${sv}" ${svtype}
done < ${svtype}_per_cohort.DAG3.txt


# get results

while read line
do 
	sv=$line
	h=`grep "Sum of V(G)/Vp" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${svtype}_per_cohort.DAG3.txt


#
# dSVs
#

# submit
svtype=dSV
mkdir -p ${d}/scripts/scripts_${svtype}/heritability/logs/
cd ${d}/scripts/scripts_${svtype}/heritability/

grep "DAG3" ${d}/data/${svtype}_per_cohort.txt | cut -f1 | tail -n+2 > ${svtype}_per_cohort.DAG3.txt

#submit scripts
while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/h2_${sv}.out \
    -e logs/h2_${sv}.err \
    -J h2_${sv} \
    ${cur_script_dir}/run_h2_SV_TEMPLATE.sh "${sv}" ${svtype}
done < ${svtype}_per_cohort.DAG3.txt

# get results

while read line
do 
	sv=$line
	h=`grep "Sum of V(G)/Vp" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${svtype}_per_cohort.DAG3.txt
