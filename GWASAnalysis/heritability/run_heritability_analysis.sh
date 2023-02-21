#!/usr/bin/env bash
#
# Run family-based heritability analysis of SVs using GCTA
#

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
#cur_script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cur_script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/heritability/
#
# vSVs
#
svtype=vSV
mkdir -p ${d}/scripts/scripts_${svtype}/heritability/logs/
cd ${d}/scripts/scripts_${svtype}/heritability/

grep "DAG3" ${d}/data/${svtype}_per_cohort_all_cohorts.txt | cut -f1 | tail -n+2 > ${d}/data/${svtype}_per_cohort.DAG3.txt
#submit scripts
rm ${svtype}.split.*
split -l$((`wc -l <  ${d}/data/${svtype}_per_cohort.DAG3.txt`/10))  ${d}/data/${svtype}_per_cohort.DAG3.txt ${svtype}.split. -da 2
for f in vSV.split.*
do
	sed "s:__BL__:${f}:g" ${cur_script_dir}/run_h2_SV_TEMPLATE.sh ${svtype} > run_h2_${svtype}_${f}.sh
	sbatch run_h2_${svtype}_${f}.sh
done

# get results

while read line
do 
	sv=$line
	h=`grep "Sum of V(G)/Vp" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${d}/data/${svtype}_per_cohort.DAG3.txt


#
# dSVs
#

# submit
svtype=dSV
mkdir -p ${d}/scripts/scripts_${svtype}/heritability/logs/
cd ${d}/scripts/scripts_${svtype}/heritability/

grep "DAG3" ${d}/data/${svtype}_per_cohort_all_cohorts.txt | cut -f1 | tail -n+2 > ${d}/data/${svtype}_per_cohort.DAG3.txt

#submit scripts
rm ${svtype}.split.*
split -l$((`wc -l < ${d}/data/${svtype}_per_cohort.DAG3.txt`/10)) ${d}/data/${svtype}_per_cohort.DAG3.txt ${svtype}.split. -da 2
for f in vSV.split.*
do
	sed "s:__BL__:${f}:g" ${cur_script_dir}/run_h2_SV_TEMPLATE.sh ${svtype} > run_h2_${svtype}_${f}.sh
	sbatch run_h2_${svtype}_${f}.sh
done

# get results

while read line
do 
	sv=$line
	h=`grep "Sum of V(G)/Vp" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results/${svtype}/heritability/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${d}/data/${svtype}_per_cohort.DAG3.txt
