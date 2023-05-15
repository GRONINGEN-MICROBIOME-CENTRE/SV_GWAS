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
mkdir -p ${d}/scripts/heritability/${svtype}/logs/
cd ${d}/scripts/heritability/${svtype}/

grep "DAG3" ${d}/data_fastGWA/${svtype}_per_cohort.txt | cut -f1 | tail -n+2 > ${svtype}_per_cohort.DAG3.txt

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
	h=`grep "Sum of V(G)/Vp" ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${svtype}_per_cohort.DAG3.txt


#
# dSVs
#

# submit
svtype=dSV
mkdir -p ${d}/scripts/heritability/${svtype}/logs/
cd ${d}/scripts/heritability/${svtype}/

grep "DAG3" ${d}/data_fastGWA/${svtype}_per_cohort.txt | cut -f1 | tail -n+2 > ${svtype}_per_cohort.DAG3.txt

#submit scripts
while read line
do 
  sv="${line}"
  echo $sv
  sbatch \
    -o logs/h2_${sv}.out \
    -e logs/h2_${sv}.err \
    -J h2_${sv} \
    ${cur_script_dir}/run_h2_SV.sh "${sv}" ${svtype}
done < ${svtype}_per_cohort.DAG3.txt

# get results_fastGWA

while read line
do 
	sv=$line
	h=`grep "Sum of V(G)/Vp" ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_bKsK.hsq | awk '{print $4,$5}'`
	p=`grep "Pval" ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_bKsK.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${svtype}_per_cohort.DAG3.txt > ${svtype}_h2_bKsK_res.txt

while read line
do 
	sv=$line
	h=`grep "V(G)/Vp" ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_norel.hsq | awk '{print $2,$3}'`
	p=`grep "Pval" ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_norel.hsq | awk '{print $2}'`
	n=`tail -n 1 ${d}/results_fastGWA/${svtype}/heritability_GCTA/${sv}/${sv}_norel.hsq | awk '{print $2}'`
	echo "$sv $n $h $p"
done < ${svtype}_per_cohort.DAG3.txt > ${svtype}_h2_norel_res.txt

