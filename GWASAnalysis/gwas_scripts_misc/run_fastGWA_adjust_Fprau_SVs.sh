#!/bin/bash
#SBATCH --job-name=SV
#SBATCH --output=logs/run_GWAS.out
#SBATCH --error=logs/run_GWAS_.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml Metal/2020-05-05-foss-2018b


d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir="${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/"
gcta=${d}/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
pheno_dir=${d}/data_fastGWA/

sv=$1
svtype=$2
echo "sv=$sv, svtype=$svtype"

# get species name
all_nsamples=()
all_cohorts=()
# create the output folder for meta-analysis results
meta_out_dir=${d}/results_fastGWA/${svtype}/meta/${sv}/
meta_out_filebase=${meta_out_dir}/${sv}.meta_res.adj_Fprau_SVs
mkdir -p ${d}/results_fastGWA/${svtype}/meta/${sv}/

# prepare metal script header
metal_script=${d}/scripts/scripts_fastGWA_${svtype}/metal_per_sv/${sv}.metal.txt
cat ${d}/scripts/scripts_fastGWA_${svtype}/metal_header.txt > $metal_script

# check in which cohorts the SV is present
IFS=',' read -ra cohorts_with_sv <<< `grep -w $sv ${d}/data_fastGWA/${svtype}_per_cohort.txt | cut -f6`
cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
echo "Cohorts with SV: $cohorts_joined"

if [ $svtype == "dSV" ]
then
    gcta_mode="--fastGWA-mlm-binary"
else
    gcta_mode="--fastGWA-mlm"
fi
#
# Primary (real) GWAS:
#

for cohort in ${cohorts_with_sv[@]}
do
    echo -e "\n\nRunning the analysis for ${cohort}\n\n"

    res_dir=${d}/results_fastGWA/${svtype}/${cohort}/${sv}/adj_Fprau_SVs/

    geno_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_filtered_withrel
    grm=${d}/genotypes/${cohort}/with_relatives/GCTA/GRM_${cohort}_sparse
    gender_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_gender.txt 
    
    #
    # run real GWAS analysis
    #
    $gcta \
      --bfile $geno_file \
      --grm-sparse ${grm} \
      ${gcta_mode} \
      --pheno ${res_dir}/tmp.pheno.txt \
      --qcovar ${res_dir}/tmp.qcovar.txt \
      --covar ${gender_file} \
      --out ${res_dir}/${sv} \
      --chr 9
    
    rc=$?
    echo "$sv fastGWA return code: $rc"    
    
    if [ $rc -eq 0 ]
    then
        # Number of samples with association results for this SV for this cohort
        n=`head -2 ${res_dir}/${sv}.fastGWA | tail -1 | awk '{print $6}'`
        all_cohorts+=( $cohort )
        all_nsamples+=( $n )
    fi    
    
    
    gzip -f ${res_dir}/${sv}.fastGWA
    
    # append the per cohort result location to the metal script
    echo -e "PROCESS\t${res_dir}/${sv}.fastGWA.gz" >> $metal_script


done



# Real GWAS meta-analysis
echo -e "OUTFILE\t${meta_out_filebase} .tbl\nANALYZE HETEROGENEITY\nQUIT" >> $metal_script
metal $metal_script
echo "${sv}, real analysis metal return code: $?"

