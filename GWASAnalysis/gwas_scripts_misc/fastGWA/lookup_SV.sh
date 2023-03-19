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

module load Metal


d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir="${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/"
gcta=${d}/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
pheno_dir=${d}/data_fastGWA/

sv=$1
snp=$2
svtype=$3
cohorts=$4
IFS=',' read -r -a cohorts_with_sv <<< $cohorts
echo "sv=$sv, snp=$snp, svtype=$svtype"

echo $snp > ${d}/scripts/scripts_fastGWA_${svtype}/snps/$snp

# get species name
sp=`grep -w "$sv" ${d}/data_fastGWA/${svtype}_name_conversion_table.txt | cut -f4 | uniq` 
echo "Species name: $sp"
all_nsamples=()
all_cohorts=()
all_zscores=()
all_pvals=()
# create the output folder for meta-analysis results
meta_out_dir=${d}/results_fastGWA/${svtype}/meta/${sv}/
meta_out_filebase=${meta_out_dir}/${sv}.meta_res
mkdir -p ${d}/results_fastGWA/${svtype}/meta/${sv}/

# prepare metal script header
metal_script=${d}/scripts/scripts_fastGWA_${svtype}/metal_per_sv/${sv}.metal.txt
cat ${d}/scripts/scripts_fastGWA_${svtype}/metal_header.txt > $metal_script

# check in which cohorts the SV is present
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

    res_dir=${d}/results_fastGWA/${svtype}/${cohort}/${sv}/
    mkdir -p $res_dir
    geno_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_filtered_withrel
    grm=${d}/genotypes/${cohort}/with_relatives/GCTA/GRM_${cohort}_sparse
    gender_file=${d}/genotypes/${cohort}/with_relatives/${cohort}_gender.txt 
    
    #get the column number for the SV in the SV table
    col=`head -1 ${pheno_dir}/${cohort}.${svtype}.filtered.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sv} | cut -d ":" -f1`

    # get the species column numbers from the species abundance file 
    col_cov=`head -1 ${pheno_dir}/${cohort}.covariates.txt | sed "s:\t:\n:g"  | grep -w -n ${sp} | cut -d ":" -f1`

    # write the selected quantitative covariates to a tmp file (species abundance, age, read number)
    awk -v c=$col_cov 'BEGIN {FS=OFS="\t"}; {print $1, $2, $c, $(NF-1), $(NF) }' ${pheno_dir}/${cohort}.covariates.noheader.txt  \
    >  ${res_dir}/tmp.qcovar.txt

    #
    # run real GWAS analysis
    #
    $gcta \
      --bfile $geno_file \
      --grm-sparse ${grm} \
      ${gcta_mode} \
      --pheno ${pheno_dir}/${cohort}.${svtype}.filtered.noheader.txt \
      --mpheno $col \
      --qcovar ${res_dir}/tmp.qcovar.txt \
      --covar ${gender_file} \
      --extract ${d}/scripts/scripts_fastGWA_${svtype}/snps/$snp \
      --out ${res_dir}/${sv}.${snp}
    
    rc=$?
    echo "$sv fastGWA return code: $rc"    
    
    if [ $rc -eq 0 ]
    then
        # Number of samples with association results for this SV for this cohort
        n=`head -2 ${res_dir}/${sv}.${snp}.fastGWA | tail -1 | awk '{print $6}'`
        z=`head -2 ${res_dir}/${sv}.${snp}.fastGWA | tail -1 | awk '{print $8}'`
        pval=`head -2 ${res_dir}/${sv}.${snp}.fastGWA | tail -1 | awk '{print $10}'`
        all_cohorts+=( $cohort )
        all_nsamples+=( $n )
        all_zscores+=( $z )
        all_pvals+=( $pval )
    fi    
    
    gzip -f ${res_dir}/${sv}.${snp}.fastGWA
    
    # append the per cohort result location to the metal script
    echo -e "PROCESS\t${res_dir}/${sv}.${snp}.fastGWA.gz" >> $metal_script


done



# Real GWAS meta-analysis
echo -e "OUTFILE\t${meta_out_filebase} .tbl\nANALYZE HETEROGENEITY\nQUIT" >> $metal_script
metal $metal_script
echo "${sv}, real analysis metal return code: $?"

for cohort in ${cohorts_with_sv[@]}
do 
    res_dir=${d}/results_fastGWA/${svtype}/${cohort}/${sv}/
    rm ${res_dir}/${sv}.${snp}.fastGWA.gz
done

het_pval=`cut -f11 ${meta_out_filebase}1.tbl | tail -n+2`
cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
zscores_joined=`printf -v var '%s,' "${all_zscores[@]}"; echo "${var%,}"`
pvals_joined=`printf -v var '%s,' "${all_pvals[@]}"; echo "${var%,}"`

echo -e "$zscores_joined $pvals_joined $het_pval"