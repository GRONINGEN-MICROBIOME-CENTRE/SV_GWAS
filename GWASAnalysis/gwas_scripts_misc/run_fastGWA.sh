#!/bin/bash
#SBATCH --job-name=SV
#SBATCH --output=logs/run_GWAS.out
#SBATCH --error=logs/run_GWAS_.err
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Metal


d="/groups/umcg-lifelines/tmp01/projects/${cohort}_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir="${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/"
gcta=/groups/umcg-lifelines/tmp01/projects/${cohort}_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/tools/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1
pheno_dir==${d}/data/plink19_format_01/

sv=$1
svtype=$2
echo "sv=$sv, svtype=$svtype"

# get species name
sp=`grep -w "$sv" ${d}/data/${svtype}_name_conversion_table.txt | cut -f4 | uniq` 
echo "Species name: $sp"
all_nsamples=()

# create the output folder for meta-analysis results
meta_out_dir=${d}/results_fastGWA/${svtype}/meta/${sv}/
meta_out_filebase=${meta_out_dir}/${sv}.meta_res
mkdir -p ${d}/results_fastGWA/${svtype}/meta/${sv}/

# prepare metal script header
metal_script=${d}/scripts/scripts_fastGWA_${svtype}/metal_per_sv/${sv}.metal.txt
cat ${d}/scripts/scripts_fastGWA_${svtype}/metal_header.txt > $metal_script
for p in `seq 1 $nperm`
do
   cat ${d}/scripts/scripts_fastGWA_${svtype}/metal_header.txt > ${d}/scripts/scripts_fastGWA_${svtype}/metal_per_sv/${sv}.metal.perm${p}.txt
done

# check in which cohorts the SV is present
IFS=',' read -ra cohorts_with_sv <<< `grep -w $sv ${d}/data/${svtype}_per_cohort.txt | cut -f6`
cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
echo "Cohorts with SV: $cohorts_joined"


#
# Primary (real) GWAS:
#



for cohort in ${cohorts_with_sv[@]}
do
    echo -e "\n\nRunning the analysis for ${cohort}\n\n"

    res_dir=${d}/results_fastGWA/${svtype}/${cohort}/${sv}/

    geno_file=${d}/genotypes/${cohort}/${cohort}_filtered
    grm=${d}/genotypes/${cohort}/GCTA/GRM_${cohort}_sparse
    gender_file=/groups/umcg-lifelines/tmp01/projects/${cohort}_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/${cohort}/${cohort}_gender.txt 
    
    #get the column number for the SV in the SV table
    col=`head -1 ${pheno_dir}/${cohort}.${svtype}.filtered.plink19.txt | sed "s:\t:\n:g" | tail -n+3 | grep -w -n ${sv} | cut -d ":" -f1`

    # get the species column numbers from the species abundance file 
    col_cov=`head -1 ${pheno_dir}/${cohort}.covariates.plink19.txt | sed "s:\t:\n:g"  | grep -w -n ${sp} | cut -d ":" -f1`

    # write the selected quantitative covariates to a tmp file (species abundance, age, read number)
    awk -v c=$col_cov 'BEGIN {FS=OFS="\t"}; {print $1, $2, $c, $(NF-1), $(NF) }' ${pheno_dir}/${cohort}.covariates.plink19.noheader.txt  \
    >  ${res_dir}/tmp.qcovar.txt

    #
    # run real GWAS analysis
    #
    $gcta \
      --bfile /groups/umcg-lifelines/tmp01/projects/${cohort}_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/genotypes/${cohort}/${cohort}_filtered \
      --grm-sparse ${grm} \
      --fastGWA-mlm \
      --pheno ${pheno_dir}/${cohort}.${svtype}.filtered.plink19.noheader.txt \
      --mpheno $col \
      --qcovar ${res_dir}/tmp.qcovar.z.txt \
      --covar ${gender_file} \
      --out ${res_dir}/${sv}

    
    #mkfifo ${res_dir}/${svtype}.${cohort}.${sv}
    
    echo "$sv real analysis plink return code: $?"    
    
    # Number of samples with association results for this SV for this cohort
    n=`head -2 ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.${f_ext} | tail -1 | awk '{print $9}'`
    all_nsamples+=( $n )
    
    gzip -f ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.${f_ext}
    
    # append the per cohort result location to the metal script
    echo -e "PROCESS\t${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.${f_ext}.gz\n" >> $metal_script


done



# Real GWAS meta-analysis
echo -e "OUTFILE\t${meta_out_filebase} .tbl\nANALYZE\nQUIT" >> $metal_script
metal $metal_script
echo "${sv}, real analysis metal return code: $?"

for cohort in ${cohorts_with_sv[@]}
do 
    res_dir=${d}/results_fastGWA/${svtype}/${cohort}/${sv}/
    rm ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.${f_ext}.gz
done

# Convert to EMP
cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
samplesize_joined=`printf -v var '%s,' "${all_nsamples[@]}"; echo "${var%,}"`

tail -n+2 ${meta_out_filebase}1.tbl | \
sort -k6g | \
python3 ${script_dir}/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined 0.05 | tail -n+2  | gzip -c \
> ${meta_out_filebase}.eQTLs.txt.gz

tail -n+2  ${meta_out_filebase}1.tbl | sort -k6g | awk -v c=${cohorts_joined} -v s=${samplesize_joined} -v svname=${sv} 'BEGIN {FS=OFS="\t"}; {print svname,$0, c, s}' | gzip -c \
> ${meta_out_filebase}.annot.tbl.gz

rm ${meta_out_filebase}1.tbl*

