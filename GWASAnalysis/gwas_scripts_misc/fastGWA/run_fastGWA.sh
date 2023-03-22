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
svtype=$2
echo "sv=$sv, svtype=$svtype"

# get species name
sp=`grep -w "$sv" ${d}/data_fastGWA/${svtype}_name_conversion_table.txt | cut -f4 | uniq` 
echo "Species name: $sp"
all_nsamples=()
all_cohorts=()
# create the output folder for meta-analysis results
meta_out_dir=${d}/results_fastGWA/${svtype}/meta/${sv}/
meta_out_filebase=${meta_out_dir}/${sv}.meta_res
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
      --out ${res_dir}/${sv}
    
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

# Convert to EMP
cohorts_joined=`printf -v var '%s,' "${all_cohorts[@]}"; echo "${var%,}"`
samplesize_joined=`printf -v var '%s,' "${all_nsamples[@]}"; echo "${var%,}"`

#tail -n+2 ${meta_out_filebase}1.tbl | \
#sort -k6g | \
#python3 ${script_dir}/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined 0.05 | tail -n+2  | gzip -c \
#> ${meta_out_filebase}.eQTLs.txt.gz

tail -n+2  ${meta_out_filebase}1.tbl | sort -k6,6g | awk -v c=${cohorts_joined} -v s=${samplesize_joined} -v svname=${sv} 'BEGIN {FS=OFS="\t"}; {print svname,$0, c, s}' | gzip -cf \
> ${meta_out_filebase}.annot.tbl.gz

zcat ${meta_out_filebase}.annot.tbl.gz | awk '{FS=OFS="\t"}; {if ($7 < 5e-8) print}' | gzip -cf > ${meta_out_filebase}.annot.5e-8.tbl.gz


# add per cohort Z and P
cmd=""
for cohort in ${all_cohorts[@]}
do
    cmd="$cmd | python3 ${script_dir}/add_columns_from_file_v2.py -i stdin -i_m 1 -f_m 1 -f_cols 7,9 -f ${d}/results_fastGWA/${svtype}/${cohort}/${sv}/${sv}.fastGWA.gz"
done
full_cmd="zcat ${meta_out_filebase}.annot.5e-8.tbl.gz $cmd | ${script_dir}/fastGWA/postprocess_vSV_gwas.py | gzip -c  > ${meta_out_filebase}.annot.5e-8.per_cohort.tbl.gz"
eval $full_cmd

for cohort in ${cohorts_with_sv[@]}
do 
    res_dir=${d}/results_fastGWA/${svtype}/${cohort}/${sv}/
    rm ${res_dir}/${sv}.fastGWA.gz
done

rm ${meta_out_filebase}1.tbl*

