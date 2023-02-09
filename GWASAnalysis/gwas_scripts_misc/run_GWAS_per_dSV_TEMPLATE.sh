#!/bin/bash
#SBATCH --job-name=SV___BACLIST__
#SBATCH --output=logs/run_GWAS_SV___BACLIST__.out
#SBATCH --error=logs/run_GWAS_SV___BACLIST__.err
#SBATCH --time=90:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load PLINK
module load Metal

svtype=dSV
d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir="${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/"
nperm=10


while read line
do
    sv=$line
    echo $sv
    
    # get species name
    sp=`grep -w "$sv" ${d}/data/${svtype}_name_conversion_table.txt | cut -f4 | uniq` 
    all_nsamples=()

	# create the output folder for meta-analysis results
    meta_out_dir=${d}/results/${svtype}/meta/${sv}/
    meta_out_filebase=${meta_out_dir}/${sv}.meta_res
    mkdir -p ${d}/results/${svtype}/meta/${sv}/

	# prepare metal script header
    metal_script=${d}/scripts/scripts_${svtype}/metal_per_sv/${sv}.metal.txt
    cat ${d}/scripts/scripts_${svtype}/metal_header.txt > $metal_script
    for p in `seq 1 $nperm`
    do
       cat ${d}/scripts/scripts_${svtype}/metal_header.txt > ${d}/scripts/scripts_${svtype}/metal_per_sv/${sv}.metal.perm${p}.txt
    done

    # check in which cohorts the SV is present
    IFS=',' read -ra cohorts_with_sv <<< `grep -w $sv ${d}/data/${svtype}_per_cohort.txt | cut -f6`

	# run the GWAS per cohort with SV
    for cohort in ${cohorts_with_sv[@]}
    do
        echo -e "\n\nRunning the analysis for ${cohort}\n\n"

        res_dir=${d}/results/${svtype}/${cohort}/${sv}/
        mkdir -p ${res_dir}/permutations/

        geno_file=${d}/genotypes/${cohort}/${cohort}_filtered
        perm_fam_dir=${d}/genotypes/${cohort}/permuted/
        pheno_file=${d}/data/${cohort}.${svtype}.filtered.txt
        covar_file=${d}/data/${cohort}.covariates.txt   
        
        covars="age,read_number,$sp,PC1,PC2"
        if [ $cohort == "DAG3" ]
        then
            echo "DAG3! Use PCs 1-5"
            covars="age,read_number,$sp,PC1,PC2,PC3,PC4,PC5"
        fi

        #
        # run real GWAS analysis
        #
        plink2 \
            --glm sex hide-covar \
            --bfile ${geno_file} \
            --pheno ${pheno_file} \
            --pheno-name "${sv}" \
            --covar ${covar_file} \
            --covar-name "${covars}" \
            --covar-variance-standardize \
            --out ${res_dir}/${svtype}.${cohort}.${sv} 
        
        echo "$sv real analysis plink return code: $?"
        
        # Number of samples with association results for this SV for this cohort
        n=`head -2 ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.logistic | tail -1 | awk '{print $8}'`
        all_nsamples+=( $n )
        
        # format the assoc v2/results for METAL: add A1 and A2
        awk '{OFS="\t"}; {if ($6 == $4) {oa=$5}; if ($6 == $5) {oa=$4}; if (NR == 1) {oa="A2"}; {print $1,$2,$3,$6, oa, $7,$8,$9,$10,$11,$12}}' \
        ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.logistic | gzip -c > ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.logistic.gz
        reformat_returncode=$?
        
        rm ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.logistic
        
        # append the per cohort result location to the metal script
        echo -e "PROCESS\t${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.logistic.gz\n" >> $metal_script
        
        
        
        #
        # run permuted GWAS
        #
        for i in `seq 1 $nperm`
        do
            plink2 \
                --glm sex hide-covar \
                --bed ${geno_file}.bed \
                --bim ${geno_file}.bim \
                --fam ${perm_fam_dir}/perm${i}.fam \
                --pheno ${pheno_file} \
                --pheno-name "${sv}" \
                --covar ${covar_file} \
                --covar-name "${covars}" \
                --covar-variance-standardize \
                --out ${res_dir}/permutations/${svtype}.${cohort}.${sv}.perm${i}
                

            echo "$sv permutation $i plink return code: $?"
            
            # format the assoc v2/results for METAL: add A1 and A2
            awk '{OFS="\t"}; {if ($6 == $4) {oa=$5}; if ($6 == $5) {oa=$4}; if (NR == 1) {oa="A2"}; {print $1,$2,$3,$6, oa, $7,$8,$9,$10,$11,$12}}'  \
            ${res_dir}/permutations/${svtype}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic | gzip -c > ${res_dir}/permutations/${svtype}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic.gz
            reformat_returncode=$?
            
            rm ${res_dir}/permutations/${svtype}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic 
            
            # append the per cohort result location to the metal script
            echo -e "PROCESS\t${res_dir}/permutations/${svtype}.${cohort}.${sv}.perm${i}.${sv}.glm.logistic.gz" >> ${d}/scripts/scripts_${svtype}/metal_per_sv/${sv}.metal.perm${i}.txt
        done
    done


    #
    # Run meta-analysis
    #
    
    # Real GWAS meta-analysis
    echo -e "OUTFILE\t${meta_out_filebase} .tbl\nANALYZE\nQUIT" >> $metal_script
    metal $metal_script
    echo "${sv}, real analysis metal return code: $?"

	# Permuted GWAS meta-analysis
    for i in `seq 1 $nperm`
    do
        echo -e "OUTFILE\t${meta_out_filebase}.perm${i} .tbl\nANALYZE\nQUIT" >>  ${d}/scripts/scripts_${svtype}/metal_per_sv/${sv}.metal.perm${i}.txt
        metal ${d}/scripts/scripts_${svtype}/metal_per_sv/${sv}.metal.perm${i}.txt
        echo "${sv}, permutation $i metal return code: $?"
    done
    

    # remove per cohort results
    for cohort in ${cohorts_with_sv[@]}
    do 
        res_dir=${d}/results/${svtype}/${cohort}/${sv}/
        rm ${res_dir}/${svtype}.${cohort}.${sv}.${sv}.glm.logistic.gz
        rm ${res_dir}/permutations/${svtype}.${cohort}.${sv}.perm*.${sv}.glm.logistic.gz
    done


	# format the meta-analysis results for FDR calculation by EMP
    cohorts_joined=`printf -v var '%s,' "${cohorts_with_sv[@]}"; echo "${var%,}"`
    samplesize_joined=`printf -v var '%s,' "${all_nsamples[@]}"; echo "${var%,}"`
    
    tail -n+2 ${meta_out_filebase}1.tbl | \
    sort -k6g | \
    python3 ${script_dir}/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined 0.05 | tail -n+2  | gzip -c \
    > ${meta_out_filebase}.eQTLs.txt.gz
    
    tail -n+2  ${meta_out_filebase}1.tbl | sort -k6g | awk -v c=${cohorts_joined} -v s=${samplesize_joined} -v svname=${sv} 'BEGIN {FS=OFS="\t"}; {print svname,$0, c, s}' | gzip -c \
    > ${meta_out_filebase}.annot.tbl.gz
    
    rm ${meta_out_filebase}1.tbl

    
    for p in `seq 1 $nperm`
    do
        tail -n+2 ${meta_out_filebase}.perm${p}1.tbl | \
        sort -k6g | \
        python3 ${script_dir}/metal_to_EMP.py stdin ${sv} $cohorts_joined $samplesize_joined 0.05 | tail -n+2 | gzip -c \
        > ${meta_out_filebase}.eQTLs.perm${p}.txt.gz
        rm ${meta_out_filebase}.perm${p}1.tbl
    done
    
done < __BACLIST__