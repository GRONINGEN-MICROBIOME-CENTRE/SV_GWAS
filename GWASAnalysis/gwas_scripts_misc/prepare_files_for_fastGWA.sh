d="/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/"
script_dir=${d}/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/

mkdir ${d}/data/data_fastGWA
cd ${d}/data/data_fastGWA

svtype=dSV

cohorts=("300OB" "500FG" "LLD" "DAG3")
for cohort in ${cohorts[@]}
do
    # recode phenotypes and recode them to plink 1.9
    if [ $svtype == "dSV" ]
    then
      Rscript ${script_dir}/recode_dSVs.R ../${cohort}.${svtype}.filtered.txt ${cohort}.${svtype}.filtered.txt
    else
      awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print  }' ../${cohort}.${svtype}.filtered.txt > ${cohort}.${svtype}.filtered.txt
    fi
    tail -n+2 ${cohort}.${svtype}.filtered.txt > ${cohort}.${svtype}.filtered.noheader.txt

    # Format covariates in plink 1.9, make sure they include all samples with relatives
    sed 's:"::g' ../abund/${cohort}_abundances.tsv | \
    sed '1s:.:#IID\t&:' | \
    python ${script_dir}/rename_header_based_on_file.py stdin ../${svtype}_name_conversion_table.txt 1 3 | \
    python ${script_dir}/add_columns_from_file.py -i stdin  -f ../pheno/${cohort}_pheno.txt -f_m 0 -f_cols 1,2 | \
    grep -w -v "NA" | \
    awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) $1 = "#FID\tIID"; else $1 = "0\t" $1; print }' \
    > ${cohort}.covariates.txt
    
    Rscript ${script_dir}/scale_read_count.R  ${cohort}.covariates.txt ${cohort}.covariates.scaled.txt
    mv ${cohort}.covariates.scaled.txt ${cohort}.covariates.txt
    tail -n+2 ${cohort}.covariates.txt > ${cohort}.covariates.noheader.txt


done

# make gender files

# make GRMs