svtype="dSV"

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
res_dir=${d}/results_fastGWA/${svtype}/meta_combined/
script_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/


#
# 1. merge per SV results
#

mkdir ${res_dir}
cd ${res_dir}


gunzip -f ../meta/*/*.meta_res.annot.5e-8.per_cohort.tbl.gz

echo -e "SV\tSNP\tEffect_allele\tOther_allele\tEffect_size\tSE\tPvalue\tDirection\tHet_Pvalue\tN\tCohorts\tN_per_cohort\tBeta_per_cohort\tPvalue_per_cohort\tEffect_direction_concordant" \
> ${svtype}.fastGWA.5e-08.txt

sort -m -k7,7g -T $TMPDIR -S 70G --buffer-size=1000 ../meta/*/*.meta_res.annot.5e-8.per_cohort.tbl  >> ${svtype}.fastGWA.5e-08.txt

gzip ../meta/*/*.meta_res.annot.5e-8.per_cohort.tbl

#
# 2. Add rs ids
#
python3 ${script_dir}/utils/add_rs_by_position.py \
    ${svtype}.fastGWA.5e-08.txt \
    ${d}../resources/All_20180423.vcf.gz \
    1 \
    > ${svtype}.fastGWA.5e-08.rsids.txt
#
# 3. Filter vSV results (exclude results with heterogeneity P < 0.05, with discordant directions between cohorts, with at least one )
#

awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1 || ($10 > 0.05 && $17 == 1)) print }' ${svtype}.fastGWA.5e-08.rsids.txt  > ${svtype}.fastGWA.5e-08.rsids.flt.txt


#
# 4. Clump
#
Rscript ${script_dir}/utils/clump.R ${svtype}.fastGWA.5e-08.rsids.flt.txt


#
# 5. AnnotATE
#
f=${svtype}.fastGWA.5e-08.rsids.flt.txt.clumped_0.1.txt

col=3
genes=${d}/../resources/ensembl_b37_genes.bed
genes_prot=${d}/../resources/ensembl_b37_protcoding_genes.bed

module load BEDTools


# Add overlapping genes
awk -v c=${col} '{FS=OFS="\t"}; {split($c, snp, ":"); if (NR != 1) {print snp[1], snp[2] - 1, snp[2], $c}}' $f > ${f}.tmp.bed
echo -e "SNP\tOverlapping_genes" > ${f}.tmp.genes
intersectBed -a ${f}.tmp.bed -b $genes -wa -wb | \
  python2.7 ${script_dir}/utils/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes

 python3 ${script_dir}/utils/add_columns_from_file_v2.py  \
  -i ${f} --header -f ${f}.tmp.genes -i_m 2 -f_m 0 -f_cols 1 \
  > ${f}.tmp.genes1.txt

# Add protein-coding genes within a 250kb window from the SNP
echo -e "SNP\tProtcoding_genes_within_250kb" > ${f}.tmp.genes2
windowBed -a ${f}.tmp.bed -b $genes_prot -w 250000 | \
python2.7 ${script_dir}/utils/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes2

python3 ${script_dir}/utils/add_columns_from_file_v2.py  \
  -i ${f}.tmp.genes1.txt -f ${f}.tmp.genes2 -i_m 2 -f_m 0 -f_cols 1 \
  > ${f%txt}genes.txt


rm ${f}.tmp*


# Add the original SV id and SNP-associated traits from the GWAS catalog

python3 ${script_dir}/utils/add_columns_from_file_v2.py \
-i ${f%txt}genes.txt  \
-f ${d}/data_fastGWA/${svtype}_name_conversion_table.txt \
-i_m 0 -f_m 2 -f_cols 1,0 --header |  \
python3 ${script_dir}/utils/add_columns_from_file_v2.py \
-i stdin \
-f ${d}/../resources/gwas_catalog_v1.0.2-associations_e105_r2021-12-21.cut.collated.txt \
-i_m 1 -f_m 0 -f_cols 1 --header \
 | python3 ${script_dir}/utils/add_columns_from_file_v2.py -i stdin -i_m 21 -f ${d}/data_fastGWA/${svtype}_annotation.txt -f_m 3 -f_cols 12,13 --header \
> ${f%txt}genes.gwas_annot.txt
