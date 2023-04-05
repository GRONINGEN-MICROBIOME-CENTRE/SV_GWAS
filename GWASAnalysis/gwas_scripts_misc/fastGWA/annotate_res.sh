#!/usr/bin/env bash
#
# Annotates the association results with rs id, gene names etc
#

f=$1 # eQTLs.txt.gz file path
svtype=$3

d=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/
res_dir=${d}/results_fastGWA/${svtype}/meta_combined/
script_dir=/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/v2/scripts/SV_GWAS/GWASAnalysis/gwas_scripts_misc/

col=3
genes=${d}/../resources/ensembl_b37_genes.bed
genes_prot=${d}/../resources/ensembl_b37_protcoding_genes.bed

module load BEDTools
module load Pysam

# Add overlapping genes
awk -v c=${col} '{FS=OFS="\t"}; {split($c, snp, ":"); if (NR != 1) {print snp[1], snp[2] - 1, snp[2], $c}}' $f > ${f}.tmp.bed
echo -e "SNP\tOverlapping_genes" > ${f}.tmp.genes
intersectBed -a ${f}.tmp.bed -b $genes -wa -wb | \
  python2.7 ${script_dir}/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes

 python3 ${script_dir}/add_columns_from_file_v2.py  \
  -i ${f} --header -f ${f}.tmp.genes -i_m 2 -f_m 0 -f_cols 1 \
  > ${f}.tmp.genes1.txt

# Add protein-coding genes within a 250kb window from the SNP
echo -e "SNP\tProtcoding_genes_within_250kb" > ${f}.tmp.genes2
windowBed -a ${f}.tmp.bed -b $genes_prot -w 250000 | \

python2.7 ${script_dir}/collapseIntersectBedRes.py - \
  >> ${f}.tmp.genes2

python3 ${script_dir}/add_columns_from_file_v2.py  \
  -i ${f}.tmp.genes1.txt -f ${f}.tmp.genes2 -i_m 2 -f_m 0 -f_cols 1 \
  > ${f%txt}genes.txt


rm ${f}.tmp*


# Add the original SV id and SNP-associated traits from the GWAS catalog

python3 ${script_dir}/add_columns_from_file_v2.py \
-i ${f%txt}genes.txt  \
-f ${d}/data_fastGWA/${svtype}_name_conversion_table.txt \
-i_m 0 -f_m 2 -f_cols 1,0 --header |  \
python3 ${script_dir}/add_columns_from_file_v2.py \
-i stdin \
-f ${d}/../resources/gwas_catalog_v1.0.2-associations_e105_r2021-12-21.cut.collated.txt \
-i_m 1 -f_m 0 -f_cols 1 --header \
 | python3 ${script_dir}/add_columns_from_file_v2.py -i stdin -i_m 21 -f ${d}/data_fastGWA/${svtype}_annotation.txt -f_m 3 -f_cols 12,13 --header \
> ${f%txt}genes.gwas_annot.txt
