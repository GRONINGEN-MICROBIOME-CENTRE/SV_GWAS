#!/usr/bin/env bash
#
# Annotates the association results with rs id, gene names etc
#

in_f=$1 # eQTLs.txt.gz file path
f=$2 # output file prefix
svtype=$3
basedir="/data/umcg-tifn/SV/SV_GWAS/"
script_dir="/home/umcg-dzhernakova/scripts/umcg_scripts/SV_GWAS/clean/"
col=3
genes=${base_dir}resources/ensembl_b37_genes.bed
genes_prot=${base_dir}resources/ensembl_b37_protcoding_genes.bed

module load BEDTools
module load Pysam

# Extract all associations with P < 5e-08
zcat eQTLs.txt.gz | awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1 || $1 < 5e-8) print }' | cut -f1-5,9-12,14 > ${f}.txt

# Add rs ids from dbSNP
python ${script_dir}/gwas_scripts_misc/add_rs_by_position.py \
${f}.txt ${base_dir}/resources/All_20180423.vcf.gz 1 \
> ${f}.rsids.txt

# Add overlapping genes
awk -v c=${col} '{FS=OFS="\t"}; {split($c, snp, ":"); if (NR != 1) {print snp[1], snp[2] - 1, snp[2], $c}}' $f > ${f}.rsids.tmp.bed
echo -e "SNPCoord\tOverlapping_genes" > ${f}.rsids.tmp.genes
intersectBed -a ${f}.rsids.tmp.bed -b $genes -wa -wb | \
  python2.7 ${script_dir}/gwas_scripts_misc/collapseIntersectBedRes.py - \
  >> ${f}.rsids.tmp.genes

 python2.7 ${script_dir}/gwas_scripts_misc/add_columns_from_file.py  \
  -i ${f} --header -f ${f}.rsids.tmp.genes -i_m 2 -f_m 0 -f_cols 1 \
  > ${f}.rsids.tmp.genes1.txt

# Add protein-coding genes within a 250kb window from the SNP
echo -e "SNPCoord\tProtcoding_genes_within_250kb" > ${f}.rsids.tmp.genes2
windowBed -a ${f}.rsids.tmp.bed -b $genes_prot -w 250000 | \
  python2.7 ${script_dir}/gwas_scripts_misc/collapseIntersectBedRes.py - \
  >> ${f}.rsids.tmp.genes2

 python2.7 ${script_dir}/gwas_scripts_misc/add_columns_from_file.py  \
  -i ${f}.rsids.tmp.genes1.txt -f ${f}.rsids.tmp.genes2 -i_m 2 -f_m 0 -f_cols 1 \
  > ${f%txt}genes.txt


  rm ${f}.rsids.tmp*

# Add per-cohort z-scores and p-values and heterogeneity p-value. Fix!
python ${script_dir}/gwas_scripts_misc/add_per_cohort_z_p.py ${f} ${svtype} \
 >  ${f}.rsids.tmp.hetero

# Add the original SV id and SNP-associated traits from the GWAS catalog
# TODO: FDR significant p-value - determine the p-value threshold automatically from the other file
python2.7 ${script_dir}/gwas_scripts_misc/add_columns_from_file.py \
-i ${f%txt}genes.txt  \
-f ${base_dir}data/${svtype}_name_conversion_table.txt \
-i_m 5 -f_m 2 -f_cols 1,0 --header |  \
python2.7 ${script_dir}/gwas_scripts_misc/add_columns_from_file.py \
-i stdin \
-f ${base_dir}resources/gwas_catalog_v1.0.2-associations_e105_r2021-12-21.cut.collated.txt \
-i_m 1 -f_m 0 -f_cols 1 --header \
|  awk 'BEGIN {FS=OFS="\t"}; {if (NR == 1) print $0, "FDR<0.05"; else if ($1 <= 5.179e-09) print $0, "1"; else print $0, "0"} ' |
 python2.7 ${script_dir}/gwas_scripts_misc/add_columns_from_file.py -i stdin -i_m 5 -f ${base_dir}data/${svtype}_annotation.txt -f_m 9 -f_cols 7,8 --header \
> ${f%txt}genes.gwas_annot.txt

# Clump the results:
Rscript ${script_dir}/gwas_scripts_misc/clumping_and_pleiotropy.R ${f%txt}genes.gwas_annot.txt
