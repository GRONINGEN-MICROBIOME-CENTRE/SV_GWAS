#!/usr/bin/env bash
#
# Creates new SV table with short names
# 
import os

sv_raw=$1
sv_annot=$2
svtype=$3

#sv_raw="/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_full/20211212_full_deletionStructuralVariation_12388samples.tsv"
#sv_annot="/data/umcg-tifn/SV/profile/20211212_full_v4.0_12388samples_final/SV/SV_info/20211212_full_Informative_species_information.tsv"
#svtype=dSV


d=/data/umcg-tifn/SV/SV_GWAS/data
cur_script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
script_dir=`echo "${cur_script_dir}/gwas_scripts_misc/"`

cd $d

head -1 $sv_raw  | \
sed 's:\t:\n:g' > all_${svtype}_names.txt

sort all_${svtype}_names.txt | sed 's:"::g' > all_${svtype}_names.srt.txt
rm all_${svtype}_names.txt

cut -f9,10 $sv_annot > ${svtype}_names_to_species.txt

python ${script_dir}/convert_sv_names.py all_${svtype}_names.srt.txt ${svtype}_names_to_species.txt > ${svtype}_name_conversion_table.txt
