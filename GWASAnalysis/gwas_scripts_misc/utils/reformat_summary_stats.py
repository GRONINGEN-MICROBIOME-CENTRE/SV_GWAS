#
# Reformats GWAS summary statistics for upload
#


import gzip
from collections import defaultdict
import sys

species = sys.argv[1]
svtype = sys.argv[2]
basepath = sys.argv[3]

annot_fname = basepath + "/data/" + svtype + "_annotation.txt"
annot_dict = {}
with open(annot_fname, mode="rt") as annot_f:
    annot_f.readline()
    for l in annot_f:
        spl = l.rstrip().split("\t")
        annot_dict[spl[9]] = spl[3]


def get_meta_maf(spl, freq_res, cohorts_dict):
	res_freq_lst = ["NA"]*4
	snp_cohorts = spl[8].split(",")
	cohorts_used = spl[7]
	
	AC = 0
	AN = 0
	
	for i, c in enumerate(snp_cohorts):
		if not cohorts_used[i] == '?':
			c_num = cohorts_dict[c]
			c_res = freq_res[c_num]
			if not c_res == 'NA':
				freq_res_cohort = align_alleles(spl, c_res)
				AC += freq_res_cohort[0]
				AN += sum(freq_res_cohort)
	return (AC/AN)


# TODO: replace DAG3 with DMP

def align_alleles(spl, c_res):
    ea = spl[2]
    oa = spl[3]
    if ea in c_res and oa in c_res:
        return((c_res[ea], c_res[oa]))
    elif not isAT_GCsnp(ea,oa) and complement(ea) in c_res and complement(oa) in c_res:
        return((c_res[complement(ea)], c_res[complement(oa)]))
    else:
        print ("\t".join(spl) + "\t" +"\t".join(c_res.keys()))
    return(("NA", "NA"))
    #elif ea not in c_res and oa not in c_res:
    #    if 


def complement(nucl):
    if nucl == 'A':
        return 'T'
    if nucl == 'T':
        return 'A'
    if nucl == 'G':
        return 'C'
    if nucl == 'C':
        return 'G'
 
def isAT_GCsnp(ref, alt):
    if ref == 'T' and alt == 'A':
        return True
    if ref == 'A' and alt == 'T':
        return True 
    if ref == 'G' and alt == 'C':
        return True
    if ref == 'C' and alt == 'G':
        return True
    return False



freq_dict = {}
cohorts = ["DAG3", "LLD", "300OB", "500FG"]
cohorts_dict = {name : num for num, name in enumerate(cohorts)}
for c_num, c in enumerate(cohorts):
    with gzip.open(basepath + "/genotypes/" + c + "/" + c + "_filtered.frq.counts.gz", mode="rt") as frq_f:
        frq_f.readline()
        cnt = 0
        for l in frq_f:
            spl = l.strip().split()
 
            res_dict = {spl[2] : int(spl[4]), spl[3] : int(spl[5])}
            if spl[1] in freq_dict:
                freq_dict[spl[1]][c_num] = res_dict
            else:
                res_lst = ["NA"]*4
                res_lst[c_num] = res_dict
                freq_dict[spl[1]] = res_lst


out_f = open(basepath + "/results_all_summary_stats/" + svtype + "/meta_combined/" + species+ "/" + species + "."+ svtype + ".meta-analysis.annot.txt", "w")
out_f.write("SV_id\tSNP\tchromosome\tposition\teffect_allele\tother_allele\tmeta_samplesize\tmeta_zscore\tmeta_pvalue\tmeta_EAF\tcohorts_with_SV\teffect_direction_per_cohort\tsamplesize_per_cohort\n")
with gzip.open(basepath + "/results_all_summary_stats/" + svtype + "/meta_combined/" + species+ "/" + species+ "." + svtype + ".meta-analysis.txt.gz", mode="rt") as f:
	#cnt = 0
	f.readline()
	for l in f:
		spl = l.rstrip().split()
		spl[2] = spl[2].upper()
		spl[3] = spl[3].upper()
		spl[4] = spl[4].replace(".00", "")
		sv_id = annot_dict[spl[0]]
		spl[0] = sv_id
		AF = get_meta_maf(spl, freq_dict[spl[1]], cohorts_dict)
	
		chrom, pos = spl[1].split(":")
		spl[1] = spl[1] + "\t" + chrom + "\t" + pos
		cohorts = spl[8]
		spl[8] = spl[7]
		spl[7] = cohorts
	
		_ = out_f.write ("{0}\t{1:.6f}\t{2}\n".format("\t".join(spl[:7]), AF, "\t".join(spl[7:])))

out_f.close()

