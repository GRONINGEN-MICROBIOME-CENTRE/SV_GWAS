#
# Runs plink again per cohort for each SNP - SV pair in the top results
#

import subprocess
import sys

sv_type = sys.argv[2]
script_dir = os.path.abspath(os.path.dirname(__file__)) 

with open(sys.argv[1]) as f:
    header = f.readline().rstrip()
    print(header + "\tzscores_per_cohort\tpvalues_per_cohort\theterogeneity_pvalue")
    colnames = {col : i for i,col in enumerate(header.rstrip().split("\t"))}
    snp_col = colnames["SNPCoord"]
    sv_col = colnames["ProbeName"]
    for l in f:
        spl = l.rstrip("\r\n").split("\t")
        snp = spl[snp_col]
        sv = spl[sv_col]
        if sv_type == "dSV":
		output = subprocess.check_output(['sh', script_dir + '/lookup_dSV.sh', sv, snp])
        elif sv_type == "vSV":
		output = subprocess.check_output(['sh', script_dir + '/lookup_vSV.sh', sv, snp])
	else:
		sys.stderr.write("Wrong SV type\n")
	res = output.split("\n")[-2].split(" ")
        spl.extend(res)
        print("\t".join(spl))