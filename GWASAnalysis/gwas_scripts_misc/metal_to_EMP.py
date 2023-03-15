"""
Converts the meta-analysis results created by Metal into a eQTL file format suitable for the eQTL mapping pipeline: https://github.com/molgenis/systemsgenetics/wiki/QTL-mapping-pipeline
python3 metal_toEMP.py <metal meta-analysis output> <SV name> <names of cohorts with this SV> <sample sizes of cohorts with this SV> [<first line is header>] [<max p-value to output results>]
Will write no more than 1000000 lines by default
"""


import sys
import gzip

empty_out_spl=["1","x","1","1","x","1","1","trans","A/T","A","0","x","0","0","0","0","x","0","0","0","0","0"]

fname = sys.argv[1]   # metal meta-analysis output
pheno = sys.argv[2]   # SV name
dataset_names = sys.argv[3].split(",")  # names of cohorts with this SV
dataset_samplesizes = sys.argv[4].split(",")   # sample sizes of cohorts with this SV
header = False
pval_threshold = 1
if len(sys.argv) > 5:
    
    if "header" in sys.argv[5]:
        header = True
    else:
        pval_threshold = float(sys.argv[5])

if len(sys.argv) > 6:
    if "header" in sys.argv[6]:
        header = True

max_lines = 1000000

if fname == "stdin":
    f = sys.stdin
else:
    if fname.endswith(".gz"):
        f = gzip.open(fname)
    else:
        f = open(fname)

if header:
    f.readline()
print("PValue\tSNPName\tSNPChr\tSNPChrPos\tProbeName\tProbeChr\tProbeCenterChrPos\tCisTrans\tSNPType\tAlleleAssessed\tOverallZScore\tDatasetsWhereSNPProbePairIsAvailableAndPassesQC\tDatasetsZScores\tDatasetsNrSamples\tIncludedDatasetsMeanProbeExpression\tIncludedDatasetsProbeExpressionVariance\tHGNCName\tIncludedDatasetsCorrelationCoefficient\tMeta-Beta (SE)\tBeta (SE)\tFoldChange\tFDR")

line_num = 1

for l in f:
    spl = l.strip().split()
    snp = spl[0]
    p = float(spl[5])
    if p > pval_threshold:
        continue
    if line_num > max_lines:
        break
    
    datasets = ""
    zscores = ""
    samplesizes = ""
    for i, direction in enumerate(spl[6]):
        if direction == "+":
            datasets += ";" + dataset_names[i]
            samplesizes += ";" + dataset_samplesizes[i]
            zscores += ";1"
        elif direction == "-":
            datasets += ";" + dataset_names[i]
            samplesizes += ";" + dataset_samplesizes[i]
            zscores += ";-1"
    
    # Skip results tested in less than 2 cohorts    
    iif len(datasets.split(";")) < 2:
        continue
    
    zscores = zscores.replace(";","", 1)
    datasets = datasets.replace(";","", 1)
    samplesizes = samplesizes.replace(";","", 1)
    
    zscore = spl[4]
    empty_spl_cp = empty_out_spl[:]
    empty_spl_cp[0] = str(p)
    empty_spl_cp[1] = snp
    empty_spl_cp[2] = snp.split(":")[0]
    empty_spl_cp[3] = snp.split(":")[1]
    empty_spl_cp[4] = pheno
    empty_spl_cp[8] = spl[1].upper() + "/" + spl[2].upper()
    empty_spl_cp[9] = spl[1].upper()
    empty_spl_cp[16] = pheno.replace("\:\d+$", "")
    empty_spl_cp[10] = zscore
    empty_spl_cp[11] = datasets
    empty_spl_cp[12] = zscores
    empty_spl_cp[13] = samplesizes

    print ("\t".join(empty_spl_cp))
    line_num += 1





