import sys
import gzip

"""
Add per cohort P-values, betas, Ns to the meta-analysis results. 
"""


fname = sys.argv[1]

if fname != "stdin":
    if fname.endswith(".gz"):
        f = gzip.open(fname, mode = 'rt')
    else:
        f = open(fname, mode = 'rt')
else:
    f = sys.stdin
    
for l in f:
    spl = l.rstrip().split("\t")
    snp = spl[1]
    sv = spl[0]
    directions = [*spl[7]]
    datasets = spl[13].split(",")
    Ns = spl[14].split(",")

    new_directions = []
    new_datasets = []
    new_Ns = [] # per cohort N
    betas = [] # per cohort betas
    pvals = [] # per cohort P values
    all_cohort_signif_count = 0 # number of cohorts with association P < 0.05
    all_cohort_signif = 0 # is the results nominally significant in at least 2 cohorts?
    for i, direction in enumerate(directions):
        # only keep per-cohort results for tested SNP-SV combinations
        if not direction == "?": 
            new_directions.append(direction)
            new_datasets.append(datasets[i])
            new_Ns.append(Ns[i])
            if direction == '+':
                betas.append(abs(float(spl[15 + 2*i])))
            elif direction == '-':
                betas.append(-1*abs(float(spl[15 + 2*i])))
            pval = float(spl[15 + 2*i + 1])
            pvals.append(pval)
            if pval < 0.05:
                all_cohort_signif_count += 1
    
    if all_cohort_signif_count > 1: all_cohort_signif = 1
      
    if not(len(new_directions) == len(new_datasets) == len(new_Ns)):
        print("ERROR! Resulting lists are not of equal length! " + " ".join(spl))
    
    # require that the association is nominally significant in at least 2 cohorts
    if len(new_directions) > 1:
        spl[7] = "".join(new_directions)
        spl[13] = ",".join(new_datasets)
        spl[14] = ",".join(new_Ns)
        # check if the effect direction is the same in all cohorts
        if new_directions.count(new_directions[0]) == len(new_directions):
            concordant = "1"
        else:
            concordant = "0"
    else:
        continue
    print("\t".join(spl[:8]) + "\t" + "\t".join(spl[11:15]) + "\t" + ",".join(map(str,betas)) + "\t" + ",".join(map(str,pvals)) + "\t" + concordant + "\t" + str(all_cohort_signif))
