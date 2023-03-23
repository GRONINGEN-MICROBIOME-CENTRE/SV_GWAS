import sys
import gzip

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
    #Ns = datasets

    new_directions = []
    new_datasets = []
    new_Ns = []
    betas = []
    pvals = []
    for i, direction in enumerate(directions):
        if not direction == "?":
            new_directions.append(direction)
            new_datasets.append(datasets[i])
            new_Ns.append(Ns[i])
            if direction == '+':
                betas.append(abs(float(spl[15 + 2*i])))
            elif direction == '-':
                betas.append(-1*abs(float(spl[15 + 2*i])))
            pvals.append(float(spl[15 + 2*i + 1]))
    if not(len(new_directions) == len(new_datasets) == len(new_Ns)):
        print("ERROR! Resulting lists are not of equal length! " + " ".join(spl))
    
    if len(new_directions) > 1:
        spl[7] = "".join(new_directions)
        spl[13] = ",".join(new_datasets)
        spl[14] = ",".join(new_Ns)
        if new_directions.count(new_directions[0]) == len(new_directions):
            concordant = "1"
        else:
            concordant = "0"
    else:
        continue
    print("\t".join(spl[:8]) + "\t" + "\t".join(spl[11:15]) + "\t" + ",".join(map(str,betas)) + "\t" + ",".join(map(str,pvals)) + "\t" + concordant)
