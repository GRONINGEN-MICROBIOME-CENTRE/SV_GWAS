import os
import sys
import subprocess

script_dir = os.path.abspath(os.path.dirname(__file__)) 
svtype = sys.argv[2]

with open(sys.argv[1]) as f:
    print(f.readline().rstrip() + "\tsame_direction\tbetas_per_cohort\tP_per_cohort\tHetero_P")
    for l in f:
        spl = l.rstrip().split("\t")
        snp = spl[2]
        sv = spl[0]
        directions = [*spl[8]]
        datasets = spl[9].split(",")
        Ns = spl[10].split(",")
        
        #Filter results tested in > 1 cohort, check if the effect direction is concordant
        new_directions = []
        new_datasets = []
        new_Ns = []
        for i, direction in enumerate(directions):
            if not direction == "?":
                new_directions.append(direction)
                new_datasets.append(datasets[i])
                new_Ns.append(Ns[i])
        if not(len(new_directions) == len(new_datasets) == len(new_Ns)):
            print("ERROR! Resulting lists are not of equal length! " + " ".join(spl))
        
        if len(new_directions) > 1:
            spl[7] = "".join(new_directions)
            spl[8] = ",".join(new_datasets)
            spl[9] = ",".join(new_Ns)
            if new_directions.count(new_directions[0]) == len(new_directions):
                concordant = "1"
            else:
                concordant = "0"
        else:
            continue
        if concordant == "0":
            res = 3*['NA']
        else:
            # Add per cohort beta, P, heterogeneity P
            output = subprocess.check_output(['sh', script_dir + '/lookup_SV.sh', sv, snp, svtype, spl[8]])
            res = output.decode().split("\n")[-2].split(" ")
        
            if new_directions[0] == '+':
                new_betas = ",".join(map(str,(map(abs, map(float,res[0].split(","))))))
            elif new_directions[0] == '-':
                new_betas = ",".join(map(str,[-1*i for i in map(abs, map(float,res[0].split(",")))]))
            else:
                "ERROR! Wrong direction char"
                new_betas = "NA"
            res[0] = new_betas
        print("\t".join(spl) + "\t" + concordant + "\t" + "\t".join(res))
        
