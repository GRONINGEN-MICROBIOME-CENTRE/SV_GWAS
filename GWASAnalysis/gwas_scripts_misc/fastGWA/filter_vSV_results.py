import os
import sys
import subprocess

script_dir = os.path.abspath(os.path.dirname(__file__)) 

with open(sys.argv[1]) as f:
    print(f.readline().rstrip() + "\tsame_direction")
    for l in f:
        spl = l.rstrip().split("\t")
        snp = spl[1]
        sv = spl[0]
        directions = [*spl[7]]
        datasets = spl[8].split(",")
        Ns = spl[9].split(",")
        
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
        print("\t".join(spl))
        
