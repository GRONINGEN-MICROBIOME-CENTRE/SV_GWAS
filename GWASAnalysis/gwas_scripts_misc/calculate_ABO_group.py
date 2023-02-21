#!/usr/local/bin/python
import pandas as pd

cohorts = ["LLD", "500FG", "DAG3"]
base_path = "/data/umcg-tifn/SV/SV_GWAS/"

for c in cohorts:
    fname = base_path + "genotypes/" + c + "/text_genotypes/" + c + ".abo.genotypes.txt"
    data = pd.read_csv(fname, sep = "\t", index_col = 0)
    out = open(base_path + "data/pheno/" + c + ".abo_blood_group.txt", 'w')
    data2 = data.transpose()
    out.write("sampleid\tBloodtype\tBlood_genotype\n")
    for i in range(1,len(data2)):
        geno = data2.iloc[i,:]
        if geno[0] == "T/T":
            abo_geno = "BB"
            abo = "B"
        elif geno[1] == "T/T":
            abo_geno = "O"
            abo = "O"
        elif geno[0] in ["C/T", "T/C"] and geno[1] in ["C/T", "T/C"]:
            abo_geno = "B"
            abo = "B"
        elif geno[0] in ["C/T", "T/C"] and geno[1] == "C/C":
            abo_geno = "AB"
            abo = "AB"
        elif geno[0] == "C/C" and geno[1] in ["C/T", "T/C"]:
            abo_geno = "A"
            abo = "A"
        else:
            abo_geno = "AA"
            abo = "A"
            
        out.write(data2.index[i] + "\t" + abo + "\t" + abo_geno + "\n")
    out.close()
