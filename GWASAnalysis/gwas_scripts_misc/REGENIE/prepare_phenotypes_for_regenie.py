import sys
import pandas as pd

all_svs = pd.read_csv(sys.argv[1], sep = "\t", index_col=0, dtype = 'str')
all_svs = all_svs.replace("1", "0").replace("2", "1")
name_conv = pd.read_csv(sys.argv[2], sep = "\t")
cohort = sys.argv[3]
svtype = sys.argv[4]
outpath = sys.argv[5]

# Select SVs that are called in the cohort
name_conv = name_conv[name_conv['cohorts'].str.contains(cohort)]

# Species : SV dict
name_conv = name_conv.rename(columns={name_conv.columns[0]: 'sv'})
name_conv['species'] = name_conv['sv'].str.replace(":.*", "")

sv_dict = name_conv.groupby('species')['sv'].apply(list).to_dict()

# write SVs of each species in a separate file
for sp,svs in sv_dict.items():
    sp_sv  = all_svs[svs]
    nonmissing = sp_sv.index[~sp_sv.isnull().all(1)]
    sp_sv.loc[nonmissing,:].to_csv(outpath + cohort + "/" + cohort + "." + svtype + "." + sp + ".txt", sep='\t')