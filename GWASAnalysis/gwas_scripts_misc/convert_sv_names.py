import sys
sv_name_fname = sys.argv[1]
annot_fname = sys.argv[2]

sv_sp_dict = {}
# Read the species name for each SV from the annotation file
with open(annot_fname) as annot:
	for l in annot:
        spl = l.rstrip().split("\t")
        sv_sp_dict[spl[0]] = spl[1]

# Make the new SV id by using species:N and replacing special chars like / and -
with open(sv_name_fname) as f:
    print("sv_id\tsv_species\tnew_sv_id\tnew_sv_species")
    cur_sp = ""
    idx = 0
    for l in f:
        sv = l.rstrip()
        sp = ":".join(sv.split(":")[:len(sv.split(":"))-1])
        sp_short = sv_sp_dict[sp]
        if cur_sp == sp_short:
            idx = idx + 1
        else:
            idx = 1
            cur_sp = sp_short
        sp_short = sp_short.replace("-","_").replace("/", "_").replace(" ", "_")
        new_id = sp_short + ":" + str(idx)
        print (sv + "\t" + sp + "\t" + new_id + "\t" + sp_short)