import sys
from os import path

result_dir = sys.argv[2]
with open(sys.argv[1]) as sv_list:
    cur_sv = ""
    for l in sv_list:
        sv = l.rstrip()
        sv_resdir = result_dir + "/" + sv
        if not path.exists(sv_resdir + "/" + sv + ".meta_res.annot.tbl.gz") or not path.exists(sv_resdir + "/" + sv + ".meta_res.eQTLs.txt.gz"):
            print (sv + " 0")
            cur_sv = sv
        for i in range(1,11):
            if not path.exists(sv_resdir + "/" + sv + ".meta_res.eQTLs.perm" + str(i) + ".txt.gz"):
                if sv != cur_sv:
                    print (sv + " " + str(i - 1))
                    cur_sv = sv
                print (sv + " " + str(i))
