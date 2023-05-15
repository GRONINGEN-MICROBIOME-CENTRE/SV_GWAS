#!/usr/bin/python
"""
Updates the file header according to the conversion table provided as the second argument.
The optional 3rd and 4th arguments are the column numbers in the conversion file to get the column number containing the old col name and the column number with the new column name.
By default the first column contains the old name, the second - the new name.
"""
import sys
import gzip

in_fname = sys.argv[1]
rename_fname = sys.argv[2]
col_old = 0
col_new = 1
if len(sys.argv) > 3:
    col_old = int(sys.argv[3])
    col_new = int(sys.argv[4])

if in_fname != "stdin":
    if in_fname.endswith(".gz"):
        in_file = gzip.open(in_fname)
    else:
        in_file = open(in_fname)
else:
    in_file = sys.stdin


# fill the name conversion dict from the second file
rename_dict = {}
with open (rename_fname) as f:
    for l in f:
        spl = l.rstrip().split("\t")
        rename_dict[spl[col_old]] = spl[col_new]

# Replace the old header with the new values
out_str = ""
header = in_file.readline().rstrip().split("\t")
for h in header:
    h2 = h
    if len(h) > 0:
        if h in rename_dict.keys():
             h2 = rename_dict[h]
        else:
             sys.stderr.write("WARNING! No such value in the renaming file: " + h + ". The old column name will be used!\n")
    out_str += "\t" + h2
print(out_str.replace("\t", "", 1))

for l in in_file:
    print(l.rstrip())
