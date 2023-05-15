#!/usr/bin/python

import sys
from collections import defaultdict

if len(sys.argv) > 1:
	if sys.argv[1] != "-":
		f = open(sys.argv[1])
	else:
		f = sys.stdin
else:
    f = sys.stdin

map=defaultdict(set)
for line in f:
	spl = line.strip().split("\t")
	map[spl[3]].add(spl[7])

for snp, g in map.items():
    print (snp + "\t" + ",".join(g))

