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

snp2gene = defaultdict(set)
for line in f:
	spl = line.strip().split("\t")
	snp2gene[spl[3]].add(spl[7])

for snp, g in snp2gene.items():
    print (snp + "\t" + ",".join(g))

