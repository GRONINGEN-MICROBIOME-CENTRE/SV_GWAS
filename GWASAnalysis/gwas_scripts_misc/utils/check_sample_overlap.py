import sys
import pandas as pd

f1 = pd.read_csv(sys.argv[1], sep = "\t", header = None)
col1 = int(sys.argv[2])
ids1 = list(f1.iloc[:, col1])

f2 = pd.read_csv(sys.argv[3], delimiter = '[\s\t]+', header = None, engine='python')
col2 = int(sys.argv[4])
ids2 = list(f2.iloc[:, col2])

f3 = pd.read_csv(sys.argv[5], delimiter = '[\s\t]+', header = None, engine='python')
col3 = int(sys.argv[6])
ids3 = list(f3.iloc[:, col3])

if len(sys.argv) > 6:
	print(sys.argv[7])


cnt = [0,0,0,0]

for id in ids1:
    res = [id, '0', '0']
    cnt[0] += 1
    if id in ids2:
        res[1] = '1'
        cnt[1] += 1
    if id in ids3:
        res[2] = '1'
        cnt[2] += 1
    if res[2] == '1' and res[1] == '1':
        cnt[3] += 1
    print ("\t".join(res))

sys.stderr.write(" --- ".join(map(str, cnt)) + "\n")
