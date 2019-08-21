import sys
from datetime import datetime

fname = sys.argv[1]
method = fname.split("/")[0]
tissue = fname.split("/")[1].split("_")[0]
filtex = fname.split("/")[1].split("_")[1].split(".")[0]
rtimes = []
with open(fname) as fx:
    for line in fx:
        if line.startswith("Begin PBS"):
            rtimes.append(line.strip()[19:])
s2 = rtimes[-1]
s1 = rtimes[-2]
FMT = '%a %b %d %H:%M:%S %Z %Y'
tdelta = datetime.strptime(s2, FMT) - datetime.strptime(s1, FMT)
print(method, tissue, filtex, len(rtimes), tdelta)
