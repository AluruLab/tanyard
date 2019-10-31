import sys
from datetime import datetime

TIME_FMT = '%a %b %d %H:%M:%S %Z %Y'
PATTERN = "Begin PBS"

fname = sys.argv[1]
method = fname.split("/")[0]
tissue = fname.split("/")[1].split("_")[0]
filtex = fname.split("/")[1].split("_")[1].split(".")[0]
rtimes = []
with open(fname) as fx:
    for line in fx:
        if line.startswith(PATTERN):
            rtimes.append(line.strip()[19:])
begin_times = [rtimes[i] for i in range(0, len(rtimes), 2)]
end_times = [rtimes[i] for i in range(1, len(rtimes), 2)]
#s2 = rtimes[-1]
#s1 = rtimes[-2]
total_time = None
max_time = None
for s1, s2 in zip(begin_times, end_times):
   tdelta = datetime.strptime(s2, TIME_FMT) - datetime.strptime(s1, TIME_FMT)
   if total_time is None:
       total_time = tdelta
       max_time = tdelta
   else:
       total_time += tdelta
       max_time = max(max_time, tdelta)

print(method, tissue, filtex, len(begin_times), max_time, total_time)
