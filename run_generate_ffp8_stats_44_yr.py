import matplotlib
matplotlib.use("agg")
import sys
sys.path.append('/global/project/projectdirs/planck/software/zonca/software/testenv')
from testenv import cluster

#years
plotdir='/global/homes/p/peterm/notebooks/plotting/nulls/years/'
s1list=['full','full','full','full','full','yr1+yr3','full','full']
s2list=['full','yr1','yr2','yr3','yr4','yr2+yr4','yr1+yr3','yr2+yr4']
hr1list=['_hr1','null','null','null','null','null','null','null']
hr2list=['_hr2','null','null','null','null','null','null','null']
q1list=['null','null','null','null','null','null','null','null']
q2list=['null','null','null','null','null','null','null','null']
freqs=['044']
for freq in freqs:
    for q1,q2,s1,s2,hr1,hr2 in zip(q1list,q2list,s1list,s2list,hr1list,hr2list):
        cmd = [
            "python generate_ffp8_cls_cmd.py %s %s %s %s %s %s %s" % (freq,q1,s1,hr1,q2,s2,hr2)
            ]
        cluster.run_serial("stats%s%s%s" % (freq,s1,s2), "\n".join(cmd), mem=20)