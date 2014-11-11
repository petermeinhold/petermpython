import matplotlib
matplotlib.use("agg")
import sys
sys.path.append('/global/project/projectdirs/planck/software/zonca/software/testenv')
from testenv import cluster

#years
plotdir='/global/homes/p/peterm/notebooks/plotting/nulls/years/'
s1list=['full','full','full']
s2list=['full','full','full']
hr1list=['null','null','null']
hr2list=['null','null','null']
q1list=['_18_23','_18_23','_19_22']
q2list=['_19_22','_20_21','_20_21']
freqs=['070']
for freq in freqs:
    for q1,q2,s1,s2,hr1,hr2 in zip(q1list,q2list,s1list,s2list,hr1list,hr2list):
        cmd = [
            "python generate_ffp8_cls_cmd.py %s %s %s %s %s %s %s" % (freq,q1,s1,hr1,q2,s2,hr2)
            ]
        print (cmd)
        cluster.run_serial("mcstats%s%s%s" % (freq,q1,q2), "\n".join(cmd), mem=20)