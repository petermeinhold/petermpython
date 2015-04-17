import matplotlib
matplotlib.use("agg")
import sys
sys.path.append('/global/project/projectdirs/planck/software/zonca/software/testenv')
from testenv import cluster
    

s1list=['s1','s3','s5','s7']
s2list=['s2','s4','s6','s8']
hr1list=['null','null','null','null']
hr2list=['null','null','null','null']
q1list=['null','null','null','null']
q2list=['null','null','null','null']
freqs=['070']
for freq in freqs:
    for q1,q2,s1,s2,hr1,hr2 in zip(q1list,q2list,s1list,s2list,hr1list,hr2list):
        cmd = [
            "python generate_ffp8_cls_cmd.py %s %s %s %s %s %s %s" % (freq,q1,s1,hr1,q2,s2,hr2)
            ]
        cluster.run_serial("cl70_%s%s" % s1,s2, "\n".join(cmd), mem=20)
