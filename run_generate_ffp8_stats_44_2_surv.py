import matplotlib
matplotlib.use("agg")
import sys
sys.path.append('/global/project/projectdirs/planck/software/zonca/software/testenv')
from testenv import cluster
    

s1list=['full','full','full','full']
s2list=['s1','s2','s5','s6']
s3list=['s3','s4','s7','s8']

hr1list=['null','null','null','null']
hr2list=['null','null','null','null']
q1list=['null','null','null','null']
q2list=['null','null','null','null']
freqs=['044']
for freq in freqs:
    for q1,q2,s1,s2,hr1,hr2,s3 in zip(q1list,q2list,s1list,s2list,hr1list,hr2list,s3list):
        cmd = [
            "python generate_ffp8_cls_2_cmd.py %s %s %s %s %s %s %s %s" % (freq,q1,s1,hr1,q2,s2,hr2,s3)
            ]
        cluster.run_serial("cl_44_%s_%s" % (s2,s3), "\n".join(cmd), mem=20, hours=8)