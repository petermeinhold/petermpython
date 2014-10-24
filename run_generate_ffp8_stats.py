from testenv import cluster
    

s1list=['full','full','full','full','full','full','full','full']
s2list=['s1','s2','s3','s4','s5','s6','s7','s8']
hr1list=['','','','','','','','']
hr2list=['','','','','','','','']
q1list=['','','','','','','','']
q2list=['','','','','','','','']
freqs=['070']
for freq in freqs:
    for q1,q2,s1,s2,hr1,hr2 in zip(q1list,q2list,s1list,s2list,hr1list,hr2list):
        cmd = [
            "python generate_ffp8_cls_cmd.py %s %s %s %s %s %s %s" % (freq,q1,s1,hr1,q2,s2,hr2)
            ]
        cluster.run_serial("pipe_%s_%s_%s_%s_%s_%s_%s" % (freq,q1, s1,hr1,q2,s2,hr2), "\n".join(cmd), mem=5)