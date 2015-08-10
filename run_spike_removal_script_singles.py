import sys
sys.path.append('/global/project/projectdirs/planck/software/zonca/software/testenv')
from testenv import cluster
horns=['18','19','20','21']
if sys.argv[1] != 'all':
    horns=[sys.argv[1]]
diodes=['00','01','10','11']
for horn in horns:
	for diode in diodes:
		cmd = ["python SPIKE_correct_all_script.py %s %s" % (horn,diode)]
		cluster.run_serial("spikrem_%s_%s" % (horn,diode), "\n".join(cmd), mem=20, hours=4)
