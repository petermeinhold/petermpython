import sys
sys.path.append('/global/project/projectdirs/planck/software/zonca/software/testenv')
from testenv import cluster
horns=['18','19','20','21','22','23','24','25','26','27','28']
horns=['22']
for horn in horns:
    cmd = ["source activate py27"," python ADC_correct_all_script.py  %s" % (horn)]
    cluster.run_serial("adc_correct_%s" % horn, "\n".join(cmd), mem=20, hours=1)
