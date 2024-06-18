import subprocess
import os
import glob
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np
import time
from optimize_polyN import first_opti

sep = os.sep
polyN_cali_dir = os.path.dirname(os.path.abspath(__file__))
exec_dir = polyN_cali_dir + sep + "running"
sys.path.append(polyN_cali_dir)

from read_param import read_param

p = read_param()



def framework():

    t0 = time.time()

    #first_optimization = "python optimize_polyN_cluster.py"
    #subprocess.call(first_optimization, shell=True, cwd=polyN_cali_dir)
    first_opti()

    



framework()