import subprocess
import os
import glob
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np
import time
from optimize_polyN_cluster import first_opti

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

    run_large = f"python run_large_sim_cluster.py"
    subprocess.call(run_large, shell=True, cwd=exec_dir)

    time.sleep(20)

    lck_files = glob.glob(exec_dir + sep + "*.lck")
    while len(lck_files) > 0:
        t = time.time() - t0
        print(t)
        time.sleep(5)
        lck_files = glob.glob(exec_dir + sep + "*.lck")
    



framework()