import numpy as np
import pandas as pd
import os

material = "DP780"
file_dir = os.path.dirname(os.path.abspath(__file__))
sep = os.sep

def analyze_exp_data(material):
    d = {}
    folder = file_dir + sep + "results_exp" + sep + material
    files = os.listdir(folder)
    for f in files:
        if not(f[0] == "_"):
            test, ori, n = f.strip(".csv").split("_")
            if d.get(test)==None:
                d[test] = {}
            if d[test].get(ori) == None :
                d[test][ori] = 1
            else:
                d[test][ori] = d[test][ori] + 1
    return(d)
        

analyze_exp_data(material)