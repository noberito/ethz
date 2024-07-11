#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import sys
from get_calibration_data import analyze_exp_data

# In[2]:


file_dir = os.path.dirname(os.path.abspath(__file__))
sep = os.sep


# In[3]:


# In[4]:


"""--------------------------------------------EXTRACTING DATA TO CALIBRATE PLASTIC RULE AND YOUND MODULUS------------------------------------------------------------"""

n_ut_00_mat = {"DP780" : 3, "DP600" : 3, "AA7020-T6" : 2}


def get_hardening_law(material):
    """
        Generate the csv file "data_plasticlaw_{material}" used to calibrate the hardening law.
        Based on the UT_00 tests from the given material
        Input :
            - material : string
    """
    
    mat_tests = analyze_exp_data(material)
    n_ut_00 = mat_tests["UT"]["00"]
    df_out = pd.DataFrame(columns = ["PlasticStrain", "PlasticStress", "YoungModulus"])
    ymod = 0
    for i in range(1, n_ut_00 + 1):

        filename_in = "UT_00_{}.csv".format(i)
        foldername_in = file_dir + sep + "results_exp" + sep + material + sep

        filepath = foldername_in + filename_in
        
        df = pd.read_csv(filepath, index_col=False)

        ymod = ymod + df[" Young's Modulus [MPa]"][0]
        df = df[df["PlasticStrain_longi"] > 0.002]

        df["PlasticStrain"] = df["PlasticStrain_longi"] - df["PlasticStrain_longi"].iloc[0]
        df["PlasticStress"] = df["PlasticStress[MPa]"]
        df["YoungModulus"] = 0
        df_out = pd.concat([df_out, df[["PlasticStrain", "PlasticStress", "YoungModulus"]]])

    df_out.iloc[0, df_out.columns.get_loc('YoungModulus')] = ymod / n_ut_00

    foldername_out = file_dir + sep + "calibration_data" + sep + material
    filename_out = foldername_out + sep + f"data_plasticlaw_{material}.csv"
    df_out.to_csv(filename_out, index=False)

