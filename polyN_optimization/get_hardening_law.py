#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np


# In[2]:


current_dir = os.getcwd()
dir = os.sep

# In[3]:


# In[4]:


"""--------------------------------------------EXTRACTING DATA TO CALIBRATE PLASTIC RULE AND YOUND MODULUS------------------------------------------------------------"""

n_ut_00_mat = {"DP780" : 3, "DP600" : 3, "AA7020-T6" : 2}


def get_hardening_law(material):
    n_ut_00 = n_ut_00_mat[material]
    df_out = pd.DataFrame(columns = ["PlasticStrain", "PlasticStress", "YoungModulus"])
    ymod = 0
    for i in range(1, n_ut_00 + 1):

        filename_in = "UT_00_{}.csv".format(i)
        foldername_in = current_dir + dir + "results_exp" + dir + material + dir

        filepath = foldername_in + filename_in
        
        df = pd.read_csv(filepath, index_col=False)

        ymod = ymod + df[" Young's Modulus [MPa]"][0]
        df = df[df["PlasticStrain_longi"] > 0.002]

        df["PlasticStrain"] = df["PlasticStrain_longi"] - df["PlasticStrain_longi"].iloc[0]
        df["PlasticStress"] = df["PlasticStress[MPa]"]
        df["YoungModulus"] = 0
        df_out = pd.concat([df_out, df[["PlasticStrain", "PlasticStress", "YoungModulus"]]])

    df_out.iloc[0, df_out.columns.get_loc('YoungModulus')] = ymod / n_ut_00

    foldername_out = current_dir + dir + "calibration_data" + dir + material
    filename_out = foldername_out + dir + "data_plasticlaw.csv"
    df_out.to_csv(filename_out, index=False)

