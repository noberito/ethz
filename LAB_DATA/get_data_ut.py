import numpy as np
import pandas as pd
import os

exp = "UT_00"
nexp = "1"
material = "DP600"

def get_data(exp, nexp, material):
    filename = exp + "_" + nexp + ".csv"
    foldername = material + "_results/DATA/"
    filepath = "./" + foldername + filename

    db = pd.read_csv(filepath)
    return(db)

data = get_data(exp, nexp, material)
print(data.columns)

args = 

def data_to_plot(data, args)

