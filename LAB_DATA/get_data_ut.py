import numpy as np
import pandas as pd
import os

exp = "UT_00"
nexp = "1"
material = "DP600"

filename = exp + "_" + nexp + ".csv"
foldername = material + "_results/DATA/"
filepath = "./" + foldername + filename

db = pd.read_csv(filepath)
db["PlasticStrain"] = np.square(db["PlasticStrain_longi"]) + np.square(db["PlasticStrain_thick"]) + np.square(db["PlasticStrain_width"])

yield_stress = db[db["PlasticStrain_longi"] > 0.002]["Time[s]"]
print(yield_stress)


