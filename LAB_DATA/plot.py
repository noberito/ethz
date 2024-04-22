import numpy as np
import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

exp = "UT_00"
nexp = "1"
material = "DP600"

filename = exp + "_" + nexp + ".csv"
foldername = material + "_results/DATA/"
filepath = "./" + foldername + filename

sns.set_theme()
sns.set_style('whitegrid')

db = pd.read_csv(filepath)
print(db.columns)
db["PlasticStrain"] = np.square(db["PlasticStrain_longi"]) + np.square(db["PlasticStrain_thick"]) + np.square(db["PlasticStrain_width"])
sns.lineplot(data=db, x="PlasticStrain_longi", y="TrueStress[MPa]")
plt.show()