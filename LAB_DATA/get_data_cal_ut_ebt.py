import numpy as np
import pandas as pd
import os
import openpyxl

nexp = 1
 
materials = ["AA7020-T6", "DP600"]

for material in materials:
    ut_data = pd.DataFrame(columns=["q", "LoadAngle", "YieldStress", "Rval", "Type"])
    foldername = material + "_results/DATA/"


    exp = "UT_EBT"
    ys_exp = np.array([])

    for i in range(nexp):
        test = exp + "_" + str(i + 1)
        filename = test + ".csv"
        filepath = "./" + foldername + filename

        db = pd.read_csv(filepath, names=["PlasticStrain_longi", "PlasticStress[MPa]"])
        yield_stress = db[db["PlasticStrain_longi"] > 0.002]["PlasticStress[MPa]"].iloc[0]
            
        ys_exp = np.append(ys_exp,yield_stress)

    ys = np.mean(ys_exp)
    rval = "*"
    q = 1
    pttype = "e"
    theta = 0.0

    row = pd.DataFrame({"q": [q], "LoadAngle" : [theta], "YieldStress" : [ys], "Rval" : [rval], "Type" : [pttype]})
    ut_data = pd.concat([ut_data, row])

    ut_data.to_csv(foldername + 'UT_EBT_cal_{}.csv'.format(material), index=False)