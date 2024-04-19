import numpy as np
import pandas as pd
import os
import openpyxl

thetas = ["00", "15", "30", "45", "60", "75", "90"]
materials = ["AA7020-T6", "DP600", "DP780"]

for material in materials:
    ut_data = pd.DataFrame(columns=["q", "LoadAngle", "YieldStress", "Rval", "Type"])
    foldername = material + "_results/DATA/"

    if material == "AA7020-T6" :
        nexp = 2
    else :
        nexp = 3

    for theta in thetas:
        exp = "UT_" + theta
        ys_exp = np.array([])
        rval_exp = np.array([])

        for i in range(nexp):
            test = exp + "_" + str(i + 1)
            filename = test + ".csv"
            filepath = "./" + foldername + filename

            db = pd.read_csv(filepath)
            yield_stress = db[db["PlasticStrain_longi"] > 0.002]["PlasticStress[MPa]"].iloc[0]
            rval = db["R-Value"].iloc[0]
            
            ys_exp = np.append(ys_exp,yield_stress)
            rval_exp = np.append(rval_exp,rval)

        ys = np.mean(ys_exp)
        rval = np.mean(rval_exp)
        q = 0.0
        pttype = "e"
        theta = int(theta)

        row = pd.DataFrame({"q": [q], "LoadAngle" : [theta], "YieldStress" : [ys], "Rval" : [rval], "Type" : [pttype]})
        ut_data = pd.concat([ut_data, row])

    ut_data.to_csv('UT_cal_{}.csv'.format(material), index=False)