import numpy as np
import pandas as pd
import os
import openpyxl

thetas = ["00", "15", "30", "45", "60", "75", "90"]
materials = ["AA7020-T6", "DP600", "DP780"]

for material in materials:
    foldername = material + "_results/DATA/"

    if material == "AA7020-T6" :
        nexp = 2
    else :
        nexp = 3

    ut_yf_data = pd.DataFrame()
    ut_rval_data = pd.DataFrame()

    for theta in thetas:
        exp = "UT_" + theta
        ys_exp = []
        rval_exp = []
        for i in range(nexp):
            
            test = exp + "_" + str(i + 1)
            filename = test + ".csv"
            filepath = "./" + foldername + filename

            db = pd.read_csv(filepath)
            yield_stress = db[db["PlasticStrain_longi"] > 0.002]["PlasticStress[MPa]"].iloc[0]
            rval = db["R-Value"].iloc[0]
            
            ys_exp.append(yield_stress)
            rval_exp.append(rval)

        ut_yf_data[exp] = ys_exp
        ut_rval_data[exp] = rval_exp


    with pd.ExcelWriter('UT_exp_{}.xlsx'.format(material)) as writer:
        ut_yf_data.to_excel(writer, sheet_name='Yield stress (MPa)', index=False)
        ut_rval_data.to_excel(writer, sheet_name='R-Values', index=False)