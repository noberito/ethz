import numpy as np
import pandas as pd
import os

file_dir = os.path.dirname(os.path.abspath(__file__)) 
dir = os.sep


thetas = ["00", "15", "30", "45", "60", "75", "90"]
materials = ["AA7020-T6", "DP600", "DP780"]

def export_exp_data(material, thetas):
    ut_data = pd.DataFrame(columns=["q", "LoadAngle", "YieldStress", "Rval", "Type", "YoungMod", "Width", "Thickness"])
    foldername_in = f"{file_dir}{dir}results_exp{dir}{material}"

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
            filepath = foldername_in + dir + filename

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

    exp = "UT_EBT"
    ys_exp = np.array([])

    if material == "DP780" :
        nexp = 2
        for i in range(nexp):
            test = exp + "_" + str(i + 1)
            filename = test + ".csv"
            filepath = foldername_in + dir + filename

            db = pd.read_csv(filepath)
            db = db.rename(columns={db.columns[0]:"PlasticStrain", db.columns[1]:"PlasticStress[MPa]"})
            yield_stress = db[db["PlasticStrain"] > 0.002]["PlasticStress[MPa]"].iloc[0]
                
            ys_exp = np.append(ys_exp,yield_stress)

        ys = np.mean(ys_exp)
        rval = 0
        q = 1
        pttype = "e"
        theta = 0.0

    else:
        nexp=1
        for i in range(nexp):
            test = exp + "_" + str(i + 1)
            filename = test + ".csv"
            filepath = foldername_in + dir + filename

            db = pd.read_csv(filepath)
            db = db.rename(columns={db.columns[0]:"PlasticStrain", db.columns[1]:"PlasticStress[MPa]"})
            yield_stress = db[db["PlasticStrain"] > 0.002]["PlasticStress[MPa]"].iloc[0]
                
            ys_exp = np.append(ys_exp,yield_stress)

        ys = np.mean(ys_exp)
        rval = 0
        q = 1
        pttype = "e"
        theta = 0.0

    row = pd.DataFrame({"q": [q], "LoadAngle" : [theta], "YieldStress" : [ys], "Rval" : [rval], "Type" : [pttype]})
    ut_data = pd.concat([ut_data, row])

    foldername_out = f"{file_dir}{dir}calibration_data{dir}{material}"
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    print(foldername_out)
    ut_data.to_csv(foldername_out + dir + 'data_exp_{}.csv'.format(material), index=False)
