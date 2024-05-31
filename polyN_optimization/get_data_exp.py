import numpy as np
import pandas as pd
import os

file_dir = os.path.dirname(os.path.abspath(__file__)) 
sep = os.sep
material = "DP780"

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
        

def export_exp_data(material, thetas):

    ut_data = pd.DataFrame(columns=["q", "LoadAngle", "YieldStress", "Rval", "Type", "YoungMod", "Width", "Thickness"])
    foldername_in = f"{file_dir}{sep}results_exp{sep}{material}"

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
            filepath = foldername_in + sep + filename

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
            filepath = foldername_in + sep + filename

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
            filepath = foldername_in + sep + filename

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

    foldername_out = f"{file_dir}{sep}calibration_data{sep}{material}"
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    print(foldername_out)
    ut_data.to_csv(foldername_out + sep + 'data_exp_{}.csv'.format(material), index=False)

def export_exp_data_2(material):
    mat_dic = analyze_exp_data(material)
    ut_data = pd.DataFrame(columns=["q", "LoadAngle", "YieldStress", "Rval", "Type", "YoungMod", "Width", "Thickness"])
    dirname_in = f"{file_dir}{sep}results_exp{sep}{material}"

    for test in mat_dic.keys():
        test_dic = mat_dic[test]
        if test == "UT":

            for ori in test_dic.keys():
                n_exp = test_dic[ori]
                ys_exp = np.array([])
                rval_exp = np.array([])

                if ori == "EBT":
                    for i in range(n_exp):
                        filename = f"{test}_{ori}_{i+1}.csv"
                        filepath = dirname_in + sep + filename
                        df = pd.read_csv(filepath)
                        df = df.rename(columns={df.columns[0]:"PlasticStrain", df.columns[1]:"PlasticStress[MPa]"})
                        yield_stress = df[df["PlasticStrain"] > 0.002]["PlasticStress[MPa]"].iloc[0]
                        ys_exp = np.append(ys_exp,yield_stress)
                        rval_exp = np.append(rval_exp,0)
                    q = 1
                    theta = 0.0  
                else:
                    for i in range(n_exp):
                        filename = f"{test}_{ori}_{i+1}.csv"
                        filepath = dirname_in + sep + filename
                        df = pd.read_csv(filepath)
                        yield_stress = df[df["PlasticStrain_longi"] > 0.002]["PlasticStress[MPa]"].iloc[0]
                        rval = df["R-Value"].iloc[0]
                        ys_exp = np.append(ys_exp,yield_stress)
                        rval_exp = np.append(rval_exp,rval)
                    q = 0.0
                    theta = float(ori)
                rval = np.mean(rval_exp)
                ys = np.mean(ys_exp)
                pttype = "e"
                row = pd.DataFrame({"q": [q], "LoadAngle" : [theta], "YieldStress" : [ys], "Rval" : [rval], "Type" : [pttype]})
                ut_data = pd.concat([ut_data, row])

    foldername_out = f"{file_dir}{sep}calibration_data{sep}{material}"
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    print(foldername_out)
    ut_data.to_csv(foldername_out + sep + 'data_exp_{}.csv'.format(material), index=False)

export_exp_data_2(material)