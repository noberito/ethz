import numpy as np
import pandas as pd
import os
from bezier import bezier_3D

polyN_dir = os.path.dirname(os.path.abspath(__file__)) 
sep = os.sep

tol_yf = 0.01
itermax = 20


def analyze_exp_data(material):
    """
        Returns a dictionnary of the different tests led for the given material /(in the results_exp/material folder, files
        not to read must start with "_")
        Input :
            - material : string
        Output :
            - d : dictionnary (d[test][orientation][number])

    """
    d = {}
    folder = polyN_dir + sep + "results_exp" + sep + material
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

def export_exp_data(material):
    """
        Generate the csv file data_exp_{material} needed to optimize.
        For the moment, only the UT files are checked: Must have a column "PlasticStrain_longi" and "PlasticStress[MPa]"
        and if available a column "R-value" containing the r-val
        Input :
            - material : string
    """
    mat_dic = analyze_exp_data(material)
    ut_data = pd.DataFrame(columns=["q", "LoadAngle", "YieldStress", "Rval", "Type", "YoungMod", "Width", "Thickness"])
    dirname_in = f"{polyN_dir}{sep}results_exp{sep}{material}"

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
                        yield_stress = df[df["PlasticStrain"] > 0.06]["PlasticStress[MPa]"].iloc[0]
                        ys_exp = np.append(ys_exp,yield_stress)
                        rval_exp = np.append(rval_exp,0)
                    q = 1
                    theta = 0.0  
                else:
                    for i in range(n_exp):
                        filename = f"{test}_{ori}_{i+1}.csv"
                        filepath = dirname_in + sep + filename
                        df = pd.read_csv(filepath)
                        yield_stress = df[df["PlasticStrain_longi"] > 0.06]["PlasticStress[MPa]"].iloc[0]
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

    ut_data = ut_data.sort_values(by=['q', 'LoadAngle'], ascending=[True, True])
    foldername_out = f"{polyN_dir}{sep}calibration_data{sep}{material}"
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    print(foldername_out)
    ut_data.to_csv(foldername_out + sep + 'data_exp_{}.csv'.format(material), index=False)

def mises(S):
    """
    Vectorized Mises function
    Input :
        - S : np.array of shape (N,6) or (6,), stress
    Output :
        - res : np.array of shape (N, 1)
    """
    S = np.atleast_2d(S)
    if S.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    s11 = S[:, 0]
    s22 = S[:, 1]
    s33 = S[:, 2]
    s12 = S[:, 3]
    s13 = S[:, 4]
    s23 = S[:, 5]

    res = np.sqrt(0.5 * (s22 - s33)**2 + 0.5 * (s33 - s11)**2 + 0.5 * (s11 - s22)**2 + 3 * (s23**2 + s13**2 + s12**2))
    return res


def generate_dir(nb_virtual_pt, ndim):
    """
        Generate random nb_virtual_pt directions in dimension ndim
        Input :
            - nb_virtual_pt : integer
            - ndim : integer
        Output :
            - u : ndarray of shape (nb_virtual_pt, ndim)
    """
    u = np.random.normal(0, 1, (nb_virtual_pt, ndim))
    return(u)

def data_yf_sphere(f, itermax, nb_pt_check):
    """
        For a given function f, returns nb_pt_check points on the surface of equation f(S) = 1
        Input :
            - f : nd.array of shape (6,) -> float
            - itermax : integer, number of maximum iterations in the trial and error
            - nb_pt_check : integer
        Output :
            - data : ndarray of shape (nb_pt_check, 6)
    """
    data = np.zeros((nb_pt_check, 6))
    j = 0

    while j < nb_pt_check:

        E = np.zeros(6)
        m = -1
        u = generate_dir(1, 6)[0]
        it = 0
        lamb = 1

        while m < tol_yf and it < itermax :

            E = E + lamb * u
            m = f(E) - 1
            lamb = lamb * 2
            it = it + 1

        if it != itermax :

            S = np.zeros(6)
            res = f(E) - 1
            it = 0
            while abs(res) > tol_yf:
                I = (S+E)/2
                res = f(I) - 1
                if res > 0 :
                    E = I
                else:
                    S = I
                it = it + 1
            
            data[j] = I
            j = j + 1
    
    return(data)

def export_virtual_data(protomodel, material, nb_virtual_pt):
    """
        Generate the csv file data_exp_{material} needed to optimize the polyN function.
        Based on the protomodel given.
        Input :
            - protomodel : string, must be in ["mises", "tresca", "Hill48", "Yld2000"]
            - material : string
            - nb_virtual_pt : integer
    """

    if protomodel == "mises":
        data = data_yf_sphere(mises, itermax, nb_virtual_pt)

    if protomodel == "bezier":
        data = bezier_3D(material, nb_virtual_pt)

    df = pd.DataFrame(data, columns=["s11", "s22", "s33", "s12", "s13", "s23"])
    df["Rval"] = np.zeros(len(data))
    df["Type"] = ["v"] * len(data)

    filename = f"data_virtual_{material}_{protomodel}.csv"
    folderpath = f"{polyN_dir}{sep}calibration_data{sep}{material}"
    filepath = folderpath + sep + filename

    if not os.path.exists(folderpath):
        os.makedirs(folderpath)
    df.to_csv(filepath, index=False)
