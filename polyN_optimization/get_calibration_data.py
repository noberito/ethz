import numpy as np
import pandas as pd
import os

file_dir = os.path.dirname(os.path.abspath(__file__)) 
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
        

def export_exp_data_old(material, thetas):
    """
        Not used
    """

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

def tresca(S):
    """
    Vectorized Tresca function
    Input :
        - S : np.array of shape (N,6) or (6,), stress
    Output :
        - res : np.array of shape (N, 1)"""
    
    S = np.atleast_2d(S)
    N = len(S)
    if S.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    S_matrix = np.zeros((N, 3, 3))
    
    index_x = np.array([0, 1, 2, 0, 0, 1])
    index_y = np.array([0, 1, 2, 1, 2, 2])
    S_matrix[:, index_x, index_y] = S
    S_matrix = S_matrix + np.transpose(S_matrix, (0, 2, 1))
    S_matrix[:,np.arange(3),np.arange(3)] = S_matrix[:,np.arange(3),np.arange(3)] // 2

    ps = np.linalg.eigvals(S_matrix)
    dps = np.zeros((N, 3))

    dps[:,0] = np.abs(ps[:, 0] - ps[:, 1])
    dps[:,1] = np.abs(ps[:, 0] - ps[:, 2])
    dps[:,2] = np.abs(ps[:, 1] - ps[:, 2])

    return(np.max(dps, axis=1))

def get_coeff_hill48(material):
    pass

def hill48(sigma, coeff_mat):
    """
    Vectorized Hill48 function for the given coefficients
    Input :
        - S : np.array of shape (N,6) or (6,), stress
        - coeff_mat : np.array of shape (6,), coefficients of the hill48 function for the material studied
    Output :
        - res : np.array of shape (N, 1)"""
    
    sigma = np.atleast_2d(sigma)
    if sigma.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    s11 = sigma[:, 0]
    s22 = sigma[:, 1]
    s33 = sigma[:, 2]
    s12 = sigma[:, 3]
    s13 = sigma[:, 4]
    s23 = sigma[:, 5]

    F = coeff_mat[0]
    G = coeff_mat[0]
    H = coeff_mat[0]
    L = coeff_mat[0]
    M = coeff_mat[0]
    N = coeff_mat[0]

    res = np.sqrt(F * (s22 - s33)**2 + G * (s33 - s11)**2 + H * (s11 - s22)**2 + 2 * L * s23**2 + 2 * M * s13**2 + 2 * N * s12**2)
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
    if protomodel == "tresca":
        data = data_yf_sphere(tresca, itermax, nb_virtual_pt)
    if protomodel == "Hill48":
        coeff_hill = get_coeff_hill48(material)
        def hill48_mat(sigma):
            return(hill48(sigma, coeff_hill))
        data = data_yf_sphere(hill48_mat, itermax, nb_virtual_pt)
        pass
    if protomodel == "Yld2000":
        pass

    df = pd.DataFrame(data, columns=["s11", "s22", "s33", "s12", "s13", "s23"])
    df["Rval"] = np.zeros(len(data))
    df["Type"] = ["v"] * len(data)

    filename = f"data_virtual_{material}_{protomodel}.csv"
    folderpath = f"{file_dir}{sep}calibration_data{sep}{material}"
    filepath = folderpath + sep + filename

    if not os.path.exists(folderpath):
        os.makedirs(folderpath)
    df.to_csv(filepath, index=False)
