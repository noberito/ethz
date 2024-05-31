import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

tol_yf = 0.01
itermax = 20

file_dir = os.path.dirname(os.path.abspath(__file__))  # Assuming current directory
dir = os.sep

def mises(sigma):
    """
    Vectorized Mises function
    """
    sigma = np.atleast_2d(sigma)
    if sigma.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    s11 = sigma[:, 0]
    s22 = sigma[:, 1]
    s33 = sigma[:, 2]
    s12 = sigma[:, 3]
    s13 = sigma[:, 4]
    s23 = sigma[:, 5]

    res = np.sqrt(0.5 * (s22 - s33)**2 + 0.5 * (s33 - s11)**2 + 0.5 * (s11 - s22)**2 + 3 * (s23**2 + s13**2 + s12**2))
    return res

def tresca(sigma):
    sigma = np.atleast_2d(sigma)
    N = len(sigma)
    if sigma.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    sigma_matrix = np.zeros((N, 3, 3))
    
    index_x = np.array([0, 1, 2, 0, 0, 1])
    index_y = np.array([0, 1, 2, 1, 2, 2])
    sigma_matrix[:, index_x, index_y] = sigma
    sigma_matrix = sigma_matrix + np.transpose(sigma_matrix, (0, 2, 1))
    sigma_matrix[:,np.arange(3),np.arange(3)] = sigma_matrix[:,np.arange(3),np.arange(3)] // 2

    ps = np.linalg.eigvals(sigma_matrix)
    dps = np.zeros((N, 3))

    dps[:,0] = np.abs(ps[:, 0] - ps[:, 1])
    dps[:,1] = np.abs(ps[:, 0] - ps[:, 2])
    dps[:,2] = np.abs(ps[:, 1] - ps[:, 2])

    return(np.max(dps, axis=1))

def generate_dir(nb_virtual_pt, ndim):
    u = np.random.normal(0, 1, (nb_virtual_pt, ndim))
    return(u)

def data_yf_sphere(f, itermax, nb_pt_check):
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

    if protomodel == "mises":
        data = data_yf_sphere(mises, itermax, nb_virtual_pt)
    if protomodel == "tresca":
        data = data_yf_sphere(tresca, itermax, nb_virtual_pt)
    df = pd.DataFrame(data, columns=["s11", "s22", "s33", "s12", "s13", "s23"])
    df["Rval"] = np.zeros(len(data))
    df["Type"] = ["v"] * len(data)

    filename = f"data_virtual_{material}_{protomodel}.csv"
    folderpath = f"{file_dir}{dir}calibration_data{dir}{material}"
    filepath = folderpath + dir + filename

    if not os.path.exists(folderpath):
        os.makedirs(folderpath)
    df.to_csv(filepath, index=False)