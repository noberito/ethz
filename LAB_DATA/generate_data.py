import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

tol_yf = 0.01
itermax = 20

current_dir = "./"  # Assuming current directory
dir = "/"
if os.name == "nt":
    current_dir = ".\\"
    dir = "\\"

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


def export_virtual_data(proto_yf, protomodel, material, nb_virtual_pt):

    data = data_yf_sphere(proto_yf, itermax, nb_virtual_pt)
    db = pd.DataFrame(data, columns=["s11", "s22", "s33", "s12", "s13", "s23"])
    db["Type"] = ["v"] * len(data)

    filename = "data_virtual_" + material + "_" + protomodel + ".csv"
    filepath = current_dir + material + "_results" + dir + "DATA" + dir + filename

    db.to_csv(filepath, index=False)