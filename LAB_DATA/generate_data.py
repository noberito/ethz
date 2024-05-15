import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import sklearn

current_dir = "./"  # Assuming current directory
dir = "/"
if os.name == "nt":
    current_dir = ".\\"
    dir = "\\"

def generate_dir(nb_virtual_pt):
    u = np.random.normal(0, 1, (nb_virtual_pt,6))
    return(u)

#plot_implicit(mises_plane)
#TODO ca peut planter etre faux ici si 0 n'est pas compris dans l'intervalle des valeurs prises
def generate_data(proto_yf, nb_virtual_pt, check_valid_data=False, precision=10000):
    np.random.seed(98)
    data = np.zeros((nb_virtual_pt, 6))
    us = generate_dir(nb_virtual_pt)
    alpha = np.linspace(0, 30, precision)
    alpha = alpha[np.newaxis, :]
    for i in range(nb_virtual_pt):
        if i % 1000 == 0:
            print(i)
        u = np.expand_dims(us[i], axis=1)
        sigmas = np.dot(u, alpha).T
        yss = proto_yf(sigmas) - np.ones(precision)
        k = np.argmin(np.abs(yss))
        data[i] = sigmas[k]
    
    if check_valid_data:
        yf_res = proto_yf(data)
        print("Le max est : ", np.max(yf_res), ", Le min est : ", np.min(yf_res))

        tol = 0.01
        print("Le nombre de points où l'écart à 1 est supérieur à ", tol," est : ", np.count_nonzero(np.where(np.abs(yf_res - np.ones(nb_virtual_pt)) > tol)))
    return(data)


def export_virtual_data(proto_yf, protomodel, material, nb_virtual_pt):

    data = generate_data(proto_yf, nb_virtual_pt)

    db = pd.DataFrame(data, columns=["s11", "s22", "s33", "s12", "s13", "s23"])
    db["Type"] = ["v"] * len(data)

    filename = "data_virtual_" + material + "_" + protomodel + ".csv"
    filepath = current_dir + material + "_results" + dir + "DATA" + dir + filename

    db.to_csv(filepath, index=False)