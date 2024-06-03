import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

file_dir = os.path.dirname(os.path.abspath(__file__))  # Assuming current directory
sep = os.sep

polyN_cali_dir = os.path.dirname(file_dir)
sys.path.append(polyN_cali_dir)

from read_param import read_param
from get_calibration_data import analyze_exp_data

p = read_param()
material = p["material"]
degree = int(p["degree"])

def compare(material):

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    ut_tests_mat = analyze_exp_data(material)["UT"]
    n = len(ut_tests_mat.keys())
    nrows = n // 3 + 1

    fig, ax = plt.subplots(nrows, 3)
    i = 0
    for ori in ut_tests_mat.keys():
        j = i // 3
        k = i % 3
        sim_res_path = results_sim_dir + sep + "UT_" + ori + "_polyN.csv"
        df_sim = pd.read_csv(sim_res_path)
        df_sim["S21"] = df_sim["S12"]
        df_sim["S31"] = df_sim["S13"]
        df_sim["S32"] = df_sim["S23"]
        df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S21", "S13", "S31", "S23", "S32"]], axis=1)

        ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S"],label = "abaqus")

        for m in range(ut_tests_mat[ori]):
            exp_res_path = results_exp_dir + sep + "UT_" + ori + f"_{m+1}.csv"
            df_exp = pd.read_csv(exp_res_path)
            if ori == "EBT":
                df_exp = df_exp.rename(columns={df_exp.columns[0]:"TrueStrain", df_exp.columns[1]:"TrueStress[MPa]"})
                e = df_exp["TrueStrain"]
                s = df_exp["TrueStress[MPa]"]
                ax[j,k].plot(e, s, label = f"exp. {m+1}")
            else:
                imax = df_exp["EngStress[MPa]"].idxmax()
                e = df_exp["TrueStrain"].values[:imax]
                s = df_exp["TrueStress[MPa]"].values[:imax]
                ax[j,k].plot(e, s, label = f"exp. {m+1}") 

        ax[j,k].set_title(f"UT_{ori}")
        ax[j,k].set_xlabel(r"$\epsilon$")
        ax[j,k].set_ylabel(r"$\sigma$")
        ax[j,k].grid(1)
        ax[j,k].legend()
        i = i + 1

    for i in range(n, nrows * 3):
        fig.delaxes(ax.flatten()[i])

    fig.suptitle(f"{material} with poly{degree} : Check Experiments vs Abaqus results", fontsize=12)
    plt.legend()
    plt.show()

compare(material)