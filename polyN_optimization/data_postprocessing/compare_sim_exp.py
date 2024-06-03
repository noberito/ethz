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

p = read_param()
material = p["material"]

postprocess_dir = os.path.dirname(os.path.abspath(__file__))
polyN_cali_dir = os.path.dirname(postprocess_dir)


tests = ["UT_00", "UT_15", "UT_30", "UT_45", "UT_60", "UT_75","UT_90", "UT_EBT"]
""""CH_00", "CH_45",
"NT6_00", "NT6_45", "NT6_90",
"NT20_00", "NT20_45",
"SH_000", "SH_090", "SH_p45"]"""

def compare(tests, material):

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    n = len(tests)

    fig, ax = plt.subplots(n // 3 + 1, 3)

    for i in range(n):
        test = tests[i]
        j = i // 3
        k = i % 3

        exp_res_path = results_exp_dir + sep + test + "_1.csv"
        sim_res_path = results_sim_dir + sep + test + "_polyN.csv"

        df_sim = pd.read_csv(sim_res_path)
        df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S13", "S23"]], axis=1)
        df_sim["E"] = np.linalg.norm(df_sim[["E11", "E22", "E33", "E12", "E13", "E23"]], axis=1)
        df_exp = pd.read_csv(exp_res_path)
        if test == "UT_EBT":
            df_exp = df_exp.rename(columns={df_exp.columns[0]:"TrueStrain", df_exp.columns[1]:"TrueStress[MPa]"})
        ax[j,k].plot(df_sim["E"], df_sim["S"], marker = "x",label = "sim")
        ax[j,k].plot(df_exp["TrueStrain"], df_exp["TrueStress[MPa]"], marker = "x", label = "exp")   
        ax[j,k].set_title(test)
        ax[j,k].legend()

    plt.legend()
    plt.tight_layout()
    plt.show()

compare(tests, material)