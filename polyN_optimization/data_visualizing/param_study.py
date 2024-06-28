import numpy as np
import pandas as pd
import os
import sys
import glob
import matplotlib.pyplot as plt

file_dir = os.path.dirname(os.path.abspath(__file__))
polyN_dir = os.path.dirname(file_dir)
sep = os.sep

sys.path.append(polyN_dir)

from get_calibration_data import analyze_exp_data

def test_evolution(test, input_type, material, var_optim):
    results_sim_dir = polyN_dir + sep + "results_sim" + sep + material + sep + "param_study"
    n_opti = len(glob.glob(results_sim_dir + sep + f"{test}_{input_type}_{var_optim}_*"))
    fig, ax = plt.subplots(1)
    for i in range(n_opti):
        filename = f"{test}_{input_type}_{var_optim}_{i}.csv"
        df = pd.read_csv(results_sim_dir + sep + filename)
        X = df["U2"]
        Y = df["RF2"]
        ax.plot(X, Y, label=f"{i}", linestyle="dashed")
    
    tests_mat = analyze_exp_data(material)
    type_test = test.split("_")[0]
    ori = test.split("_")[1]
    type_tests_mat = tests_mat[type_test]
    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material

    for m in range(type_tests_mat[ori]):
        exp_res_path = results_exp_dir + sep + test + f"_{m+1}.csv"
        df_exp = pd.read_csv(exp_res_path)
        e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
        s = df_exp["Force[kN]"]
        ax.plot(e, s, label=f"exp. {m+1}")

    ax.set_xlabel("Displacement[mm]")
    ax.set_ylabel("Force[kN]")
    ax.grid(True)

    ax.set_title(f"{material} parameter study\n Effects of the variable {var_optim} on the {test} test ")
    ax.legend()

    plt.grid(1)
    plt.show()

test_evolution("NT6_90", "UMAT", "DP780", 4)