import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

current_dir = "./"  # Assuming current directory
dir = "/"
if os.name == "nt":
    current_dir = ".\\"
    dir = "\\"

material = "DP780"

postprocess_dir = os.path.dirname(os.path.abspath(__file__))
polyN_cali_dir = os.path.dirname(postprocess_dir)



sim_params = [["UT_00", 0], ["UT_15", 1], ["UT_30",2], ["UT_45",3], ["UT_60",4], ["UT_75",5],["UT_90",6], ["UT_EBT",7]]
""""CH_00", "CH_45",
"NT6_00", "NT6_45", "NT6_90",
"NT20_00", "NT20_45",
"SH_000", "SH_090", "SH_p45"]"""

tests = ["UT_00", "UT_15", "UT_30", "UT_45", "UT_60", "UT_75","UT_90", "UT_EBT"]
""""CH_00", "CH_45",
"NT6_00", "NT6_45", "NT6_90",
"NT20_00", "NT20_45",
"SH_000", "SH_090", "SH_p45"]"""

def compare(tests, material):

    results_exp_dir = polyN_cali_dir + dir + "results_exp" + dir + material
    results_sim_dir = polyN_cali_dir + dir + "results_sim" + dir + material

    n = len(tests)

    fig, ax = plt.subplots(n // 3 + 1, 3)

    for i in range(n):
        test = tests[i]
        j = i // 3
        k = i % 3

        exp_res_path = results_exp_dir + dir + "DATA" + dir + test + "_1.csv"
        sim_res_path = results_sim_dir + dir + test + "_polyN.csv"

        df_sim = pd.read_csv(sim_res_path)
        df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S13", "S23"]], axis=1)
        df_sim["E"] = np.linalg.norm(df_sim[["E11", "E22", "E33", "E12", "E13", "E23"]], axis=1)
        df_exp = pd.read_csv(exp_res_path)
        if test == "UT_EBT":
            df_exp = df_exp.rename(columns={df_exp.columns[0]:"TrueStrain", df_exp.columns[1]:"TrueStress[MPa]"})
        ax[j,k].scatter(df_sim["E"], df_sim["S"], marker = "x",label = "sim")
        ax[j,k].scatter(df_exp["TrueStrain"], df_exp["TrueStress[MPa]"], marker = "x", label = "exp")   
        ax[j,k].set_title(test)  

    plt.legend()
    plt.tight_layout()
    plt.show()

compare(tests, material)