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

facs_displ = {"UT": 20, "NT6" : 2., "NT20": 2., "CH": 2., "SH": 1.}
facs_thick = {"UT": 20, "NT6" : 2., "NT20": 2., "CH": 1., "SH": 8/5}
facs_width = {"UT": 20, "NT6" : 10, "NT20": 10, "CH": 20, "SH": 20}

def compare_ut_fd(material, degree, input_type, p=0, m=0):
    """
        Plot Force Displacement curves for UT except for EBT
    """

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = analyze_exp_data(material)

    n = 0
    for type_test in tests_mat:
        if type_test == "UT":
            for ori in tests_mat[type_test]:
                if ori != "EBT":
                    n += 1
    
    ncols = 3
    nrows = (n + ncols - 1) // ncols 

    fig, ax = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows))
    ax = ax.flatten()  

    i = 0
    for type_test in tests_mat.keys():
        type_tests_mat = tests_mat[type_test]
        if type_test == "UT":

            for ori in type_tests_mat.keys():
                if ori != "EBT":
                    sim_res_path = results_sim_dir + sep + f"{type_test}_{ori}_{input_type}_{p}_{m}.csv"
                    print(sim_res_path)
                    plot = 1

                    if not os.path.exists(sim_res_path):
                        plot = 0
                    for k in range(type_tests_mat[ori]):
                        exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                        if not os.path.exists(exp_res_path):
                            plot = 0
                    
                    if plot:
                        df_sim = pd.read_csv(sim_res_path)
                        ax[i].plot(df_sim["U2"], df_sim["RF2"], label="abaqus", c="red")

                        for k in range(type_tests_mat[ori]):
                            exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                            df_exp = pd.read_csv(exp_res_path)
                            e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
                            s = df_exp["Force[kN]"]
                            ax[i].plot(e, s, label=f"exp. {k+1}")

                        ax[i].set_title(f"{type_test}_{ori}")
                        ax[i].set_xlabel("Displacement[mm]")
                        ax[i].set_ylabel("Force[kN]")
                        ax[i].grid(True)
                        ax[i].legend()
                        i = i + 1

    for j in range(i, nrows * ncols):
        fig.delaxes(ax[j])

    fig.suptitle(f"{material} with poly{degree} : Check Experiments vs Abaqus results\n variable {p}", fontsize=12)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  
    plt.subplots_adjust(hspace=0.5)

    figdir = file_dir + sep + material + sep + "opti"
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    
    filename = f"{material}_largestrain_{p}_{m}.png"
    filepath = figdir + sep + filename

    plt.savefig(filepath)

def compare_ut_s(material, func, degree, input_type, p=0, m=0):
    """
        Plot Stress Strain curves for each UT except for EBT
    """

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = analyze_exp_data(material)

    n = 0

    for ori in tests_mat["UT"]:
        n = n + 1
    nrows = n // 3 + 1
    fig, ax = plt.subplots(nrows, 3)
    i = 0
    for type_test in tests_mat.keys():
        type_tests_mat = tests_mat[type_test]
        if type_test == "UT":
            for ori in type_tests_mat.keys():
                if ori != "EBT":
                    j = i // 3
                    k = i % 3
                    sim_res_path = results_sim_dir + sep + f"UT_{ori}_{input_type}_{p}_{m}.csv"
                    df_sim = pd.read_csv(sim_res_path)
                    df_sim["S21"] = df_sim["S12"]
                    df_sim["S31"] = df_sim["S13"]
                    df_sim["S32"] = df_sim["S23"]
                    df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S21", "S13", "S31", "S23", "S32"]], axis=1)

                    ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S"],label = "abaqus")

                    for k in range(type_tests_mat[ori]):
                        exp_res_path = results_exp_dir + sep + "UT_" + ori + f"_{k+1}.csv"
                        df_exp = pd.read_csv(exp_res_path)
                        if ori == "EBT":
                            df_exp = df_exp.rename(columns={df_exp.columns[0]:"TrueStrain", df_exp.columns[1]:"TrueStress[MPa]"})
                            e = df_exp["TrueStrain"]
                            s = df_exp["TrueStress[MPa]"] * np.sqrt(2) #In the EBT file, True stress corresponds to the stress applied in one direction which is 1/sqrt(2) * S
                            ax[j,k].plot(e, s, label = f"exp. {k+1}")
                        else:
                            imax = df_exp["EngStress[MPa]"].idxmax()
                            e = df_exp["TrueStrain"].values[:imax]
                            s = df_exp["TrueStress[MPa]"].values[:imax]
                            ax[j,k].plot(e, s, label = f"exp. {k+1}") 

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

    figdir = file_dir + sep + material + sep + "var_" + str(p)
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    
    filename = f"{material}_ut_{func}_{p}_{m}.png"
    filepath = figdir + sep + filename
    print(filepath)
    plt.savefig(filepath)

def compare_ut_s_2(material, func, degree, input_type, p=0, m=0):
    """
        Plot Stress Strain curves for UT on one and only plot
    """

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = analyze_exp_data(material)

    n = 0

    colors = plt.cm.viridis(np.linspace(0, 1, 9))
    fig, ax = plt.subplots(1)
    i = 0
    for type_test in tests_mat.keys():
        type_tests_mat = tests_mat[type_test]
        if type_test == "UT":
            for ori in type_tests_mat.keys():
                sim_res_path = results_sim_dir + sep + f"UT_{ori}_{input_type}_{p}_{m}.csv"
                df_sim = pd.read_csv(sim_res_path)
                df_sim["S21"] = df_sim["S12"]
                df_sim["S31"] = df_sim["S13"]
                df_sim["S32"] = df_sim["S23"]
                df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S21", "S13", "S31", "S23", "S32"]], axis=1)


                

                for k in range(min(type_tests_mat[ori], 1)):
                    exp_res_path = results_exp_dir + sep + "UT_" + ori + f"_{k+1}.csv"
                    df_exp = pd.read_csv(exp_res_path)
                    if ori == "EBT":
                        df_exp = df_exp.rename(columns={df_exp.columns[0]:"PlasticStrain", df_exp.columns[1]:"PlasticStress[MPa]"})
                        e = df_exp["PlasticStrain"]
                        s = df_exp["PlasticStress[MPa]"] * np.sqrt(2) + 2000 #In the EBT file, True stress corresponds to the stress applied in one direction which is 1/sqrt(2) * S
                        ax.plot(e, s, color=colors[i])
                        ax.plot(df_sim["SDV_EPBAR"], df_sim["S"] + 2000,linewidth=1,linestyle="--", label = "UT_" + ori, color=colors[i])
                    else:
                        imax = df_exp["PlasticStrain_longi"].idxmax()
                        e = df_exp["PlasticStrain_longi"].values[:imax]
                        s = df_exp["PlasticStress[MPa]"].values[:imax] + (float(ori))/90 * 1500
                        ax.plot(e, s, color=colors[i]) 
                        ax.plot(df_sim["SDV_EPBAR"], df_sim["S"] + (float(ori))/90 * 1500,linewidth=1, linestyle="--", label = "UT_" + ori, color=colors[i])

                ax.set_xlabel(r"$\epsilon^{p}$")
                ax.set_ylabel(r"$\sigma[MPa]$")
                ax.grid(1)
                i = i + 1

    fig.suptitle(f"{material} with poly{degree} : Check UTs Experiments vs Abaqus", fontsize=12)
    ax.set_title("-- : abaqus")
    plt.legend()

    figdir = file_dir + sep + material + sep + "var_" + str(p)
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    
    filename = f"{material}_ut_{func}_{p}_{m}.png"
    filepath = figdir + sep + filename
    plt.savefig(filepath)

def compare_all(material, degree, input_type):
    """
        Plot Stress Strain curves for UT and force displacement for large strain test
    """

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = analyze_exp_data(material)

    n = 0

    for type_test in tests_mat:
        for ori in tests_mat[type_test]:
            n = n + 1
    nrows = n // 3 + 1
    fig, ax = plt.subplots(nrows, 3)
    i = 0
    for type_test in tests_mat.keys():
        type_tests_mat = tests_mat[type_test]
        if type_test == "UT":
            for ori in type_tests_mat.keys():
                j = i // 3
                k = i % 3
                sim_res_path = results_sim_dir + sep + "UT_" + ori + "_polyN.csv"
                df_sim = pd.read_csv(sim_res_path)
                df_sim["S21"] = df_sim["S12"]
                df_sim["S31"] = df_sim["S13"]
                df_sim["S32"] = df_sim["S23"]
                df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S21", "S13", "S31", "S23", "S32"]], axis=1)

                ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S"],label = "abaqus")

                for k in range(type_tests_mat[ori]):
                    exp_res_path = results_exp_dir + sep + "UT_" + ori + f"_{k+1}.csv"
                    df_exp = pd.read_csv(exp_res_path)
                    if ori == "EBT":
                        df_exp = df_exp.rename(columns={df_exp.columns[0]:"TrueStrain", df_exp.columns[1]:"TrueStress[MPa]"})
                        e = df_exp["TrueStrain"]
                        s = df_exp["TrueStress[MPa]"] * np.sqrt(2) #In the EBT file, True stress corresponds to the stress applied in one direction which is 1/sqrt(2) * S
                        ax[j,k].plot(e, s, label = f"exp. {k+1}")
                    else:
                        imax = df_exp["EngStress[MPa]"].idxmax()
                        e = df_exp["TrueStrain"].values[:imax]
                        s = df_exp["TrueStress[MPa]"].values[:imax]
                        ax[j,k].plot(e, s, label = f"exp. {k+1}") 

                ax[j,k].set_title(f"UT_{ori}")
                ax[j,k].set_xlabel(r"$\epsilon$")
                ax[j,k].set_ylabel(r"$\sigma$")
                ax[j,k].grid(1)
                ax[j,k].legend()
                i = i + 1
        else :
            for ori in type_tests_mat.keys():
                j = i // 3
                k = i % 3
                sim_res_path = results_sim_dir + sep + type_test + "_" + ori + "_" + input_type + ".csv"

                if os.path.exists(sim_res_path):
                    df_sim = pd.read_csv(sim_res_path)
                    ax[j,k].plot(df_sim["U2"], df_sim["RF2"],label = "abaqus")

                for k in range(type_tests_mat[ori]):
                    exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                    if os.path.exists(exp_res_path):
                        df_exp = pd.read_csv(exp_res_path)
                        if type_test=="SH":
                            e = df_exp["Displacement longi[mm]"].values
                        else : 
                            e = df_exp["Displacement[mm]"].values
                        s = df_exp["Force[kN]"]
                        ax[j,k].plot(e, s, label = f"exp. {k+1}") 

                ax[j,k].set_title(f"{type_test}_{ori}")
                ax[j,k].set_xlabel("Displacement[mm]")
                ax[j,k].set_ylabel("Force[kN]")
                ax[j,k].grid(1)
                ax[j,k].legend()
                i = i + 1

    for i in range(n, nrows * 3):
        fig.delaxes(ax.flatten()[i])

    fig.suptitle(f"{material} with poly{degree} : Check Experiments vs Abaqus results", fontsize=12)
    plt.legend()
    plt.show()

import os
import pandas as pd
import matplotlib.pyplot as plt

def compare_large_strain(material, func, degree, input_type, p=0, m=0):
    """
        Plot Force Displacement for large strain test
    """
    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = analyze_exp_data(material)

    n = 0
    for type_test in tests_mat:
        if type_test != "UT":
            for ori in tests_mat[type_test]:
                n += 1
    
    ncols = 3
    nrows = (n + ncols - 1) // ncols 

    fig, ax = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows))
    ax = ax.flatten()  

    i = 0
    for type_test in tests_mat.keys():
        type_tests_mat = tests_mat[type_test]
        if type_test != "UT":

            for ori in type_tests_mat.keys():
                sim_res_path = results_sim_dir + sep + f"{type_test}_{ori}_{input_type}_{p}_{m}.csv"
                plot = 1

                if not os.path.exists(sim_res_path):
                    plot = 0
                for k in range(type_tests_mat[ori]):
                    exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                    if not os.path.exists(exp_res_path):
                        plot = 0
                
                if plot:
                    df_sim = pd.read_csv(sim_res_path)
                    ax[i].plot(df_sim["U2"], df_sim["RF2"], label="abaqus", c="red")
                    if "Strain_ext" in df_sim.columns:
                        ax2 = ax[i].twinx()
                        ax2.plot(df_sim["U2"], df_sim["Strain_ext"], label="abaqus", c="red")
                    colors = plt.cm.viridis(np.linspace(0, 0.2, type_tests_mat[ori]))
                    for k in range(type_tests_mat[ori]):
                        exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                        df_exp = pd.read_csv(exp_res_path)
                        e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
                        s = df_exp["Force[kN]"]
                        ax[i].plot(e, s, color=colors[k])
                        if "Strain_ext" in df_sim.columns:
                            s = df_exp["AxStrain_1"]
                            ax2.plot(e,s, color=colors[k])
                            ax2.set_ylim(top= 1.5 * np.max([np.max(s), np.max(df_sim["Strain_ext"])]))

                    ax[i].set_title(f"{type_test}_{ori}")
                    ax[i].set_xlabel("Displacement[mm]")
                    ax[i].set_ylabel("Force[kN]")
                    ax[i].grid(True)
                    ax[i].legend()
                    i = i + 1

    for j in range(i, nrows * ncols):
        fig.delaxes(ax[j])

    fig.suptitle(f"{material} with poly{degree} : Check Experiments vs Abaqus results\n variable {p}", fontsize=12)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  
    plt.subplots_adjust(hspace=0.5)

    figdir = file_dir + sep + material + sep + "var_" + str(p)
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    
    filename = f"{material}_largestrain_{func}_{p}_{m}.png"
    filepath = figdir + sep + filename
    print(filepath)
    plt.savefig(filepath)

if __name__ == "__main__":
    p = read_param()
    material = p["material"]
    degree = int(p["degree"])
    input_type = p["input_type"]
    func = p["func"]
    compare_ut_s_2(material, func, degree, input_type, p=10, m=10)