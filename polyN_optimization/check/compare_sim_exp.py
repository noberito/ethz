import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt

check_dir = os.path.dirname(os.path.abspath(__file__))  # Assuming current directory
sep = os.sep

polyN_cali_dir = os.path.dirname(check_dir)
sys.path.append(polyN_cali_dir)

from read import read_param
from get_calibration_data import get_tests_ori, get_ori_tests

facs_displ = {"UT": 20, "NT6" : 2., "NT20": 2., "CH": 2., "SH": 1.}
facs_thick = {"UT": 20, "NT6" : 2., "NT20": 2., "CH": 1., "SH": 8/5}
facs_width = {"UT": 20, "NT6" : 10, "NT20": 10, "CH": 20, "SH": 20}

def compare_ut_fd(material, degree, input_type, p=0, m=0):
    """
        Plot Force Displacement curves for UT except for EBT
    """

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = get_tests_ori(material)

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

    figdir = polyN_cali_dir + sep + "plots" + sep + material
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    filename = f"{material}_largestrain_{p}_{m}.png"
    filepath = figdir + sep + filename

    plt.savefig(filepath)


def compare_ut_s_2(material, degree, input_type, p=0, m=0):
    """
        Plot Stress Strain curves for UT on one and only plot
    """

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = get_tests_ori(material)

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

    figdir = polyN_cali_dir + sep + "plots" + sep + material
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    filename = f"{material}_poly{degree}_ut_s_{p}_{m}.png"
    filepath = figdir + sep + filename
    plt.savefig(filepath)


def compare_large_strain(material, degree, input_type, p=0, m=0):
    """
        Plot Force Displacement for large strain test ordered by test name
    """
    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    tests_mat = get_tests_ori(material)

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
    for type_test in tests_mat:
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
                    ax[i].plot(df_sim["U2"], df_sim["RF2"], c="red")
                    
                    if "Strain_ext" in df_sim.columns:
                        ax2 = ax[i].twinx()
                        ax2.plot(df_sim["U2"], df_sim["Strain_ext"], c="red")
                    colors = plt.cm.viridis(np.linspace(0, 0.2, type_tests_mat[ori]))
                    for k in range(1):
                        exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                        df_exp = pd.read_csv(exp_res_path)
                        e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
                        s = df_exp["Force[kN]"]
                        ax[i].plot(e, s, color=colors[k])
                        
                        indexes = {"CH": 1300, "SH":300, "NT6" : 1000, "NT20":1000}
                        ax[i].text(e[indexes[type_test]], 0.92 * s[indexes[type_test]], s="Exp.", color=colors[k])
                        if "Strain_ext" in df_sim.columns:
                            s = df_exp["AxStrain_1"]
                            ax2.plot(e,s, color=colors[k])
                            ax2.set_ylim(top= 1.5 * np.max([np.max(s), np.max(df_sim["Strain_ext"])]))
                    if "Strain_ext" in df_sim.columns:
                        ax2.set_ylabel(r"$\epsilon$ [-]")
                    ax[i].text(df_sim["U2"].iloc[20], 1.08 * df_sim["RF2"].iloc[20], s="FEA", color="red")
                    ax[i].set_title(f"{type_test}_{ori}")
                    ax[i].set_xlabel("Displacement[mm]")
                    ax[i].set_ylabel("Force[kN]")
                    ax[i].grid(True)
                    i = i + 1

    for j in range(i, nrows * ncols):
        fig.delaxes(ax[j])

    fig.suptitle(f"{material} with poly{degree} : Check Experiments vs Abaqus results\n variable {p}", fontsize=12)
    rect = np.array([0, 0.03, 1, 0.95])
    
    plt.tight_layout(rect=rect)  
    plt.subplots_adjust(hspace=0.5)

    figdir = polyN_cali_dir + sep + "plots" + sep + material
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    
    filename = f"{material}_poly{degree}_largestrain_{p}_{m}.png"
    filepath = figdir + sep + filename
    print(filepath)
    plt.savefig(filepath, dpi=600)

def compare_large_strain2(material, degree, input_type, p=0, m=0):
    """
        Plot Force Displacement for large strain test ordered by orientation
    """
    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    ori_mat = get_ori_tests(material)

    n = 0
    for ori in ori_mat:
        for type_test in ori_mat[ori]:
            if type_test != "UT":
                n += 1
    
    ncols = 3
    nrows = (n + ncols - 1) // ncols 

    fig, ax = plt.subplots(nrows, ncols, figsize=(15, 5 * nrows))
    ax = ax.flatten()  

    i = 0
    ori_mat_keys = sorted(list(ori_mat.keys()))
    print(ori_mat_keys)
    for ori in ori_mat_keys:
        for type_test in ori_mat[ori]:
            if type_test != "UT":
                sim_res_path = results_sim_dir + sep + f"{type_test}_{ori}_{input_type}_{p}_{m}.csv"
                plot = 1

                if not os.path.exists(sim_res_path):
                    plot = 0
                for k in range(ori_mat[ori][type_test]):
                    exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                    if not os.path.exists(exp_res_path):
                        plot = 0
                
                if plot:
                    df_sim = pd.read_csv(sim_res_path)
                    ax[i].plot(df_sim["U2"], df_sim["RF2"], c="red")
                    
                    if "Strain_ext" in df_sim.columns:
                        ax2 = ax[i].twinx()
                        ax2.plot(df_sim["U2"], df_sim["Strain_ext"], c="red", linestyle="dashed")
                    colors = plt.cm.viridis(np.linspace(0, 0.2, 1))
                    for k in range(1):
                        exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
                        df_exp = pd.read_csv(exp_res_path)
                        e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
                        s = df_exp["Force[kN]"]
                        ax[i].plot(e, s, color=colors[k])
                        
                        indexes = {"CH": 1300, "SH":300, "NT6" : 1000, "NT20":1000}
                        
                        facs_exp = {"CH": 0.92, "SH":1.08, "NT6" : 0.92, "NT20":0.92}
                        ax[i].text(e[indexes[type_test]], facs_exp[type_test] * s[indexes[type_test]], s="Exp.", color=colors[k])
                        if "Strain_ext" in df_sim.columns:
                            s = df_exp["AxStrain_1"]
                            ax2.plot(e,s, color=colors[k], linestyle="dashed")
                            ax2.set_ylim(top= 1.5 * np.max([np.max(s), np.max(df_sim["Strain_ext"])]))
                    if "Strain_ext" in df_sim.columns:
                        ax2.set_ylabel(r"$\epsilon$ [-]")
                    facs_sim = {"CH": 1.08, "SH": 0.92, "NT6" : 1.08, "NT20": 1.08}
                    ax[i].text(df_sim["U2"].iloc[20], facs_sim[type_test] * df_sim["RF2"].iloc[20], s="FEA", color="red")
                    ax[i].set_title(f"{type_test}_{ori}")
                    ax[i].set_xlabel("Displacement[mm]")
                    ax[i].set_ylabel("Force[kN]")
                    ax[i].grid(True)
                    i = i + 1

    for j in range(i, nrows * ncols):
        fig.delaxes(ax[j])

    fig.suptitle(f"{material} with poly{degree} : Check Experiments vs Abaqus results\n variable {p}", fontsize=12)
    rect = np.array([0, 0.03, 1, 0.95])
    
    plt.tight_layout(rect=rect)  
    plt.subplots_adjust(hspace=0.5)

    figdir = polyN_cali_dir + sep + "plots" + sep + material
    if not(os.path.exists(figdir)):
        os.makedirs(figdir)
    
    filename = f"{material}_poly{degree}_largestrain_{p}_{m}.png"
    filepath = figdir + sep + filename
    plt.savefig(filepath, dpi=600)

if __name__ == "__main__":
    p = read_param()
    material = p["material"]
    degree = int(p["degree"])
    input_type = p["input_type"]
    
    compare_large_strain2(material, degree, input_type, p=10, m=12)