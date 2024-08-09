import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import sklearn
import sklearn.preprocessing
import scipy.optimize as sc

file_dir = os.path.dirname(os.path.abspath(__file__))  # Assuming current directory
sep = os.sep

polyN_cali_dir = os.path.dirname(file_dir)
sys.path.append(polyN_cali_dir)

from read import read_param, get_coeff_mini, get_coeff_law
from get_calibration_data import analyze_exp_data
from optimize_polyN_mini import get_param_polyN_mini, f_min_squared


def swift(a, b, c, ep):
    return(a * np.power((b + ep), c))

def voce(a, b, c, ep):
    return(a - b * (1- np.exp(- c * ep)))

def curve_theo_norm(degree, law, ori, coeff_mini, coeff_law, ymod):
    n = 1000

    powers = get_param_polyN_mini(degree)
    def f(S):
        return np.sqrt(f_min_squared(S, coeff_mini, powers))
    
    a, b, c = coeff_law
    theta = float(ori) / 360 * 2 * np.pi

    if law == "swift":
        def hlaw(ep):
            return swift(a, b, c, ep)
    else:
        def hlaw(ep):
            return voce(a, b, c, ep)
        
    hlaw = np.vectorize(hlaw)
        
    ees = np.linspace(0, 0.003, 1000)
    S_elastic = ymod * ees
    ys = hlaw(0)
    imax = np.searchsorted(S_elastic, ys, side="right")
    ees = ees[:imax]
    S_elastic = S_elastic[:imax]

    eemax = ees[-1]

    eps = np.linspace(0, 1 - eemax, n)
    S_plastic = np.zeros(n)

    u = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
    for i in range(n):
        ep = eps[i]
        val = hlaw(ep) 
        lamb = 2
        m = 1000
        epsilon = 1e-5
        while (f(m * u) < val):
            m = m * lamb
        s = 0
        e = m
        m = (s + e) / 2
        res = f(m * u) - val
        while (abs(res) > epsilon):
            if res > 0 :
                e = m
            else :
                s = m
            m = (s + e) / 2
            res = f(m * u) - val
        S_plastic[i] = m
    
    e = np.concatenate((ees, eps + eemax))
    S = np.concatenate((S_elastic, S_plastic))
    return(e, S)

def curve_theo_ebt_norm(degree, law, coeff_mini, coeff_law, ymod):

    n = 1000

    powers = get_param_polyN_mini(degree)
    def f(S):
        return np.sqrt(f_min_squared(S, coeff_mini, powers))
    
    a, b, c = coeff_law

    if law == "swift":
        def hlaw(ep):
            return swift(a, b, c, ep)
    else:
        def hlaw(ep):
            return voce(a, b, c, ep)
        
    hlaw = np.vectorize(hlaw)
        
    ees = np.linspace(0, 0.003, 1000)
    S_elastic = ymod * ees
    ys = hlaw(0)
    imax = np.searchsorted(S_elastic, ys, side="right")
    ees = ees[:imax]
    S_elastic = S_elastic[:imax]

    eemax = ees[-1]

    eps = np.linspace(0, 1 - eemax, n)
    S_plastic = np.zeros(n)

    u = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0, 0, 0])
    for i in range(n):
        ep = eps[i]
        val = hlaw(ep) 
        lamb = 2
        m = 10
        epsilon = 1e-5
        while (f(m * u) < val):
            m = m * lamb
        s = 0
        e = m
        m = (s + e) / 2
        res = f(m * u) - val
        while (abs(res) > epsilon):
            if res > 0 :
                e = m
            else :
                s = m
            m = (s + e) / 2
            res = f(m * u) - val
        S_plastic[i] = m
    
    e = np.concatenate((ees, eps + eemax))
    S = np.concatenate((S_elastic, S_plastic))
    return(e, S)

def compare_norm(material, degree, input_type, law, coeff_mini, coeff_law, ymod, var_optim=0, n_try=0):
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
    ut_tests_mat = analyze_exp_data(material)["UT"]
    n = len(ut_tests_mat.keys())
    nrows = n // 3 + 1

    fig, ax = plt.subplots(nrows, 3)
    i = 0
    for ori in ut_tests_mat.keys():

        #PLOT SIMULATION
        j = i // 3
        k = i % 3
        sim_res_path = results_sim_dir + sep + "UT_" + ori + "_" + input_type + "_" + str(var_optim) + "_" + str(n_try) + ".csv"
        df_sim = pd.read_csv(sim_res_path)
        df_sim["S21"] = df_sim["S12"]
        df_sim["S31"] = df_sim["S13"]
        df_sim["S32"] = df_sim["S23"]
        df_sim["S"] = np.linalg.norm(df_sim[["S11", "S22", "S33", "S12", "S21", "S13", "S31", "S23", "S32"]], axis=1)
        
        ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S"],label = "Abaqus")
        ax[j,k].set_title(f"UT_{ori}")
        ax[j,k].grid(1)
        
        i = i + 1
        if ori != "EBT":
            e, Y = curve_theo_norm(degree, law, ori, coeff_mini, coeff_law, ymod)
        else:
            e, Y = curve_theo_ebt_norm(degree, law, coeff_mini, coeff_law, ymod)

        ax[j,k].plot(e, Y, label = "Analytical model")
        ax[j,k].legend()

    for i in range(n, nrows * 3):
        fig.delaxes(ax.flatten()[i])

    fig.suptitle(f"{material} with poly{degree} : Check Analytical model vs Abaqus")
    plt.legend()
    plt.show()

def curve_theo_comp(degree, law, ori, coeff_mini, coeff_law, ymod):
    n = 1000

    powers = get_param_polyN_mini(degree)
    def f(S):
        return np.sqrt(f_min_squared(S, coeff_mini, powers))
    
    a, b, c = coeff_law
    theta = float(ori) / 360 * 2 * np.pi

    if law == "swift":
        def hlaw(ep):
            return swift(a, b, c, ep)
    else:
        def hlaw(ep):
            return voce(a, b, c, ep)
        
    hlaw = np.vectorize(hlaw)
        
    ees = np.linspace(0, 0.003, 1000)
    S_elastic = ymod * ees
    ys = hlaw(0)
    imax = np.searchsorted(S_elastic, ys, side="right")
    ees = ees[:imax]
    S_elastic = S_elastic[:imax]

    eemax = ees[-1]

    eps = np.linspace(0, 1 - eemax, n)
    S_plastic = np.zeros(n)

    u = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
    for i in range(n):
        ep = eps[i]
        val = hlaw(ep) 
        lamb = 2
        m = 1000
        epsilon = 1e-5
        while (f(m * u) < val):
            m = m * lamb
        s = 0
        e = m
        m = (s + e) / 2
        res = f(m * u) - val
        while (abs(res) > epsilon):
            if res > 0 :
                e = m
            else :
                s = m
            m = (s + e) / 2
            res = f(m * u) - val
        S_plastic[i] = m
    
    e = np.concatenate((ees, eps + eemax))
    u = u[np.newaxis, :]
    S = np.concatenate((S_elastic, S_plastic))[:, np.newaxis]
    S = np.matmul(S, u)
    print(S)
    return(e, S)

def curve_theo_ebt_comp(degree, law, coeff_mini, coeff_law, ymod):
    n = 1000

    powers = get_param_polyN_mini(degree)
    def f(S):
        return np.sqrt(f_min_squared(S, coeff_mini, powers))
    
    a, b, c = coeff_law

    if law == "swift":
        def hlaw(ep):
            return swift(a, b, c, ep)
    else:
        def hlaw(ep):
            return voce(a, b, c, ep)
        
    hlaw = np.vectorize(hlaw)
        
    ees = np.linspace(0, 0.003, 1000)
    S_elastic = ymod * ees
    ys = hlaw(0)
    imax = np.searchsorted(S_elastic, ys, side="right")
    ees = ees[:imax]
    S_elastic = S_elastic[:imax]

    eemax = ees[-1]

    eps = np.linspace(0, 1 - eemax, n)
    S_plastic = np.zeros(n)

    u = np.array([1/np.sqrt(2), 1/np.sqrt(2), 0, 0, 0, 0])
    
    for i in range(n):
        ep = eps[i]
        val = hlaw(ep) 
        lamb = 2
        m = 10
        epsilon = 1e-5
        while (f(m * u) < val):
            m = m * lamb
        s = 0
        e = m
        m = (s + e) / 2
        res = f(m * u) - val
        while (abs(res) > epsilon):
            if res > 0 :
                e = m
            else :
                s = m
            m = (s + e) / 2
            res = f(m * u) - val
        S_plastic[i] = m
    
    e = np.concatenate((ees, eps + eemax))
    u = u[np.newaxis, :]
    S = np.concatenate((S_elastic, S_plastic))[:, np.newaxis]
    S = np.matmul(S, u)
    return(e, S)

def compare_comp(material, degree, input_type, law, coeff_mini, coeff_law, ymod, var_optim=0, n_try=0):
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
    ut_tests_mat = analyze_exp_data(material)["UT"]
    n = len(ut_tests_mat.keys())
    nrows = n // 3 + 1

    fig, ax = plt.subplots(nrows, 3)
    i = 0
    for ori in ut_tests_mat.keys():

        #PLOT SIMULATION
        j = i // 3
        k = i % 3
        sim_res_path = results_sim_dir + sep + "UT_" + ori + "_" + input_type + "_" + str(var_optim) + "_" + str(n_try) + ".csv"
        df_sim = pd.read_csv(sim_res_path)

        s11 = df_sim["S11"]
        s22 = df_sim["S22"]
        s12 = df_sim["S12"]
        
        ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S11"] , label = "Abaqus sxx")
        ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S22"] , label = "Abaqus syy")
        ax[j,k].plot(df_sim["SDV_EPBAR"], df_sim["S12"] , label = "Abaqus sxy")
        ax[j,k].set_title(f"UT_{ori}")
        ax[j,k].grid(1)
        
        i = i + 1
        if ori != "EBT":
            e, Y = curve_theo_comp(degree, law, ori, coeff_mini, coeff_law, ymod)
        else:
            e, Y = curve_theo_ebt_comp(degree, law, coeff_mini, coeff_law, ymod)

        ax[j,k].plot(e, Y[:,0], label = "Analytical model sxx")
        ax[j,k].plot(e, Y[:,1], label = "Analytical model syy")
        ax[j,k].plot(e, Y[:,3], label = "Analytical model sxy")
        ax[j,k].legend()

    for i in range(n, nrows * 3):
        fig.delaxes(ax.flatten()[i])

    fig.suptitle(f"{material} with poly{degree} : Check Analytical model vs Abaqus")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    p = read_param()
    material = p["material"]
    degree = int(p["degree"])
    law = p["law"]
    input_type = p["input_type"]

    coeff_mini = get_coeff_mini(material, degree)
    coeff_law, ymod = get_coeff_law(material, law)
    
    compare_norm(material, degree, input_type, law, coeff_mini, coeff_law, ymod)
    compare_comp(material, degree, input_type, law, coeff_mini, coeff_law, ymod)