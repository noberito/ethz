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

from read_param import read_param
from get_calibration_data import analyze_exp_data

p = read_param()
material = p["material"]
degree = int(p["degree"])
law = p["law"]

polyN_coeff_path = polyN_cali_dir + sep + "polyN_coeff.npy"
polyN_coeff = np.load(polyN_coeff_path)
law_coeff_path = polyN_cali_dir + sep + f"{law}_coeff.npy"
law_coeff = np.load(law_coeff_path)

def get_param_polyN(degree):
    """
        Returns the parameters of polyN according to the degree
        Input :
            - degree : integer, degree of the polyN function
        Output :
            - nmon : integer, number of monomials im the polyN function of degree n
            - nmon_abq : integer, number of monomials of degree n in a homogeneous function of degree n
            - powers : ndarray of shape (nmon, 5), powers[i, j] is the power of the variable j in the monomial i
    """
    x0 = np.zeros((2,5))
    polyN = sklearn.preprocessing.PolynomialFeatures((degree, degree), include_bias=False)
    X = polyN.fit_transform(x0)
    powers = polyN.powers_
    nmon_degree = len(powers)

    selected_indices = []
    for i in range(nmon_degree):
        k,l,m,n,p = powers[i]
        if ((m==n) and (n==p)) or ((m%2==0) and (n%2==0) and (p%2==0)):
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
    powers = powers.T[sorted_indices]
    X = X[:,sorted_indices]

    ndata, nmon = X.shape

    return(powers, nmon)

powers, nmon = get_param_polyN(degree)

def dev(S):
    """
        Returns the deviatoric stress.
        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
        Output :
            - D : ndarray of shape (n,6), deviatoric stress components in 3D
    
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = S.copy()
    trace_dev = S[:,0] + S[:,1] + S[:,2]
    D[:,0] = S[:,0] - (1/3) * trace_dev
    D[:,1] = S[:,1] - (1/3) * trace_dev
    D[:,2] = S[:,2] - (1/3) * trace_dev
    return(D)

def polyN(S, coeff):
    """
        Compute the polyN function.
        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff : ndarray of shape (nmon,), coefficients of the polyN function
        Output :
            - res : float, result of polyN
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = dev(S)
    X = D[:,[0, 1, 3, 4, 5]]
    res = np.zeros(len(X))
    for i in range(nmon) :
        p = powers[i]
        res = res + coeff[i] * np.prod(X ** p, axis=1)
    return(res)

def swift(a, b, c, ep):
    return(a * np.power((b + ep), c))

def voce(a, b, c, ep):
    return(a - b * (1- np.exp(- c * ep)))

def curve_theo(ori, polyN_coeff, law_coeff):
    n = 1000
    def f(S):
        return polyN(S, polyN_coeff)**(1/degree)
    
    a, b, c, ymod = law_coeff
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
        epsilon = 1e-4
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

def curve_theo_ebt(polyN_coeff, law_coeff):
    n = 1000
    def f(S):
        return polyN(S, polyN_coeff)**(1/degree)
    
    a, b, c, ymod = law_coeff

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

    u = np.array([1, 1, 0, 0, 0, 0])
    for i in range(n):
        ep = eps[i]
        val = hlaw(ep) 
        lamb = 2
        m = 10
        epsilon = 1e-4
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

def compare(material, polyN_coeff, law_coeff):
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
        sim_res_path = results_sim_dir + sep + "UT_" + ori + "_polyN.csv"
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
            e, Y = curve_theo(ori, polyN_coeff, law_coeff)
        else:
            e, Y = curve_theo_ebt(polyN_coeff, law_coeff)

        ax[j,k].plot(e, Y, label = "Analytical model")
        ax[j,k].legend()

    for i in range(n, nrows * 3):
        fig.delaxes(ax.flatten()[i])

    fig.suptitle(f"{material} with poly{degree} : Check Analytical model vs Abaqus")
    plt.legend()
    plt.show()

compare(material, polyN_coeff, law_coeff)