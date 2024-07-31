import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize
from scipy.signal import savgol_filter
import os
import sklearn
import sklearn.preprocessing
from sklearn.linear_model import LinearRegression
from read_param import read_param
from get_calibration_data import export_exp_data, export_virtual_data, analyze_exp_data
from get_hardening_law import get_hardening_law


def phi1(S1, m):
    if S1.ndim==1:
        S1 = np.expand_dims(S1, axis=0)

    s11 = S1[:,0]
    s22 = S1[:,1]
    s12 = S1[:,2]

    res = np.power(np.square(s11 - s22) + 4 * np.power(s12, 4), m//2)

    return(res)

def phi2(S2, m):
    if S2.ndim==1:
        S2 = np.expand_dims(S2, axis=0)

    s11 = S2[:,0]
    s22 = S2[:,1]
    s12 = S2[:,2]

    res = np.power((3/2) * np.square(s11 + s22) + (1/2) + np.sqrt(np.square(s11 - s22) + 4 * np.power(s12, 4)), m)
    res = res + np.power((3/2) * np.square(s11 + s22) - (1/2) + np.sqrt(np.square(s11 - s22) + 4 * np.power(s12, 4)), m)

    return(res)

def yld2000(S, coeff, m):

    a1 = coeff[0]
    a2 = coeff[1]
    a3 = coeff[2]
    a4 = coeff[3]
    a5 = coeff[4]
    a6 = coeff[5]
    a7 = coeff[6]
    a8 = coeff[7]

    L1 = np.zeros((3,3))
    L2 = np.zeros((3,3))

    L1[0,0] = 2 * a1
    L1[0,1] = - a1
    L1[1,0] = - a2
    L1[1,1] = 2 * a2
    L1[2,2] = 3 * a7

    L1 = (1/3) * L1

    L2[0,0] = - 2 * a3 + 2 * a4 + 8 * a5 - 2 * a6
    L2[0,1] = - 4 * a4 + 4 * a6 + a3 - 4 * a5
    L2[1,0] = 4 * a3 - 4 * a4 - 4 * a5 + a6
    L2[1,1] = - 2 * a3 + 8 * a4 + 2 * a5 - 2 * a6
    L2[2,2] = 9 * a8

    L2 = (1/9) * L2

    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    
    S1 = np.einsum("ij,ik->ik", L1, S)
    S2 = np.einsum("ij,ik->ik", L2, S)

    p1 = phi1(S1, m)
    p2 = phi2(S2, m)

    res = 1/np.float_power(2, 1/m) * np.float_power((p1 + p2), 1/m)

    return(res)

def ys_ut_mini(thetas, coeff, m):
    """
        For given loading angles, returns the stress value for a uniaxial tensile test in the direction theta
        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function
            - powers : ndarray of shape (nmon, 3), powers of the polyN function
        Output :
            - ms : ndarray of shape (n_theta,)
    """
    n_theta = len(thetas)
    us = np.zeros((n_theta, 6))
    ms = []
    for i in range(n_theta):
        theta = thetas[i]
        us[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
    
    def f(S):
        return(yld2000(S, coeff, m))

    for u in us:
        ys = 0.1
        lamb = 1.1
        eps = 1e-7

        while f(u * ys) < 1:
            ys = ys * lamb
        s = 0.1
        e = ys
        m = (s + e) / 2
        res = f(u * m) - 1

        while (np.abs(res) > eps):
            if res > 0:
                e = m
            else:
                s = m
            m = (s + e) / 2
            res = f(u * m) - 1
        ms.append(m)
    return(np.array(ms))


def plot_check_yld2000(df, coeff, material, m):
    """
        Plot the rvalues and ys for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["YieldStress"], ["Rval"] and ["LoadAngle"] for UTs
            - coeff : ndarray of shape (nmon + 2,), coeff of the polyN mini function
            - powers : ndarray of shape (nmon, 3), powers of the polyN mini function
            - rest of the usual parameters for plot's title
    """
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi / 2, n_thetas)

    ys_model = ys_ut_mini(thetas_theo, coeff, m)
    #r_vals_model = rval_ut_mini(thetas_theo, coeff[:-2], powers)

    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    ys_exp = df_exp["YieldStress"]

    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    thetas_exp = df["LoadAngle"].iloc[index_rval].values
    r_vals_exp = df["Rval"].iloc[index_rval].values

    #plt.plot(thetas_theo, r_vals_model, c="red", label="R_val model")
    plt.plot(thetas_theo, ys_model, color="blue", label="YS model")
    plt.scatter(thetas_exp, r_vals_exp, c="red", marker="x", label="R_val exp")
    plt.scatter(thetas_exp, ys_exp/sigma0, color="blue", label="YS exp", marker="x")
    plt.title("Check poly")
    plt.xlabel(r"$\theta$[rad]", size=12)
    plt.ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.grid(1)
    plt.title(f"YLD2000 model on {material} : Yield Stresses and R-values from UT tests", size=12)

    
    plt.show()
