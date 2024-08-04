from traceback import print_tb
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
from optimize_polyN_mini import readData_2d

polyN_dir = os.path.dirname(os.path.abspath(__file__))
sep = os.sep

def phi1(S1, m):
    if S1.ndim==1:
        S1 = np.expand_dims(S1, axis=0)

    s11 = S1[:,0]
    s22 = S1[:,1]
    s12 = S1[:,2]

    res = np.power(np.square(s11 - s22) + 4 * np.square(s12), m/2)

    return(res)

def grad_phi1(S1, m):
    if S1.ndim==1:
        S1 = np.expand_dims(S1, axis=0)

    n = len(S1)

    s11 = S1[:,0]
    s22 = S1[:,1]
    s12 = S1[:,2]

    grad_p1 = np.zeros((n, 3))

    f = np.square(s11 - s22) + 4 * np.square(s12)

    grad_f = np.zeros((n, 3))
    grad_f[:,0] = 2 * (s11 - s22)
    grad_f[:,1] = - 2 * (s11 - s22)
    grad_f[:,2] = 8 * s12

    for i in range(3):
        grad_p1[:,i] = m/2 * grad_f[:,i] * np.power(f, m/2 - 1)

    return(grad_p1)

def phi2(S2, m):
    if S2.ndim==1:
        S2 = np.expand_dims(S2, axis=0)

    s11 = S2[:,0]
    s22 = S2[:,1]
    s12 = S2[:,2]

    a = (3/2) * (s11 + s22)
    b = (1/2) * np.sqrt(np.square(s11 - s22) + 4 * np.square(s12))

    res = np.power(a + b, m) + np.power(a - b, m)
    return(res)

def grad_phi2(S2, m):
    if S2.ndim==1:
        S2 = np.expand_dims(S2, axis=0)

    n = len(S2)

    s11 = S2[:,0]
    s22 = S2[:,1]
    s12 = S2[:,2]

    a = (3/2) * (s11 + s22)
    grad_a = np.zeros((n, 3))
    grad_a[:,0] = 3/2
    grad_a[:,1] = 3/2
    grad_a[:,2] = 0

    c = np.square(s11 - s22) + 4 * np.square(s12)
    grad_c = np.zeros((n, 3))
    grad_c[:,0] = 2 * (s11 - s22)
    grad_c[:,1] = - 2 * (s11 - s22)
    grad_c[:,2] = 8 * s12

    b = (1/2) * np.sqrt(c)
    grad_b = np.zeros((n, 3))

    for i in range(3):
        grad_b[:,i] = 1/2 * 1/2 * grad_c[:,i] / np.sqrt(c) 
    
    f1 = a + b
    f2 = a - b

    grad_f1 = grad_a + grad_b
    grad_f2 = grad_a - grad_b

    grad_p2 = np.zeros((n, 3))
    
    for i in range(3):
        grad_p2[:,i] = m * grad_f1[:,i] * np.power(f1, m - 1) + m * grad_f2[:,i] * np.power(f2, m - 1)

    return(grad_p2)

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
    L2[1,1] = - 2 * a3 + 8 * a4 + 2 * a5 - 2 * a6
    L2[0,1] = a3 - 4 * a4 - 4 * a5 + 4 * a6
    L2[1,0] = 4 * a3 - 4 * a4 - 4 * a5 + a6
    L2[2,2] = 9 * a8

    L2 = (1/9) * L2

    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    
    S1 = np.einsum("ik,jk->ji", L1, S)
    S2 = np.einsum("ik,jk->ji", L2, S)

    #S1 = np.matmul(L1, S.T).T
    #S2 = np.matmul(L2, S.T).T

    p1 = phi1(S1, m)
    p2 = phi2(S2, m)

    res = np.float_power((p1 + p2)/2, 1/m)
    
    return(res)

def grad_yld2000(S, coeff, m):
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
    L2[1,1] = - 2 * a3 + 8 * a4 + 2 * a5 - 2 * a6
    L2[0,1] = a3 - 4 * a4 - 4 * a5 + 4 * a6
    L2[1,0] = 4 * a3 - 4 * a4 - 4 * a5 + a6
    L2[2,2] = 9 * a8

    L2 = (1/9) * L2

    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    
    S1 = np.einsum("ik,jk->ji", L1, S)
    S2 = np.einsum("ik,jk->ji", L2, S)

    p1 = phi1(S1, m)
    p2 = phi2(S2, m)
    grad_p1 = grad_phi1(S1, m)
    grad_p2 = grad_phi2(S2, m)

    f = p1 + p2
    f = np.expand_dims(f, axis=1)

    grad_f1 = np.einsum("ik,jk->ji", L1.T, grad_p1)
    grad_f2 = np.einsum("ik,jk->ji", L2.T, grad_p2)
    grad_f = grad_f1 + grad_f2

    grad = 1/(np.float_power(2, 1/m) * m) * grad_f * np.float_power(f, 1/m - 1)

    return(grad)

def load_coeff_yld2000(material):
    filename = "_YLD2000_pre.csv"
    filepath = polyN_dir + sep + "results_exp" + sep + material + sep + filename
    try :
        df = pd.read_csv(filepath)
        coeff = df.loc[3].values
        return coeff
    except:
        print("YLD2000 not calibrated")

def test():
    material = "DP780"
    coeff = load_coeff_yld2000(material)
    S = np.zeros(3)
    S[2] = 2
    print(grad_yld2000(S, coeff, 8))