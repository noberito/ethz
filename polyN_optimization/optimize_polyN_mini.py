
#!/usr/bin/env python
# coding: utf-8

# In[200]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import os
import time
import tensorflow as tf
from tensorflow import keras
from keras import Sequential
from tensorflow.keras.layers import Dense, Input
from keras.models import load_model
import multiprocessing
from pickle import load
import sklearn
import sklearn.preprocessing 
from read_param import read_param
from get_calibration_data import export_exp_data, export_virtual_data
import shutil
from get_hardening_law import get_hardening_law
from convexity_domain.generating_convexity_data import check_convexity_lpm
import random


# In[201]:


file_dir = os.path.dirname(os.path.abspath(__file__))  
sep = os.sep
exec_dir = file_dir + sep + "running"


# In[202]:


"---------------------------------------------------------------- PARAMETERS ---------------------------------------------------------------------------------------------"

p = read_param()

material = p["material"]
gseed = int(p["gseed"])
enu = float(p["enu"])
density = float(p["density"])
nb_virtual_pt = int(p["nb_virtual_pt"])
degree = int(p["degree"])
weight_exp = float(p["weight_exp"])
weight_rval = float(p["weight_rval"])
protomodel = p["protomodel"]
law = p["law"]

gen_v_data = int(p["gen_v_data"])
gen_e_data = int(p["gen_e_data"])

opti = int(p["opti"])
loadcoeff = int(p["loadcoeff"])
adapt = int(p["adapt"])

export_coeff_abq = int(p["export_coeff_abq"])
export_coeff_user = int(p["export_coeff_user"])


# In[203]:


"""----------------------------------------------------------------- GENERATING DATA ----------------------------------------------------------------------------------"""
thetas = ["00", "15", "30", "45", "60", "75", "90"]
materials = ["AA7020-T6", "DP600", "DP780"]

def mises(sigma):
    """
    Vectorized Mises function
    """
    sigma = np.atleast_2d(sigma)
    if sigma.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    s11 = sigma[:, 0]
    s22 = sigma[:, 1]
    s33 = sigma[:, 2]
    s12 = sigma[:, 3]
    s13 = sigma[:, 4]
    s23 = sigma[:, 5]

    res = np.sqrt(0.5 * (s22 - s33)**2 + 0.5 * (s33 - s11)**2 + 0.5 * (s11 - s22)**2 + 3 * (s23**2 + s13**2 + s12**2))
    return res

# In[204]:

""" --------------------------------------------------------COPY LAB DATA----------------------------------------------"""

def copy_lab_data(material):
    copy_dir = f"{file_dir}{sep}results_exp{sep}{material}"
    if not os.path.exists(copy_dir):
        parent_dir = os.path.dirname(file_dir)
        lab_data_dir = f"{parent_dir}{sep}lab_data"
        mat_lab_data_dir = f"{lab_data_dir}{sep}{material}_results{sep}DATA"
        shutil.copytree(mat_lab_data_dir, copy_dir)

""" ----------------------------------------------------------- READ DATA ----------------------------------------------------------------------------------------------"""
def readData_2d(material, protomodel):
    """
        Input :
            - material : string
            - protomodel : string
    """

    folderpath = f"{file_dir}{sep}calibration_data{sep}{material}"
    filename_e = "data_exp_" + material + ".csv"
    filename_v = "data_virtual_" + material + "_" + protomodel + ".csv"

    filepath_e = folderpath + sep + filename_e
    filepath_v = folderpath + sep + filename_v

    df_e = pd.read_csv(filepath_e)
    df_v = pd.read_csv(filepath_v)

    df_e["LoadAngle"] = 2 * np.pi / 360 * df_e["LoadAngle"]

    sigma0 = df_e["YieldStress"].iloc[0]
    
    df_e["s11"] = df_e["YieldStress"] / sigma0 * (np.square(np.cos(df_e["LoadAngle"])) + df_e["q"] * np.square(np.sin(df_e["LoadAngle"])))
    df_e["s22"] = df_e["YieldStress"] / sigma0 * (np.square(np.sin(df_e["LoadAngle"])) + df_e["q"] * np.square(np.cos(df_e["LoadAngle"])))
    df_e["s33"] = df_e["YieldStress"] * 0
    df_e["s12"] = df_e["YieldStress"] / sigma0 * (1 - df_e["q"]) * np.sin(df_e["LoadAngle"]) * np.cos(df_e["LoadAngle"])
    df_e["s13"] = df_e["YieldStress"] * 0
    df_e["s23"] = df_e["YieldStress"] * 0
    
    df = pd.concat([df_e, df_v])
    
    df["Norm"] = np.linalg.norm(df[["s11", "s22", "s33", "s12", "s13", "s23"]].values, axis=1)
    df["d11"] = df["s11"] - df["s33"]
    df["d22"] = df["s22"] - df["s33"]

    return(df)




# In[205]:


"""------------------------ MODEL AND SCALER LOADING -------------------------------------------------"""




"""------------------------------------POLYN PARAMETERS---------------------------------------------------------"""

def get_param_polyN_mini(degree):
    """
        Returns the parameters of polyN according to the degree
        Input :
            - degree : integer, degree of the polyN function
        Output :
            - powers : ndarray of shape (nmon, 5), powers[i, j] is the power of the variable j in the monomial i
    """
    x0 = np.zeros((2,3))
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X = polyN.fit_transform(x0)
    powers = polyN.powers_
    nmon_degree = len(powers)

    selected_indices = []
    for i in range(nmon_degree):
        k,l,m = powers[i]
        if (k + l + m == degree) and (m%2==0) :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2]))
    powers = powers.T[sorted_indices]
    X = X[:,sorted_indices]

    return(powers)

powers = get_param_polyN_mini(degree)
print(powers)


#degree = int(p["degree"])
#powers, nmon = get_param_polyN(degree)
# In[207]:


""" ---------------------------------------------------POLYN DEFINITION------------------------------------------------------------------------------"""
def dev_min(S):
    """
        Returns the deviatoric stress.
        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
        Output :
            - D : ndarray of shape (n,6), deviatoric stress components in 3D
    
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = np.zeros((len(S), 5))
    D[:,0] = S[:,0] - S[:,2]
    D[:,1] = S[:,1] - S[:,2]
    D[:,2] = S[:,3]
    D[:,3] = S[:,4]
    D[:,4] = S[:,5]
    return(D)

def polyN_2d(S, coeff_min, powers):
    """
        Compute the polyN function.
        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff : ndarray of shape (nmon,), coefficients of the polyN function
        Output :
            - res : float, result of polyN
    """
    degree = np.sum(powers[0])
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = dev_min(S)[:,[0,1,2]]
    res = np.zeros(len(X))
    nmon = len(powers)
    for i in range(nmon) :
        p = powers[i]
        res = res + coeff_min[i] * np.prod(X ** p, axis=1)
    return(res)

"""---------------------------------------------- GRADIENT AND HESSIAN OF POLYN DEFINITION -----------------------------------------------------------"""

#True we lose time here but ok

def jac_dev_min():
    """
        Returns the jacobian of the deviatoric minimalistic operator
        Output :
            - jac : ndarray of shape (5,6)
    """
    jac = np.zeros((3, 6))
    jac[0] = np.array([1, 0, -1, 0, 0, 0])
    jac[1] = np.array([0, 1, -1, 0, 0, 0])
    jac[2] = np.array([0, 0, 0, 1, 0, 0])
   
    return(jac)

def jac_polyN_2d_param(coeff_min, powers):
    """
        Compute the different parameters and coefficients to compute the Jacobian of polyN
        Input :
            - coeff (float ndarray of shape (nmon)) : Coefficients of the polyN function
            ordered by power[1] asc, power[2] asc, power[3] asc, power[4] asc
            - powers (float ndarray of shape (nmon, 5)) : Powers for each monomial of the 
            PolyN function
        
        Output :
            - coeff_grad (float ndarray of shape (5, nmon)) : coeff_grad[i,j] contains the 
            coefficients of the monomial j derived with respect to dev[i]
            - powers_grad (float ndarray of shape (5, nmon, 5)) : powers_grad[i,j] contains
            the monomial j derived with respect to dev[i]
    """

    coeff_grad = np.zeros((3, len(coeff_min)))
    powers_grad = np.zeros((3, powers.shape[0], powers.shape[1]))
    for i in range(3):
        coeff_grad[i] = coeff_min * powers[:,i]
        subs = np.zeros(powers_grad[i].shape)
        subs[:,i] = -1
        powers_grad[i] = np.maximum(0,powers + subs)
    
    return(coeff_grad, powers_grad)

def grad_polyN_2d(S, coeff_grad, powers_grad):
    """
        Input :
            - S : float ndarray of shape (n, 6) or (6,), stress components in 3d
            - coeff_grad : float ndarray of shape (5, nmon),coeff_grad[i] contains the 
            coefficients of each monomial derived with respect to dev[i]
            - powers_grad : float ndarray of shape (5, nmon, 5) : powers_grad[i][j] contains
            the monomial j derived with respect to i-th deviatoric component
        
        Output :
            - grad : ndarray of shape (n, 6), grad_polyN[i, j] = dpolyN/dsj of the i-th data
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = dev_min(S)[:,[0,1,2]]

    grad_polyN = np.zeros(6)
    grad_f = np.zeros((len(X),3))

    nmon = len(powers_grad[0])

    for i in range(3):
        for j in range(nmon):
            p = powers_grad[i][j]
            grad_f[:,i] = grad_f[:,i] + coeff_grad[i][j] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev_min()

    grad_polyN = np.dot(grad_f, jac_d)
    return(grad_polyN)

def hessian_polyN_2d_param(coeff_grad, powers_grad):
    """
    Compute the different parameters and coefficients to compute the Hessian of polyN
        Input :
            - coeff_grad : float ndarray of shape (5, nmon), coeff_grad[i] contains the 
            coefficients of dpolyN/ddev[i]
            - powers_grad : float ndarray of shape (5, nmon, 5), powers_grad[i,j] contains
            the powers of dmon[j]/ddev[i]
        
        Output :
            - coeff_hessian : float ndarray of shape (5, 5, nmon)), coeff_hessian[i,j,k] contains
            the coefficients of d2mon[k]/ddev[j]ddev[i] 
            - powers_hessian : float ndarray of shape (5, 5, nmon, 5), powers_hessian[i,j,k]contains
            the powers of d2mon[k]/ddev[j]ddev[i] 

    """
    coeff_hessian = np.zeros((3, coeff_grad.shape[0], coeff_grad.shape[1]))
    powers_hessian = np.zeros((3, powers_grad.shape[0], powers_grad.shape[1], powers_grad.shape[2]))

    for i in range(3):
        for j in range(3):
            coeff_hessian[i][j] = coeff_grad[i] * powers_grad[i,:,j]
            subs = np.zeros(powers_hessian[i][j].shape)
            subs[:,j] = -1
            powers_hessian[i][j] = np.maximum(0,powers_grad[i] + subs)
    
    return(coeff_hessian, powers_hessian)

def hessian_polyN_2d(S, coeff_hessian, powers_hessian):
    """
        Compute the hessian of polyN.
        Input :
            - S : float ndarray of shape (n, 6) or (6,), stress components in 3d
            - coeff_hessian : float ndarray of shape (5, 5, nmon)), coeff_hessian[i,j,k] contains
            the coefficients of d2mon[k]/ddev[j]ddev[i] 
            - powers_hessian : float ndarray of shape (5, 5, nmon, 5), powers_hessian[i,j,k]contains
            the powers of d2mon[k]/ddev[j]ddev[i] 
        Output :
            - hessian : float ndarray of shape (n, 6, 6), hessian[i,j,k] = dpolyN/ds[k]ds[j] of the ith data pt

    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = dev_min(S)[:,[0,1,2]]

    hessian_polyN = np.zeros((6,6))
    jac_grad_f = np.zeros((len(X), 3, 3))
    nmon = len(powers_hessian)

    for i in range(3):
        for j in range(3):
            for k in range(nmon):
                p = powers_hessian[i][j][k]
                jac_grad_f[:,i,j] = jac_grad_f[:,i,j] + coeff_hessian[i][j][k] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev_min()
    hessian_polyN = np.transpose(np.dot(jac_d.T,np.dot(jac_grad_f[:], jac_d)), (1, 2, 0))
    return(hessian_polyN)

def f1(S, coeff_mini, powers):
    degree = np.sum(powers[0])
    res = np.float_power(polyN_2d(S, coeff_mini, powers),(2/degree))
    return(res)

    
def f2(S, coeff):
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = dev_min(S)
    res = coeff[-2] * np.square(X[:,3]) + coeff[-1] * np.square(X[:,4])
    return(res)


# In[208]:

def f_min_squared(S, coeff, powers):
    return(f1(S, coeff[:-2], powers) + f2(S, coeff[-2:]))

""" ---------------------------------------------PARAMETERS OPTIMIZATION-----------------------------------------------------------------------------"""


def rval_ut(thetas, coeff_mini, powers):
    n_theta = len(thetas)
    S = np.zeros((n_theta, 6))
    v2 = np.zeros((n_theta, 3))
    r_val = np.zeros(n_theta)
    for i in range(n_theta):
        theta = thetas[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        v2[i] = np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])

    degree = np.sum(powers[0])
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    grad_P = grad_polyN_2d(S, coeff_grad, powers_grad)[:,[0,1,3]]
    grad_f_plane = np.zeros((n_theta,3))
    for i in range(n_theta):
        for j in range(3):
            grad_f_plane[i,j] = (2 / degree) * grad_P[i,j] * np.float_power(polyN_2d(S[i], coeff_mini, powers),(2/degree - 1))[0]
    for i in range(n_theta):
        r_val[i] = - np.dot(grad_f_plane[i], v2[i]) / (grad_f_plane[i,0] + grad_f_plane[i,1])
    return(r_val)

def eq_rval_ut(thetas, r_val_theo, coeff_mini, powers):
    n_theta = len(thetas)
    S = np.zeros((n_theta, 6))
    v2 = np.zeros((n_theta, 3))
    for i in range(n_theta):
        theta = thetas[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        v2[i] = np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])

    degree = np.sum(powers[0])
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    grad_P = grad_polyN_2d(S, coeff_grad, powers_grad)[:, [0,1,3]]
    grad_f_plane = np.zeros((n_theta, 3))
    for i in range(n_theta):
        for j in range(3):
            grad_f_plane[i,j] = (2 / degree) * grad_P[i,j] * np.float_power(polyN_2d(S[i], coeff_mini, powers),(2/degree - 1))[0]
    eq = np.zeros(n_theta)
    for i in range(n_theta):
        rval = r_val_theo[i]
        theta = thetas[i]
        eq[i] = (rval + np.square(np.sin(theta))) * grad_f_plane[i,0] + (rval + np.square(np.cos(theta))) * grad_f_plane[i,1] - np.cos(theta) * np.sin(theta) * grad_f_plane[i,2]
    return(eq)

def compute_d2f_s_a(S, coeff_mini, powers):
    X_mini = dev_min(S)[:,:3]
    n_data = X_mini.shape[0]
    n_coeff = len(powers)
    degree = np.sum(powers[0])
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    dP_a = np.zeros((n_data, n_coeff))
    for i in range(n_coeff):
        dP_a[:,i] = np.prod(X_mini ** powers[i], axis=1)
    dP_s = grad_polyN_2d(S, coeff_grad, powers_grad)[:,[0,1,3]]
    d2P_s_a = np.zeros((n_data, 3, n_coeff))
    #HERE DDEV/DSIGMA = 1 FOR SX SY SXY
    for i in range(3):
        for j in range(n_coeff):
            d2P_s_a[:,i,j] = coeff_grad[i, j] * np.prod(X_mini ** powers_grad[i, j], axis = 1)
    d2f_s_a = np.zeros((n_data, 3, n_coeff))
    p = polyN_2d(S, coeff_mini, powers)
    k = 2/degree
    for i in range(3):
        for j in range(n_coeff):
            d2f_s_a[:,i,j] = k * np.float_power(p, k - 2) * (d2P_s_a[:,i,j] * p + (k - 1) * dP_a[:,j] * dP_s[:, i])
    return(d2f_s_a)

def grad_rval_ut(thetas, r_val_theos, coeff_mini, powers):
    n_theta = len(thetas)
    S = np.zeros((n_theta, 6))
    v2 = np.zeros((n_theta, 3))
    for i in range(n_theta):
        theta = thetas[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        v2[i] = np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])
    
    n_coeff = len(powers)
    d2f_s_a = compute_d2f_s_a(S, coeff_mini, powers)
    deq_a = np.zeros((n_theta, n_coeff))
    for i in range(n_theta):
        theta = thetas[i]
        for j in range(n_coeff):
            deq_a[i,j] = np.dot(d2f_s_a[i,:,j], v2[i]) + r_val_theos[i] * (d2f_s_a[i,0,j] + d2f_s_a[i,1,j])
    #print(np.linalg.norm(deq_a))
    return(deq_a)


def optiCoeff_polyN_mini(df, degree, weight_exp, weight_rval):
    """
        Returns the optimized coefficients of polyN on experimental data (UTs) and virtual data from a protomodel
        Input :
            - df : pd.Dataframe, must contain columns ["d11", "d22", "s12", "s13", "s23"] of yield stresses points. And if available ["Rval"] with ["LoadAngle"].
            - degree : integer, degree of polyN
            - weight_exp : float, weight for the experimental data
            - weight_rval : float, weight for the rvalue data
        Output :
            - coeff : ndarray of shape (nmon), coeff of polyN model
    """
    data = df[["s11", "s22", "s33", "s12", "s13", "s23"]].values
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X_stress = polyN.fit_transform(dev_min(data)[:,[0,1,2]])

    powers = polyN.powers_
    nmon_abq = len(powers)

    selected_indices = []
    for i in range(nmon_abq):
        k,l,m  = powers[i]
        if (k + l + m == degree) and ((m%2==0)) :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X_stress = X_stress[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2]))
    powers = powers.T[sorted_indices]
    X_stress = X_stress[:,sorted_indices]
    nmon = len(powers)
    ndata = len(data)
    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    ndata_rval = np.count_nonzero(index_rval)
    r_vals = df["Rval"].iloc[index_rval].values
    thetas = df["LoadAngle"].iloc[index_rval].values

    weight_s = np.where(df["Type"] == "e", 1 ,0)
    n_data_exp = np.sum(weight_s)
    weight_s = np.where(weight_s == 1, weight_exp / n_data_exp, (1 - weight_exp) / (ndata - n_data_exp))
    weight_s = (1 - weight_rval) * weight_s / np.sum(weight_s)
    weight_r = np.ones(ndata_rval)
    weight_r = weight_rval *  weight_r / np.sum(weight_r)

    def J(a):
        b = np.ones(len(a) + 1)
        b[1:] = a
        res = 0.5 * (np.sum(weight_s * np.square(f_min_squared(data, b, powers) - 1)))
        return(res)

    def grad_J(a):
        b = np.ones(len(a) + 1)
        b[1:] = a
        n_coeff = len(a)
        grad = np.zeros(n_coeff)
        for i in range(len(a) - 2):
            grad[i] = np.sum(weight_s * (f_min_squared(data, b, powers) - 1) * (2/degree) * X_stress[:, i + 1] * np.float_power(polyN_2d(data, b[:-2], powers),(2/degree - 1)))
        grad[-2] = np.sum(weight_s * (f_min_squared(data, b, powers) - 1) * 2 * np.square(data[:,4]))
        grad[-1] = np.sum(weight_s * (f_min_squared(data, b, powers) - 1) * 2 * np.square(data[:,5]))

        return(grad)
    
    a0 = np.ones(nmon + 2 - 1) * 4
    opt = scipy.optimize.minimize(J, x0=a0, jac=grad_J)
    return(opt.x)

def optiCoeff_polyN_mini2(df, degree, weight_exp, weight_rval):
    """
        Returns the optimized coefficients of polyN on experimental data (UTs) and virtual data from a protomodel
        Input :
            - df : pd.Dataframe, must contain columns ["d11", "d22", "s12", "s13", "s23"] of yield stresses points. And if available ["Rval"] with ["LoadAngle"].
            - degree : integer, degree of polyN
            - weight_exp : float, weight for the experimental data
            - weight_rval : float, weight for the rvalue data
        Output :
            - coeff : ndarray of shape (nmon), coeff of polyN model
    """
    data = df[["s11", "s22", "s33", "s12", "s13", "s23"]].values
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X_stress = polyN.fit_transform(dev_min(data)[:,[0,1,2]])

    powers = polyN.powers_
    nmon_abq = len(powers)

    selected_indices = []
    for i in range(nmon_abq):
        k,l,m  = powers[i]
        if (k + l + m == degree) and ((m%2==0)) :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X_stress = X_stress[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2]))
    powers = powers.T[sorted_indices]
    X_stress = X_stress[:,sorted_indices]
    nmon = len(powers)
    ndata = len(data)
    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    ndata_rval = np.count_nonzero(index_rval)
    r_vals = df["Rval"].iloc[index_rval].values
    thetas = df["LoadAngle"].iloc[index_rval].values

    weight_s = np.where(df["Type"] == "e", 1 ,0)
    n_data_exp = np.sum(weight_s)
    weight_s = np.where(weight_s == 1, weight_exp / n_data_exp, (1 - weight_exp) / (ndata - n_data_exp))
    weight_s = (1 - weight_rval) * weight_s / np.sum(weight_s)
    weight_r = np.ones(ndata_rval)
    weight_r = weight_rval *  weight_r / np.sum(weight_r)

    def J(a):
        b = np.ones(len(a) + 1)
        b[1:] = a
        eq_rval = eq_rval_ut(thetas, r_vals, b[:-2], powers)
        res = 0.5 * (np.sum(weight_s * np.square(f_min_squared(data, b, powers) - 1)) + np.sum(weight_r * np.square(eq_rval)))
        return(res)

    def grad_J(a):
        b = np.ones(len(a) + 1)
        b[1:] = a
        n_coeff = len(a)
        grad = np.zeros(n_coeff)
        grad_rval = grad_rval_ut(thetas, r_vals, b[:-2], powers)
        for i in range(len(a) - 2):
            c = np.sum(weight_r * eq_rval_ut(thetas, r_vals, b[:-2], powers) * grad_rval[:,i + 1])
            d = np.sum(weight_s * (f_min_squared(data, b, powers) - 1) * (2/degree) * X_stress[:, i + 1] * np.float_power(polyN_2d(data, b[:-2], powers),(2/degree - 1)))
            print(i, c, d)
            grad[i] = d + c
        grad[-2] = np.sum(weight_s * (f_min_squared(data, b, powers) - 1) * 2 * np.square(data[:,4]))
        grad[-1] = np.sum(weight_s * (f_min_squared(data, b, powers) - 1) * 2 * np.square(data[:,5]))
        #print(np.linalg.norm(grad))
        return(grad)
    
    if degree==2:
        a0 = np.array([-1, 1, 3, 3, 3])
    if degree==4:
        a0 = np.array([-2, 3, -2, 1, 6, -6, 6, 9, 3, 3])
    if degree==6:
        a0 = np.ones(nmon + 2 - 1)

    opt = scipy.optimize.minimize(J, x0=a0, jac=grad_J)
    return(opt.x)

if gen_v_data :
    print("Generation of virtual data")
    export_virtual_data(protomodel, material, nb_virtual_pt)
    print("Generation ended")

if gen_e_data:
    print("Processing experimental data")
    export_exp_data(material)
    print("Processing ended")

df = readData_2d(material, protomodel)
data = df[["s11", "s22", "s33", "s12", "s13", "s23"]].values
index_rval = np.where(df["Rval"]< 0.00001, False, True)
ndata_rval = np.count_nonzero(index_rval)
r_vals = df["Rval"].iloc[index_rval].values
thetas = df["LoadAngle"].iloc[index_rval].values
coeff_temp = optiCoeff_polyN_mini2(df, degree, weight_exp, weight_rval)
coeff = np.ones(len(coeff_temp) + 1)
coeff[1:] = coeff_temp
print(coeff)
plt.scatter(thetas,rval_ut(thetas,coeff[:-2], powers) , marker="x", label="model")
plt.scatter(thetas, r_vals, marker="x", label="exp")
plt.legend()
plt.show()

# In[209]:


"""----------------------------------------------------FIXING PARAMETERS-----------------------------------------------------------------------------"""



#