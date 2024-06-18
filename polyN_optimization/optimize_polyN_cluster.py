
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


# In[201]:


file_dir = os.path.dirname(os.path.abspath(__file__))  
sep = os.sep
exec_dir = file_dir + sep + "running"


# In[202]:


"---------------------------------------------------------------- PARAMETERS ---------------------------------------------------------------------------------------------"






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
def readData(material, protomodel):
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
    df["d11"] = (2/3) * df["s11"] - (1/3) * df["s22"] - (1/3) * df["s33"]
    df["d22"] = - (1/3) * df["s11"] + (2/3) * df["s22"] - (1/3) * df["s33"]

    return(df)




# In[205]:


"""------------------------ MODEL AND SCALER LOADING -------------------------------------------------"""




"""------------------------------------POLYN PARAMETERS---------------------------------------------------------"""

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
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X = polyN.fit_transform(x0)
    powers = polyN.powers_
    nmon_degree = len(powers)

    selected_indices = []
    for i in range(nmon_degree):
        k,l,m,n,p = powers[i]
        if (k + l + m + n + p == degree) and (((m==n) and (n==p)) or ((m%2==0) and (n%2==0) and (p%2==0))) :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
    powers = powers.T[sorted_indices]
    X = X[:,sorted_indices]


    return(powers)

#degree = int(p["degree"])
#powers, nmon = get_param_polyN(degree)
# In[207]:


""" ---------------------------------------------------POLYN DEFINITION------------------------------------------------------------------------------"""
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

def polyN(S, coeff, powers):
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
    nmon = len(powers)
    for i in range(nmon) :
        p = powers[i]
        res = res + coeff[i] * np.prod(X ** p, axis=1)
    return(res)

"""---------------------------------------------- GRADIENT AND HESSIAN OF POLYN DEFINITION -----------------------------------------------------------"""

#True we lose time here but ok

def jac_dev():
    """
        Returns the jacobian of the deviatoric operator
        Output :
            - jac : ndarray of shape (5,6)
    """
    jac = np.zeros((5, 6))
    jac[0] = np.array([2/3, -1/3, -1/3, 0, 0, 0])
    jac[1] = np.array([-1/3, 2/3, -1/3, 0, 0, 0])
    jac[2] = np.array([0, 0, 0, 1, 0, 0])
    jac[3] = np.array([0, 0, 0, 0, 1, 0])
    jac[4] = np.array([0, 0, 0, 0, 0, 1])
   
    return(jac)

def jac_polyN_param(coeff, powers):
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

    coeff_grad = np.zeros((5, coeff.shape[0]))
    powers_grad = np.zeros((5, powers.shape[0], powers.shape[1]))
    for i in range(5):
        coeff_grad[i] = coeff * powers[:,i]
        subs = np.zeros(powers_grad[i].shape)
        subs[:,i] = -1
        powers_grad[i] = np.maximum(0,powers + subs)
    
    return(coeff_grad, powers_grad)

def grad_polyN(S, coeff_grad, powers_grad):
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
    D = dev(S)
    X = D[:,[0, 1, 3, 4, 5]]

    grad_polyN = np.zeros(6)
    grad_f = np.zeros((len(X),5))

    nmon = len(powers_grad[0])

    for i in range(5):
        for j in range(nmon):
            p = powers_grad[i][j]
            grad_f[:,i] = grad_f[:,i] + coeff_grad[i][j] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev()

    grad_polyN = np.dot(grad_f, jac_d)
    return(grad_polyN)

def hessian_polyN_param(coeff_grad, powers_grad):
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
    coeff_hessian = np.zeros((5, coeff_grad.shape[0], coeff_grad.shape[1]))
    powers_hessian = np.zeros((5, powers_grad.shape[0], powers_grad.shape[1], powers_grad.shape[2]))

    for i in range(5):
        for j in range(5):
            coeff_hessian[i][j] = coeff_grad[i] * powers_grad[i,:,j]
            subs = np.zeros(powers_hessian[i][j].shape)
            subs[:,j] = -1
            powers_hessian[i][j] = np.maximum(0,powers_grad[i] + subs)
    
    return(coeff_hessian, powers_hessian)

def hessian_polyN(S, coeff_hessian, powers_hessian):
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
    D = dev(S)
    X = D[:,[0, 1, 3, 4, 5]]

    hessian_polyN = np.zeros((6,6))
    jac_grad_f = np.zeros((len(X), 5, 5))
    nmon = len(powers_hessian)

    for i in range(5):
        for j in range(5):
            for k in range(nmon):
                p = powers_hessian[i][j][k]
                jac_grad_f[:,i,j] = jac_grad_f[:,i,j] + coeff_hessian[i][j][k] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev(S)
    hessian_polyN = np.transpose(np.dot(jac_d.T,np.dot(jac_grad_f[:], jac_d)), (1, 2, 0))
    return(hessian_polyN)


# In[208]:


""" ---------------------------------------------PARAMETERS OPTIMIZATION-----------------------------------------------------------------------------"""




def optiCoeff_polyN(df, degree, weigth_exp, weigth_rval):
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

    data = df[["d11", "d22", "s12", "s13", "s23"]].values
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X_stress = polyN.fit_transform(data)
    powers = polyN.powers_
    nmon_abq = len(powers)

    selected_indices = []
    for i in range(nmon_abq):
        k,l,m,n,p = powers[i]
        if (k + l + m + n + p == degree) and (((m==n) and (n==p)) or ((m%2==0) and (n%2==0) and (p%2==0))) :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X_stress = X_stress[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
    powers = powers.T[sorted_indices]
    X_stress = X_stress[:,sorted_indices]
    nmon = len(powers)

    ndata = len(data)
    #RVALUE
    #TO UNDERSTAND HOW R-VALUES ARE CALCULATED, REFER TO ARTICLE OF SOARE 2023 ABOUT HOMOGENEOUS YIELD FUNCTIONS
    index_rval = np.where(df["Rval"]< 0.00001, False, True) #RVALUE SIGNIFICANT IF NOT ZERO
    ndata_rval = np.count_nonzero(index_rval)
    dX_stress_dev = np.zeros((ndata_rval, 5, nmon))
    X_stress_rval = data[index_rval]
    param = df[["Rval", "LoadAngle"]].iloc[index_rval].values
    coeff_mon = np.ones(nmon)
    coeff_dmon, powers_dmon = jac_polyN_param(coeff_mon, powers)

    for j in range(5):
        for k in range(nmon):
            p = powers_dmon[j, k]
            dX_stress_dev[:, j, k] = dX_stress_dev[:, j, k] + coeff_dmon[j, k] * np.prod(X_stress_rval ** p, axis = 1)

    dX_stress_stress = np.zeros((ndata, 6, nmon))

    jac_d = jac_dev()

    for i in range(ndata_rval):
        dX_stress_stress[i] = np.dot(jac_d.T, dX_stress_dev[i])

    X_rval = np.zeros((ndata_rval, nmon))
    for i in range(ndata_rval):
        rval, theta = param[i]
        v = np.array([rval + np.square(np.sin(theta)), rval + np.square(np.cos(theta)), 0,  -np.cos(theta) * np.sin(theta), 0, 0])
        X_rval[i] = np.dot(v, dX_stress_stress[i])
    
    weigth_s = np.where(df["Type"] == "e", weigth_exp, 1-weigth_exp)
    weigth_r = np.ones(ndata_rval) * weigth_rval

    weigth = np.concatenate((weigth_s, weigth_r), axis = 0)
    v = np.concatenate((X_stress, X_rval), axis = 0)
    d = np.concatenate((np.ones(ndata), np.zeros(ndata_rval)))
    M = np.zeros((nmon, nmon))

    for j in range(v.shape[0]):
        vj = np.expand_dims(v[j], axis=0)
        M = M + weigth[j] * np.dot(vj.T, vj)
    M = M / 2

    V = np.zeros(v.shape[1])
    
    for j in range(v.shape[0]):
        vj = v[j]
        V = V + weigth[j] * d[j] * vj

    D = np.sum(weigth * np.square(d)) / 2

    def J(a):
        return np.dot(np.dot(a, M), a) - np.dot(V, a) + D

    def Grad_J(a):
        return 2 * np.dot(M, a) - V
    
    if degree==2:
        C = np.array([3, 3, 3, 3, 3, 3], dtype=float)
    if degree==4:
        C = np.array([9, 18, 27, 18, 9, 18, 18, 18, 9, 18, 18, 18, 18, 9, 0, 0, 18, 18, 18, 18, 18, 9], dtype=float)
    if degree==6:
        C = np.array([27, 81, 162, 189, 162, 81, 27, 81, 162, 243, 162, 81, 81, 81, 81, 27, 81, 162, 243, 162, 81, 162, 162, 162, 81, 81, 81, 81, 81, 27, 0, 0, 0, 0, 81, 162, 243, 162, 81, 162, 162, 162, 81, 162, 162, 162, 162, 81, 81, 81, 81, 81, 81, 27,])

    opt = scipy.optimize.minimize(J, C, method='BFGS', jac=Grad_J)
    a_test = opt.x
    return(a_test)


    #We work with scaled data to optimize using neural network's gradient and loss' gradient in the same space of coefficients


# In[209]:


"""----------------------------------------------------FIXING PARAMETERS-----------------------------------------------------------------------------"""

a1 = 3
a2 = 3
a3 = 3
a4 = 3
a5 = 3
a6 = 3
b1 = 9
b2 = 18
b3 = 27
b4 = 18
b5 = 9
b6 = 18
b7 = 18
b8 = 18
b9 = 9
b10 = 18
b11 = 18
b12 = 18
b13 = 18
b14 = 9
b15 = 0
b16 = 0
b17 = 18
b18 = 18
b19 = 18
b20 = 18
b21 = 18
b22 = 9

coeff_deg2 = np.array([a1, a2, a3, a4, a5, a6])
coeff_deg4 = np.array([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, b12, b13, b14, b15, b16, b17, b18, b19, b20, b21, b22])


def adapt_coeff(adapt, degree, coeff):
    """
        Allow to modify the coefficients if the user wants it.
        Input :
            - adapt : bool, 1 if coeff has to be modified
            - degree : integer
            - coeff : ndarray of shape (nmon), coeff of polyN model
    """
    if adapt :
        if degree == 2:
            return coeff_deg2
        elif degree == 4:
            return coeff_deg4
    return coeff

# In[210]:


"""---------------------------------------------------------CALIBRATE SWIFT OR VOCE FUNCTION----------------------------------------------------------------------------"""

def optiCoeff_pflow(law, coeff_polyN, material, degree, powers):
    """
        Returns the a, b, c coefficients of the hardening law defined by the user and the young modulus of the material.
        Input :
            - law : string, law in ["swift", "voce"]
            - coeff_polyN : ndarray of shape (nmon,), polyN coefficients
    """
    
    def f(S):
        return np.power(polyN(S, coeff_polyN, powers), 1/degree)

    foldername = file_dir + sep + "calibration_data" + sep + material
    filename_out = "data_plasticlaw.csv"
    filepath = foldername + sep + filename_out

    df_law = pd.read_csv(filepath)

    ymod = df_law.iloc[0, df_law.columns.get_loc("YoungModulus")]
    ep = df_law["PlasticStrain"].values

    sigma = df_law["PlasticStress"].values
    S = np.expand_dims(sigma, axis = 1)
    S = np.concatenate((S, np.zeros((sigma.shape[0],5))), axis = 1)

    sig0 = f(S[0])[0]
    print(sig0)
    
    if law == "swift" :

        if(0):

            def J(x):
                a, b, c = x
                return np.sum(np.square(a * np.float_power((b + ep),c) - f(S)))
            
            def Grad_J(x):
                a, b, c = x
                rho = a * np.power((ep + b), c) - f(S)
                da = 2 * np.sum(np.float_power((ep + b), c) * rho)
                db = 2 * np.sum(c * a * np.float_power((ep + b), (c - 1)) * rho)
                dc = 2 * np.sum(a * np.log(ep + b) * np.float_power((ep + b), c) * rho)
                return(np.array([da, db, dc]))

            a0, b0, c0 = 1200, 0.4, 0.2
            x0 = np.array([a0, b0, c0])
            opt = scipy.optimize.minimize(J, x0, method="SLSQP", jac=Grad_J)

            a, b, c = opt.x

        if(1):
            #Function with a smooth gradient and enforcing a = sig0
            a = sig0

            def J(x):
                b, c = x
                return np.sum(np.square(np.float_power((1 + b * ep),c) - f(S)/a))
            
            def Grad_J(x):
                b, c = x
                rho = np.float_power((1 + b * ep), c) - f(S)/a
                db = 2 * np.sum(c * ep * np.float_power((1 + b * ep), c - 1) * rho)
                dc = 2 * np.sum(np.log(1 + ep * b) * np.float_power((1 + ep * b), c) * rho)
                return(np.array([db, dc]))
        
            b0, c0 = 1, 1
            x0 = np.array([b0, c0])
            print(a,b0,c0)
            opt = scipy.optimize.minimize(J, x0, method="SLSQP", jac=Grad_J)

            b, c = opt.x
            a, b, c = a / np.float_power(1/b,c), 1/b, c
            
            print(a,b,c)

    
    if law == "voce":

        a = sig0

        def J(x):
            b, c = x
            return np.sum(np.square(a - b * (1 - np.exp(- c * ep)) - f(S)))

        def Grad_J(x):
            b, c = x
            rho = a - b * (1 - np.exp(- c * ep)) - f(S)
            db = - 2 * np.sum((1 - np.exp(- c * ep)) * rho)
            dc = - 2 * np.sum(ep * b * np.exp(- c * ep) * rho)
            return np.array([db, dc])

        b0, c0 = 1, 1
        x0 = np.array([b0, c0])
        opt = scipy.optimize.minimize(J, x0, method="BFGS", jac=Grad_J,)

        b, c = opt.x
    
    return(a, b, c, ymod)




# In[211]:


"""--------------------------------------CHECK OF COEFFICIENTS----------------------------------------------"""

def check_coeff(a, b, c, law, material):
    """
        Plot the hardening law found against the experimental plastic curve (plastic stress according to plastic strain)
        Input :
            - a : float
            - b : float
            - c : float
            - law : string
    """
    foldername = file_dir + sep + "calibration_data" + sep + material
    filename_out = f"data_plasticlaw_{material}.csv"
    filepath = foldername + sep + filename_out

    df_law = pd.read_csv(filepath)

    ymod = df_law.iloc[0, df_law.columns.get_loc("YoungModulus")]
    ep = df_law["PlasticStrain"].values

    sigma = df_law["PlasticStress"].values

    if(1):
        def swift(ep):
            return(a * np.power((b + ep), c))
    if(0):
        def swift(ep):
            return(a * np.power((1 + ep/b), c))

    swift = np.vectorize(swift)

    def voce(ep):
        return(a - b * (1- np.exp(- c * ep)))

    voce = np.vectorize(voce)

    fig, ax = plt.subplots(1,1)

    X = np.linspace(0, np.max(ep), 100)

    if law == "swift" :
        Y = swift(X)
    elif law == "voce":
        Y = voce(X)

    ax.plot(X, Y)
    ax.set_title("Law rule")

    ax.scatter(ep, sigma, marker="x", linewidths=0.1 )
    plt.show()


# In[212]:

    
def plot_check(df, coeff_polyN, powers, material, weigth_exp, savefigyr, nb_virtual_pt, degree, protomodel):
    """
        Plot the yield surface in the sx, sy plane and the yield stresses and r-values according to the loading angle.
        Input :
            - df : pd.Dataframe, must contain columns ["d11", "d22", "s12", "s13", "s23"] of yield stresses points. And if available ["Rval"] with ["LoadAngle"]
            - coeff_polyN : ndarray of shape (nmon,)
    """
    def f(S):
        return(polyN(S, coeff_polyN, powers))

    coeff_grad, powers_grad = jac_polyN_param(coeff_polyN, powers)

    def grad(S):
        return grad_polyN(S, coeff_grad, powers_grad)


    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    fig, ax1 = plt.subplots()

    ys_exp = df_exp["YieldStress"]
    n_data = len(ys_exp)
    theta_exp = df_exp["LoadAngle"].values
    R_exp = df_exp["Rval"].values

    ax1.scatter(theta_exp, ys_exp/sigma0, color="blue", label="YS exp", marker="x")
    ax1.set_xlabel(r"$\theta$[rad]", size=12)
    ax1.set_ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)

    ax1.scatter(theta_exp, R_exp, color="red", label="Rval exp", marker="x")
    #plt.scatter(theta_exp, R_theo, label="Modelisation", marker="x", color="red")

    ntheta=100
    thetas=np.linspace(0, np.pi / 2, ntheta)

    def ys_theta(theta):
        """
            For a given loading angle, returns the norm of the yield stress in a ut_theta test
            Input :
                - theta : float (radians)
            Output :
                - m : float, norm of yield stress
        """
        u = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
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

        return(m)
    
    ys_theta = np.vectorize(ys_theta)
    ys_model = ys_theta(thetas)
    ax1.plot(thetas, ys_model, color="blue", label="YS model")


    def ys_component_theta(theta):
        """
            For a given loading angle, returns the yield stress (sx, sy, sz, sxy, sxz, syz) in a ut_theta test
            Input :
                - theta : float (radians)
            Output :
                - m : ndarray of shape (6,), yield stress
        """
        u = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        ys = 0.1 * u
        lamb = 1.1
        eps = 1e-7

        while f(ys) < 1:
            ys = ys * lamb

        s = 0.1 * u
        e = ys
        m = (s + e) / 2
        res = f(m) - 1

        while (np.abs(res) > eps):
            if res > 0:
                e = m
            else:
                s = m
            m = (s + e) / 2
            res = f(m) - 1

        return(m)


    #TO UNDERSTAND HOW R-VALUES ARE CALCULATED, REFER TO ARTICLE OF SOARE 2023 ABOUT HOMOGENEOUS YIELD FUNCTIONS
    X = np.array([ys_component_theta(theta) for theta in thetas])
    gradX = grad(X)
    vs = np.array([- np.sin(thetas), np.cos(thetas), np.zeros(ntheta)]).T
    vs = np.expand_dims(vs, axis=1)
    V = np.matmul(np.transpose(vs, (0, 2, 1)), vs).reshape(ntheta, 9)[:, [0, 4, 8, 1, 3, 2]]
    R_model = - np.sum(gradX * V, axis = 1)/(gradX[:,0] + gradX[:,1])

    ax1.plot(thetas, R_model, label="Rval model", color="red")
    ax1.set_title(f"Poly{degree} model on {material} : Yield Stresses and R-values from UT tests", size=12)
    plt.legend()
    plt.grid(1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if savefigyr :
        foldername_out = file_dir + sep + "plots" + sep + material
        if not os.path.exists(foldername_out):
            os.makedirs(foldername_out)
        filename = f"ysrval_{material}_poly{degree}_{weigth_exp}_{nb_virtual_pt}_{protomodel}.png"
        filepath = foldername_out + sep + filename
        plt.savefig(filepath)
    else:
        plt.show()

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

def plot_planestress(material, coeff, powers, savefigplane, weigth_exp, nb_virtual_pt, degree, protomodel):
    zs = np.linspace(0, 1, 10)
    mises_plot = False
    fig, ax = plt.subplots()


    sx = np.linspace(-1.5, 1.5, 100)
    sy = np.linspace(-1.5, 1.5, 100)
    sx, sy = np.meshgrid(sx,sy)

    for z in zs:
        def mises_plane(x,y):
            return(mises(np.array([x,y,0,z,0,0])))

        def f(x,y):
            return(polyN(np.array([x,y,0,z,0,0]), coeff, powers))

        f = np.vectorize(f)
        mises_plane = np.vectorize(mises_plane)

        ys_polyN = f(sx, sy)
        ys_mises = mises_plane(sx, sy)

        # Create the contour plot
        cs1 = ax.contour(sx, sy, ys_polyN, levels=[1], colors='blue', linewidths=1)
        if mises_plot:
            cs2 = ax.contour(sx, sy, ys_mises, levels=[1], colors='red', linewidths=1)

    # Set labels
    ax.set_xlabel(r"$\sigma_{xx}/\sigma_0$[-]", size=12)
    ax.set_ylabel(r"$\sigma_{yy}/\sigma_0$[-]", size=12)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(1)
    
    nm1, labels = cs1.legend_elements()
    if mises_plot:
        nm2, labels = cs2.legend_elements()


    ax.set_title(rf'{material} Yield surface in the $\sigma_{{xx}},\sigma_{{yy}}$ plane', size=12)
    plt.legend(nm1, ["polyN"])
    if mises_plot:
        plt.legend(nm2, ["Mises"])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    if savefigplane :
        foldername_out = file_dir + sep + "plots" + sep + material
        if not os.path.exists(foldername_out):
            os.makedirs(foldername_out)
        filename = f"planexy_{material}_poly{degree}_{weigth_exp}_{nb_virtual_pt}_{protomodel}.png"
        filepath = foldername_out + sep + filename
        plt.savefig(filepath)
    else:
        plt.show()




"""----------------------------------------------------TESTING FUNCTIONS (PLOT & CONVEXITY)-----------------------------------------------------------------"""
## WRAP FUNCTION BEFORE USING PLOT IMPLICIT
def plot_implicit(yf, bbox=(-1.5,1.5)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 50) # resolution of the contour
    B = np.linspace(xmin, xmax, 10) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    print("Début du plot sur XY")
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = yf(X,Y,z) - 1
        cset = ax.contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for this value of z
    print("Fin du plot sur XY")
    print("Début du plot sur XZ")
    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = yf(X,y,Z) - 1
        cset = ax.contour(X, Y+y, Z, [y], zdir='y')

    print("Fin du plot sur XZ")
    print("Début du plot sur YZ")
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = yf(x,Y,Z) - 1
        cset = ax.contour(X+x, Y, Z, [x], zdir='x')
    print("Fin du plot sur YZ")
    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)

    plt.show()


# In[213]:


def plot_implicit_coeff(coeff, powers, bbox=(-1.5,1.5)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''

    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig, ax = plt.subplots(nrows = 1, ncols = 3, subplot_kw={"projection":"3d"})
    A = np.linspace(xmin, xmax, 50) # resolution of the contour
    B = np.linspace(xmin, xmax, 20) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    def f(S):
        return polyN(S, coeff, powers)

    for i in range(3):

        if i == 0:
            f_plane = lambda x, y, z : f(np.array([x, y, 0, z, 0, 0]))
        elif i == 1:
            f_plane = lambda x, y, z : f(np.array([x, y, 0, 0, z, 0]))
        else :
            f_plane = lambda x, y, z : f(np.array([x, y, 0, 0, 0, z]))
        
        f_plane = np.vectorize(f_plane)

        print("Début du plot sur XY")
        for z in B: # plot contours in the XY plane
            X,Y = A1,A2
            Z = f_plane(X,Y,z) - 1
            ax[i].contour(X, Y, Z+z, [z], zdir='z')
            # [z] defines the only level to plot for this contour for this value of z
        print("Fin du plot sur XY")
        print("Début du plot sur XZ")
        for y in B: # plot contours in the XZ plane
            X,Z = A1,A2
            Y = f_plane(X,y,Z) - 1
            ax[i].contour(X, Y+y, Z, [y], zdir='y')

        print("Fin du plot sur XZ")
        print("Début du plot sur YZ")
        for x in B: # plot contours in the YZ plane
            Y,Z = A1,A2
            X = f_plane(x,Y,Z) - 1
            ax[i].contour(X+x, Y, Z, [x], zdir='x')
        print("Fin du plot sur YZ")
        # must set plot limits because the contour will likely extend
        # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
        # to encompass all values in the contour.
        ax[i].set_zlim3d(zmin,zmax)
        ax[i].set_xlim3d(xmin,xmax)
        ax[i].set_ylim3d(ymin,ymax)

    plt.show()


# In[214]:





# In[217]:


"""-------------------------------------------------OUTPUT---------------------------------------------------------------------------------------------"""
def write_coeff_user(coeff, protomodel, degree, material, nb_virtual_pt, powers):
    filename = "{}_deg{}_{}.txt".format(material, degree, protomodel)
    foldername = file_dir + sep + "for_user" 
    filepath = foldername + sep + filename

    n = len(coeff)
    nmon = len(powers)
    with open(filepath, "w") as file:
        file.write("#Coefficients poly{} for {} based on {} points from the {} protomodel\n".format(degree, material, nb_virtual_pt, protomodel))
        file.write("Number of coefficients : {}\n".format(n))
        file.write("np.array([")
        n0 = 0
        n_p = 0
        while n0 < nmon:
            for m in range(0, degree + 1):
                for l in range(0, degree + 1 - m):
                    for k in range(0, degree + 1 - m - l):
                        for j in range(0, degree + 1 - m - l - k):
                            i = degree - m - l - k - j
                            i0, j0, k0, l0, m0 = powers[n0]
                            if (i==i0) and (j==j0) and (k==k0) and (l==l0) and (m==m0):
                                file.write("{}, ".format(coeff[n0]))
                                n0 = n0 + 1
        file.write("])")
            
def write_coeff_abq(coeff, a, b, c, ymod, enu, nmon, protomodel, degree, material, law, density, powers):
    filename = "{}_abq_deg{}_{}_{}.inp".format(material, degree, law, protomodel)
    foldername = file_dir + sep + "running"
    filepath = foldername + sep + filename

    with open(filepath, "w") as file:
        file.write("*USER MATERIAL, constants={}\n".format(7 + nmon))
        file.write("{}, {}, {}, {}, {}, {}, {}, ".format(ymod, enu, a, b, c, degree, nmon))
        n0 = 0
        while n0 < nmon:
            for m in range(0, degree + 1):
                for l in range(0, degree + 1 - m):
                    for k in range(0, degree + 1 - m - l):
                        for j in range(0, degree + 1 - m - l - k):
                            i = degree - m - l - k - j
                            i0, j0, k0, l0, m0 = powers[n0]
                            if (i==i0) and (j==j0) and (k==k0) and (l==l0) and (m==m0):
                                file.write("{}, ".format(coeff[n0]))
                                n0 = n0 + 1
                                if (n0 + 7) % 8 == 0:
                                    file.write("\n")
                            
                            
        file.write("\n")
        file.write("*DENSITY\n")
        file.write("{}".format(density))
 

# In[218]:

"""-------------------------------------------------- TO RUN ----------------------------------------------------------- """



def first_opti():

    p = read_param()

    material = p["material"]
    gseed = int(p["gseed"])
    enu = float(p["enu"])
    density = float(p["density"])
    nb_virtual_pt = int(p["nb_virtual_pt"])
    degree = int(p["degree"])
    weigth_exp = float(p["weigth_exp"])
    weigth_rval = float(p["weigth_rval"])
    protomodel = p["protomodel"]
    law = p["law"]

    gen_v_data = int(p["gen_v_data"])
    gen_e_data = int(p["gen_e_data"])

    opti = int(p["opti"])
    loadcoeff = int(p["loadcoeff"])
    adapt = int(p["adapt"])

    export_coeff_abq = int(p["export_coeff_abq"])
    export_coeff_user = int(p["export_coeff_user"])


    plot = int(p["plot"])
    savefigyr = int(p["savefigyr"])
    savefigplane = int(p["savefigplane"])
    savecoeff = int(p["savecoeff"])

    powers = get_param_polyN(degree)
    nmon = len(powers)

    copy_lab_data(material)

    if gen_v_data :
        print("Generation of virtual data")
        export_virtual_data(protomodel, material, nb_virtual_pt)
        print("Generation ended")

    if gen_e_data:
        print("Processing experimental data")
        export_exp_data(material)
        print("Processing ended")


    get_hardening_law(material)
    df = readData(material, protomodel)
    nb_virtual_pt = len(df[df["Type"] == "v"])

    if degree==2:
        C = np.array([3, 3, 3, 3, 3, 3], dtype=float)
    if degree==4:
        C = np.array([9, 18, 27, 18, 9, 18, 18, 18, 9, 18, 18, 18, 18, 9, 0, 0, 18, 18, 18, 18, 18, 9], dtype=float)
    if degree==6:
        C = np.array([27, 81, 162, 189, 162, 81, 27, 81, 162, 243, 162, 81, 81, 81, 81, 27, 81, 162, 243, 162, 81, 162, 162, 162, 81, 81, 81, 81, 81, 27, 0, 0, 0, 0, 81, 162, 243, 162, 81, 162, 162, 162, 81, 162, 162, 162, 162, 81, 81, 81, 81, 81, 81, 27,])

    if opti:
        coeff = optiCoeff_polyN(df, degree, weigth_exp, weigth_rval)
    elif loadcoeff:
        coeff = np.load("polyN_coeff.npy")
    else :
        coeff = C

    coeff = adapt_coeff(adapt, degree, coeff)
    a, b, c, ymod = optiCoeff_pflow(law, coeff, material, degree, powers)
    check_coeff(a,b,c,law, material)

    if savecoeff:
        np.save(file_dir + sep + "polyN_coeff.npy", coeff)
        np.save(file_dir + sep + f"{law}_coeff.npy", np.array([a, b, c, ymod]))

    #plot_check(df, coeff)
    #plot_planestress(material, coeff)

    if plot:
        plot_implicit_coeff(coeff, powers)
    if export_coeff_user:
        write_coeff_user(coeff, protomodel, degree, material, nb_virtual_pt, powers)
    if export_coeff_abq:
        write_coeff_abq(coeff, a, b, c, ymod, enu, nmon, protomodel, degree, material, law, density, powers)