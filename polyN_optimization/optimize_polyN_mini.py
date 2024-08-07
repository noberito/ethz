
#!/usr/bin/env python
# coding: utf-8

# In[200]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize
import os
import sklearn
import sklearn.preprocessing
from sklearn.linear_model import LinearRegression
from read import read_param, readData_2d
from get_calibration_data import export_exp_data, export_virtual_data, analyze_exp_data
from get_hardening_law import get_hardening_law
import shutil
import random


# In[201]:


polyN_dir = os.path.dirname(os.path.abspath(__file__))  
sep = os.sep
exec_dir = polyN_dir + sep + "running"



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
    """
        Adds the experimental data to the results_exp folder concerning the given material

        Input :
            - material : string
    """
    copy_dir = f"{polyN_dir}{sep}results_exp{sep}{material}"
    if not os.path.exists(copy_dir):
        parent_dir = os.path.dirname(polyN_dir)
        lab_data_dir = f"{parent_dir}{sep}lab_data"
        mat_lab_data_dir = f"{lab_data_dir}{sep}{material}_results{sep}DATA"
        shutil.copytree(mat_lab_data_dir, copy_dir)


"""------------------------------------POLYN PARAMETERS---------------------------------------------------------"""

def get_param_polyN_mini(degree):
    """
        Returns polyN minimalistic powers for the given degree

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


#degree = int(p["degree"])
#powers, nmon = get_param_polyN(degree)
# In[207]:


""" ---------------------------------------------------POLYN DEFINITION------------------------------------------------------------------------------"""
def t_linear(S):
    """
        Returns the deviatoric stress.

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D

        Output :
            - D : ndarray of shape (n,6), linear combination stress components in 3D
    
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

def polyN_2d(S, coeff_mini, powers):
    """
        Compute the polyN o t_linear function.

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_mini : ndarray of shape (nmon,), coefficients of the polyN_mini function

        Output :
            - res : float, result of polyN_min(t_linear(S))
    """
    degree = np.sum(powers[0])
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = t_linear(S)[:,[0,1,2]]
    res = np.zeros(len(X))
    nmon = len(powers)
    for i in range(nmon) :
        p = powers[i]
        res = res + coeff_mini[i] * np.prod(X ** p, axis=1)
    return(res)

def polyN_2d_opti(X, coeff_mini, powers):
    """
        Compute the polyN function.
        Input :
            - X : ndarray of shape (n,3) or (3,), stress components in 3D
            - coeff_mini : ndarray of shape (nmon,), coefficients of the polyN_mini function
            - powers : ndarray of shape (nmon, 3), exponents
        Output :
            - res : ndarray of shape (n,), result of polyN_mini
    """
    if X.ndim==1:
        X = np.expand_dims(X, axis=0)
    res = np.zeros(1)
    nmon = len(powers)
    for i in range(nmon) :
        p = powers[i]
        res = res + coeff_mini[i] * np.prod(X ** p, axis=1)
    return(res)

"""---------------------------------------------- GRADIENT AND HESSIAN OF POLYN DEFINITION -----------------------------------------------------------"""

#True we lose time here but ok

def jac_t_linear():
    """
        Returns the jacobian of the linear transformation operator (Ref : Article Soare 2007)

        Output :
            - jac : ndarray of shape (3,6)
    """
    jac = np.zeros((3, 6))
    jac[0] = np.array([1, 0, -1, 0, 0, 0])
    jac[1] = np.array([0, 1, -1, 0, 0, 0])
    jac[2] = np.array([0, 0, 0, 1, 0, 0])
   
    return(jac)

def jac_polyN_2d_param(coeff_mini, powers):
    """
        Compute coeffcients and powers of polyN 2D Jacobian.

        Input :
            - coeff_mini (float ndarray of shape (nmon)) : Coefficients of the polyN 2D function
            ordered by power[1] asc, power[2] asc
            - powers (float ndarray of shape (nmon, 3)) : Powers for each monomial of the 
            PolyN 2D function
        
        Output :
            - coeff_grad (float ndarray of shape (3, nmon)) : coeff_grad[i,j] contains the 
            coefficients of the monomial j derived with respect to dev[i]
            - powers_grad (float ndarray of shape (3, nmon, 3)) : powers_grad[i,j] contains
            the monomial j derived with respect to dev[i]
    """

    coeff_grad = np.zeros((3, len(coeff_mini)))
    powers_grad = np.zeros((3, powers.shape[0], powers.shape[1]))
    for i in range(3):
        coeff_grad[i] = coeff_mini * powers[:,i]
        subs = np.zeros(powers_grad[i].shape)
        subs[:,i] = -1
        powers_grad[i] = np.maximum(0,powers + subs)
    
    return(coeff_grad, powers_grad)

def grad_polyN_2d(S, coeff_grad, powers_grad):
    """
        Returns the Jacobian of polyN2D o t_linear.

        Input :
            - S : float ndarray of shape (n, 6) or (6,), stress components in 3d
            - coeff_grad : float ndarray of shape (3, nmon),coeff_grad[i] contains the 
            coefficients of each monomial derived with respect to dev[i]
            - powers_grad : float ndarray of shape (3, nmon, 3) : powers_grad[i][j] contains
            the monomial j derived with respect to i-th deviatoric component
        
        Output :
            - grad : ndarray of shape (n, 6), grad_polyN[i, j] = dpolyN/dsj of the i-th data
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = t_linear(S)[:,[0,1,2]]

    grad_f = np.zeros((len(X),3))

    nmon = len(powers_grad[0])

    for i in range(3):
        for j in range(nmon):
            p = powers_grad[i][j]
            grad_f[:,i] = grad_f[:,i] + coeff_grad[i][j] * np.prod(X ** p, axis=1)
    
    jac_t = jac_t_linear()

    grad_polyN = np.dot(grad_f, jac_t)
    return(grad_polyN)

def grad_polyN_2d_opti(X, coeff_grad, powers_grad):
    """
        Returns the Jacobian of polyN2D.

        Input :
            - X : float ndarray of shape (n, 3) or (3,), stress components in 3d
            - coeff_grad : float ndarray of shape (3, nmon),coeff_grad[i] contains the 
            coefficients of each monomial derived with respect to var[i]
            - powers_grad : float ndarray of shape (3, nmon, 3) : powers_grad[i][j] contains
            the monomial j derived with respect to i-th component
        
        Output :
            - grad : ndarray of shape (n, 3), grad_polyN[i, j] = dpolyN/dsj of the i-th data
    """
    if X.ndim==1:
        X = np.expand_dims(X, axis=0)

    grad_P = np.zeros((len(X),3))
    nmon = len(powers_grad[0])

    for i in range(3):
        for j in range(nmon):
            p = powers_grad[i][j]
            grad_P[:,i] = grad_P[:,i] + coeff_grad[i][j] * np.prod(X ** p, axis=1)
    
    return(grad_P)

def hessian_polyN_2d_param(coeff_grad, powers_grad):
    """
        Returns coefficients and parameters of polyN2D hessian.

        Input :
            - coeff_grad : float ndarray of shape (5, nmon), coeff_grad[i] contains the 
            coefficients of dpolyN/ddev[i]
            - powers_grad : float ndarray of shape (5, nmon, 5), powers_grad[i,j] contains
            the powers of dmon[j]/ddev[i]
        
        Output :
            - coeff_hess : float ndarray of shape (5, 5, nmon)), coeff_hess[i,j,k] contains
            the coefficients of d2mon[k]/ddev[j]ddev[i] 
            - powers_hess : float ndarray of shape (5, 5, nmon, 5), powers_hess[i,j,k]contains
            the powers of d2mon[k]/ddev[j]ddev[i] 

    """
    coeff_hess = np.zeros((3, coeff_grad.shape[0], coeff_grad.shape[1]))
    powers_hess = np.zeros((3, powers_grad.shape[0], powers_grad.shape[1], powers_grad.shape[2]))

    for i in range(3):
        for j in range(3):
            coeff_hess[i][j] = coeff_grad[i] * powers_grad[i,:,j]
            subs = np.zeros(powers_hess[i][j].shape)
            subs[:,j] = -1
            powers_hess[i][j] = np.maximum(0,powers_grad[i] + subs)
    
    return(coeff_hess, powers_hess)

def hessian_polyN_2d(S, coeff_hess, powers_hess):
    """
        Returns polyN2D o t_linear hessian.

        Input :
            - S : float ndarray of shape (n, 6) or (6,), stress components in 3d
            - coeff_hess : float ndarray of shape (3, 3, nmon)), coeff_hess[i,j,k] contains
            the coefficients of d2mon[k]/ddev[j]ddev[i] 
            - powers_hess : float ndarray of shape (3, 3, nmon, 3), powers_hess[i,j,k]contains
            the powers of d2mon[k]/ddev[j]ddev[i]

        Output :
            - hessian : float ndarray of shape (n, 6, 6), hessian[i,j,k] = dpolyN/ds[k]ds[j] of the ith data pt
    """

    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = t_linear(S)[:,[0,1,2]]

    hessian_polyN = np.zeros((6,6))
    jac_grad_f = np.zeros((len(X), 3, 3))
    nmon = len(powers_hess)

    for i in range(3):
        for j in range(3):
            for k in range(nmon):
                p = powers_hess[i][j][k]
                jac_grad_f[:,i,j] = jac_grad_f[:,i,j] + coeff_hess[i][j][k] * np.prod(X ** p, axis=1)
    
    jac_t = jac_t_linear()
    hessian_polyN = np.transpose(np.dot(jac_t.T,np.dot(jac_grad_f[:], jac_t)), (1, 2, 0))
    return(hessian_polyN)

def hessian_polyN_2d_opti(X, coeff_hess, powers_hess):
    """
        Returns polyN2D hessian.

        Input :
            - X : float ndarray of shape (n, 3) or (3,), stress components in 3d
            - coeff_hess : float ndarray of shape (3, 3, nmon)), coeff_hess[i,j,k] contains
            the coefficients of d2mon[k]/ddev[j]ddev[i] 
            - powers_hess : float ndarray of shape (3, 3, nmon, 3), powers_hess[i,j,k]contains
            the powers of d2mon[k]/dvar[j]dvar[i] 
        Output :
            - hessian : float ndarray of shape (n, 3, 3), hessian[i,j,k] = dpolyN/ds[k]ds[j] of the ith data pt

    """
    if X.ndim==1:
        X = np.expand_dims(X, axis=0)

    hess_P = np.zeros((len(X), 3, 3))
    nmon = len(powers_hess[0][0])
    for i in range(3):
        for j in range(3):
            for k in range(nmon):
                p = powers_hess[i][j][k]
                hess_P[:,i,j] = hess_P[:,i,j] + coeff_hess[i][j][k] * np.prod(X ** p, axis=1)
    return(hess_P)

def f1(S, coeff_mini, powers):
    """
        Returns f1 (polyN o t_linear ** (2/n)).

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_mini : ndarray of shape (nmon,), coefficients of the polyN_mini function
            - powers : ndarray of shape (nmon, 3), exponents

        Output :
            - res : ndarray of shape (n,), result of polyN_min(t_linear(S)) **(2/degree)       
    """
    degree = np.sum(powers[0])
    p = polyN_2d(S, coeff_mini, powers)
    res = np.float_power(p,(2/degree))
    return(res)

def gradf1(S, coeff_mini, powers):
    """
        Returns f1 gradient (grad(polyN o t_linear ** (2/n))).

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_mini : ndarray of shape (nmon,), coefficients of the polyN_mini function
            - powers : ndarray of shape (nmon, 3), exponents
            
        Output :
            - gradf1 : ndarray of shape (n,6), grad of polyN_min(t_linear(S)) **(2/degree)       
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    n = np.sum(powers[0])

    p = polyN_2d(S, coeff_mini, powers)
    gradP = grad_polyN_2d(S, coeff_grad, powers_grad)


    gradf1 = (2/n) * gradP * np.float_power(p, (2/n) - 1)[:, np.newaxis]

    return(gradf1)



def hessf1(S, coeff_mini, powers):
    """
        Returns f1 hessian (grad(polyN o t_linear ** (2/n))).

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_mini : ndarray of shape (nmon,), coefficients of the polyN_mini function
            - powers : ndarray of shape (nmon, 3), exponents
            
        Output :
            - hessf1 : ndarray of shape (n,6,6), hessian of polyN_min(t_linear(S)) **(2/degree)       
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)

    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    coeff_hess, powers_hess = hessian_polyN_2d_param(coeff_grad, powers_grad)

    n = np.sum(powers[0])

    p = polyN_2d(S, coeff_mini, powers)
    gradP = grad_polyN_2d(S, coeff_grad, powers_grad)
    hessP = hessian_polyN_2d(S, coeff_hess, powers_hess)

    hessf1 = (2/n) * hessP * np.float_power(p, (2/n) - 1)[:, np.newaxis, np.newaxis]
    hessf1 = hessf1 + (2/n) * (2/n - 1) * np.einsum('ij,ik->ijk', gradP, gradP) *  np.float_power(p, (2/n) - 2)[:, np.newaxis, np.newaxis]

    return(hessf1)


def f2(S, coeff_oop):
    """
        Returns f2 (out of plane shear terms).

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_oop : ndarray of shape (nmon,), coefficients in front of out of plane shear components

        Output :
            - res : ndarray of shape (n,), k1 * sxz**2 + k2 * syz**2      
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    X = t_linear(S)
    res = coeff_oop[0] * np.square(X[:,3]) + coeff_oop[1] * np.square(X[:,4])
    return(res)

def gradf2(S, coeff_oop):
    """
        Returns f2 gradient (out of plane shear terms).

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_oop : ndarray of shape (nmon,), coefficients in front of out of plane shear components

        Output :
            - gradf2 : ndarray of shape (n,6), 2 * k1 * sxz + 2 * k2 * syz    
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)

    k1 = coeff_oop[0]
    k2 = coeff_oop[1]

    ndata = len(S)
    gradf2 = np.zeros((ndata, 6))

    gradf2[:,4] = 2 * k1 * S[:,4]
    gradf2[:,5] = 2 * k2 * S[:,5]

    return(gradf2)



def hessf2(S, coeff_oop):
    """
        Returns f2 hessian (out of plane shear terms).

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_oop : ndarray of shape (nmon,), coefficients in front of out of plane shear components

        Output :
            - hessf2 : ndarray of shape (n,6,6), 2 * sxz + 2 * syz    
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)

    k1 = coeff_oop[0]
    k2 = coeff_oop[1]

    ndata = len(S)
    hessf2 = np.zeros((ndata, 6, 6))

    hessf2[:,4, 4] = 2 * k1 
    hessf2[:,5, 5] = 2 * k2 

    return(hessf2)


# In[208]:

def f_min_squared(S, coeff, powers):
    """
        Compute the yield function squared.

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function
            - powers : ndarray of shape (nmon, 3), exponents of polyN2D

        Output :
            - res : ndarray of shape (n,), result of polyN_min(t_linear(S)) **(2/degree) + kxz * sxz**2 + kyz * syz**2
    """
    return(f1(S, coeff[:-2], powers) + f2(S, coeff[-2:]))

def grad_f_min_squared(S, coeff, powers):
    """
        Compute the yield function squared gradient

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function
            - powers : ndarray of shape (nmon, 3), exponents of polyN2D

        Output :
            - grad : ndarray of shape(n,6), gradient of polyN_min(t_linear(S)) **(2/degree) + kxz * sxz**2 + kyz * syz**2
    """
    return(gradf1(S, coeff[:-2], powers) + gradf2(S, coeff[-2:]))

def hess_f_min_squared(S, coeff, powers):
    """
        Compute the yield function squared hessian

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function
            - powers : ndarray of shape (nmon, 3), exponents of polyN2D

        Output :
            - hess : ndarray of shape (n,6,6), hessian of polyN_min(t_linear(S)) **(2/degree) + kxz * sxz**2 + kyz * syz**2
    """
    return(hessf1(S, coeff[:-2], powers) + hessf2(S, coeff[-2:]))


""" ---------------------------------------------PARAMETERS OPTIMIZATION-----------------------------------------------------------------------------"""

    
def ys_ut(thetas, coeff, powers):
    """
        Returns the stress value for a uniaxial tensile test in the direction theta
        for theta in the given thetas. Here the yield function is the polyN minimalistic
        function.

        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function
            - powers : ndarray of shape (nmon, 3), powers of the polyN2D function

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
        return(f_min_squared(S, coeff, powers))

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

def rval_ut(thetas, coeff_mini, powers):
    """
        Returns the rvalues for a uniaxial tensile test in the direction theta
        for theta in the given thetas. The yield function here is the polyN minimalistic 
        function.

        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff_mini : ndarray of shape (nmon + 2,), coefficients of the polyN function
            - powers : ndarray of shape (nmon, 3), powers of the polyN function

        Output :
            - rval : ndarray of shape (n_theta,)
    """
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


def eq_rval_ut(thetas, r_val_exp, coeff_mini, powers):
    """
        Returns the values concerning the r-values that need to be equal to 0 in 
        optimization procedure.

        Input :
            - thetas : ndarray of shape (n_theta,), thetas where r_val available
            - r_val_theo : ndarray of shape (n_theta,), r-values available
            - coeff_mini : ndarray of shape (nmon,), coefficients of the poly2D function
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
        
        Output :
            - eq : ndarray of shape (n_theta,), floats that need to be equal to 0
    """
    n_theta = len(thetas)
    S = np.zeros((n_theta, 6))
    v2 = np.zeros((n_theta, 3))

    for i in range(n_theta):
        theta = thetas[i]
        r_val = r_val_exp[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        v2[i] = np.array([r_val + np.square(np.sin(theta)), r_val + np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])

    degree = np.sum(powers[0])
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    grad_P = grad_polyN_2d(S, coeff_grad, powers_grad)[:, [0,1,3]]
    grad_f_plane = np.zeros((n_theta, 3))
    k = 2/degree

    for i in range(n_theta):
        for j in range(3):
            grad_f_plane[i,j] = k * grad_P[i,j] * np.float_power(polyN_2d(S[i], coeff_mini, powers),(k - 1))[0]

    eq = np.zeros(n_theta)

    for i in range(n_theta):
        eq[i] = np.dot(grad_f_plane[i], v2[i])

    return(eq)

def compute_d2f_s_a(S, coeff_mini, powers):
    """
        Returns the yield function derived twice with respect to s and a
        when k1, k2 null. 

        Input :
            - S : ndarray of shape (n,6) or (6,), stress components in 3D
            - coeff_mini : ndarray of shape (nmon,), coefficients of the polyN_mini function
            - powers : ndarray of shape (nmon, 3), exponents
        
        Output :
            - d2f_s_a : ndarray of shape (n, 3, nmon)
    """
    X_mini = t_linear(S)[:,:3]
    n = X_mini.shape[0]
    nmon = len(powers)
    degree = np.sum(powers[0])
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    dP_a = np.zeros((n, nmon))

    for i in range(nmon):
        dP_a[:,i] = np.prod(X_mini ** powers[i], axis=1)
    #dP_a checked

    dP_s = grad_polyN_2d(S, coeff_grad, powers_grad)[:,[0,1,3]]
    #dP_s checked

    d2P_s_a = np.zeros((n, 3, nmon))
    coeff_grad, powers_grad = jac_polyN_2d_param(np.ones(nmon), powers)
    for i in range(3):
        for j in range(nmon):
            d2P_s_a[:,i,j] = coeff_grad[i, j] * np.prod(X_mini ** powers_grad[i, j], axis = 1)

    #d2P_s_a checked

    d2f_s_a = np.zeros((n, 3, nmon))
    p = polyN_2d(S, coeff_mini, powers)
    k = 2/degree
    for i in range(3):
        for j in range(nmon):
            d2f_s_a[:,i,j] = k * np.float_power(p, k - 2) * (d2P_s_a[:,i,j] * p + (k - 1) * dP_a[:,j] * dP_s[:, i])

    return(d2f_s_a)

def grad_rval_ut(thetas, r_vals_exp, coeff_mini, powers):
    """
        Returns the gradient of r-val equations

        Input :
            - thetas : ndarray of shape (n_theta,), thetas where r_val available
            - r_val_exp : ndarray of shape (n_theta,), r-values available
            - coeff_mini : ndarray of shape (nmon,), coefficients of the poly2D function
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
        
        Output :
            - eq : ndarray of shape (n_theta,), floats that need to be equal to 0
    """
    n_theta = len(thetas)
    n_coeff = len(powers)

    S = np.zeros((n_theta, 6))
    v2 = np.zeros((n_theta, 3))

    for i in range(n_theta):
        theta = thetas[i]
        r_val = r_vals_exp[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        v2[i] = np.array([np.square(np.sin(theta)) + r_val, np.square(np.cos(theta)) + r_val, - np.cos(theta) * np.sin(theta)])
    
    d2f_s_a = compute_d2f_s_a(S, coeff_mini, powers)
    grad_rval = np.zeros((n_theta, n_coeff))

    for i in range(n_theta):
        for j in range(n_coeff):
            grad_rval[i,j] = np.dot(d2f_s_a[i,:,j], v2[i])

    return(grad_rval)


def points(Nh, Mh, Nv, Mv):
    """
        Returns a set of points on the unity sphere and its derivative according to omega and its second derivative
        with repect to omega (Ref : Article Soare 2007)

        Input :
            - Nh : int, Number of horizontal circles
            - Mh : int, Number of points on one horizontal circle
            - Nv : int, Number of vertical circles
            - Mv : int, Number of points on one vertical circle

        Output :
            - X : ndarray of shape (Nv * Mv + Nh * Mh, 3), points
            - V : ndarray of shape (Nv * Mv + Nh * Mh, 3), points derived with respect to omega
            - A : ndarray of shape (Nv * Mv + Nh * Mh, 3), points derived twice with respect to omega
    """

    X_v = np.zeros((Nv * Mv, 3))
    X_h = np.zeros((Nh * Mh, 3))
    V_v = np.zeros((Nv * Mv, 3))
    V_h = np.zeros((Nh * Mh, 3))
    A_v = np.zeros((Nv * Mv, 3))
    A_h = np.zeros((Nh * Mh, 3))

    thetas_v = np.linspace(0, np.pi, Nv)
    omegas_v = np.linspace(0, np.pi, Mv)

    for q in range(Nv):
        for i in range(Mv):
            theta_v = thetas_v[q]
            omega_v = omegas_v[i]
            X_v[q * Mv + i] = np.array([np.cos(omega_v) * np.cos(theta_v), np.cos(omega_v) * np.sin(theta_v), np.sin(omega_v)])
            V_v[q * Mv + i] = np.array([- np.sin(omega_v) * np.cos(theta_v), - np.sin(omega_v) * np.sin(theta_v), np.cos(omega_v)])
            A_v[q * Mv + i] = np.array([- np.cos(omega_v) * np.cos(theta_v), - np.cos(omega_v) * np.sin(theta_v), - np.sin(omega_v)])
    
    thetas_h = np.linspace(0, np.pi/2, Nh)
    omegas_h = np.linspace(0, np.pi, Mh)

    for q in range(Nh):
        for i in range(Mh):
            theta_h = thetas_h[q]
            omega_h = omegas_h[i]
            X_h[q * Mh + i] = np.array([np.cos(omega_h) * np.sin(theta_h), np.sin(omega_h) * np.sin(theta_h), np.cos(theta_h)])
            V_h[q * Mh + i] = np.array([- np.sin(omega_h) * np.sin(theta_h), np.cos(omega_h) * np.sin(theta_h), np.cos(theta_h)])
            A_h[q * Mh + i] = np.array([- np.cos(omega_h) * np.sin(theta_h), - np.sin(omega_h) * np.sin(theta_h), np.cos(theta_h)])
    
    X = np.concatenate((X_v, X_h))
    V = np.concatenate((V_v, V_h))
    A = np.concatenate((A_v, A_h))

    return(X, V, A)



def constraint_pos(X, powers):
    """
        Returns the positivity constraint.

        Input :
            - X : ndarray of shape (m, 3), set of points where the yield function needs to
            be positive
            - powers, ndarray of shape (n, 3), powers of polyN2d

        Output :
            - cons : scipy.LinearConstraint

    """
    n = len(powers)
    m = len(X)

    A = np.zeros((m, n - 1 + 2)) #-1 : a0 = 1, +2 no constraint on the last 2 coeff
    B = np.zeros(m)
    for i in range(m):
        x = X[i]
        for j in range(n):
            p = powers[j]
            if j==0:
                B[i] = - np.prod(x ** p)
            else:
                A[i,j - 1] = np.prod(x ** p)

    cons = scipy.optimize.LinearConstraint(A, lb=B, ub=np.inf, keep_feasible=True)
    return(cons)

def grad_P_a(x, powers):
    """
        Returns polyN2D derived with respect to a.

        Input :
            - x : ndarray of shape (3,), point where computed
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
        
        Output :
            - grad : ndarray of shape (nmon,)
    """
    nmon = len(powers)
    grad = np.zeros(nmon)

    for i in range(nmon):
        grad[i] = np.prod(x ** powers[i])
    return(grad)

def grad_dP_a(x, powers):
    """
        Returns grad(polyN2D) derived with respect to a.

        Input :
            - x : ndarray of shape (3,), point where computed
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
        
        Output :
            - grad : ndarray of shape (3, nmon)
    """
    n_coeff = len(powers)
    coeff_grad, powers_grad = jac_polyN_2d_param(np.ones(n_coeff), powers)
    grad = np.zeros((3, n_coeff))
    for i in range(3):
        for j in range(n_coeff):
            grad[i,j] = coeff_grad[i,j] * np.prod(x ** powers_grad[i,j])
    return(grad)

def grad_d2P_a(x, powers):
    """
        Returns hess(polyN2D) derived with respect to a.

        Input :
            - x : ndarray of shape (3,), point where computed
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
        
        Output :
            - hess : ndarray of shape (3, 3, nmon)
    """
    n_coeff = len(powers)
    coeff_grad, powers_grad = jac_polyN_2d_param(np.ones(n_coeff), powers)
    coeff_hess, powers_hess = hessian_polyN_2d_param(coeff_grad, powers_grad)
    grad = np.zeros((3, 3, n_coeff))
    for i in range(3):
        for j in range(3):
            for k in range(n_coeff):
                grad[i,j,k] = coeff_hess[i,j,k] * np.prod(x ** powers_hess[i,j,k])
    return(grad)

def constraint_convex(x, v, a, powers):
    """
        Returns convexity constraint for a point x. (Ref : article Soare 2007)

        Input :
            - x : ndarray of shape (3,), point where computed
            - v : ndarray of shape (3,), point derived once with respect to omega where computed
            - a : ndarray of shape (3,), point derived twice with respect to omega where computed
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
        
        Output :
            - cons : func (ndarray of shape (nmon - 1,) -> float)
            - grad_cons : func (ndarray of shape (nmon - 1,) -> float)
    """

    def cons(coeff):
        b = np.ones(len(coeff) + 1)
        b[1:] = coeff

        n = np.sum(powers[0])
        coeff_grad, powers_grad = jac_polyN_2d_param(b[:-2], powers)
        coeff_hess, powers_hess = hessian_polyN_2d_param(coeff_grad, powers_grad)
        grad_P = grad_polyN_2d_opti(x, coeff_grad, powers_grad)[0]
        hess_P = hessian_polyN_2d_opti(x, coeff_hess, powers_hess)[0]

        P = polyN_2d_opti(x, b[:-2], powers)[0]
        dP = np.dot(grad_P, v)
        d2P = np.dot(np.matmul(hess_P, v), v) + np.dot(grad_P, a)

        res = n**2 * P**2 - (n-1) * dP ** 2 + n * P * d2P
        if res <0 :
            print(res)
        return(res)

    def grad_cons(coeff):
        b = np.ones(len(coeff) + 1)
        b[1:] = coeff
        n_coeff = len(powers)
        n = np.sum(powers[0])       
        coeff_grad, powers_grad = jac_polyN_2d_param(b[:-2], powers)
        coeff_hess, powers_hess = hessian_polyN_2d_param(coeff_grad, powers_grad)
        grad_P = grad_polyN_2d_opti(x, coeff_grad, powers_grad)[0]
        hess_P = hessian_polyN_2d_opti(x, coeff_hess, powers_hess)[0]

        P = polyN_2d_opti(x, b[:-2], powers)[0]
        dP = np.dot(grad_P, v)
        d2P = np.dot(np.matmul(hess_P, v), v) + np.dot(grad_P, a)
        
        dgrad_Pda = grad_dP_a(x, powers)
        dhess_Pda = grad_d2P_a(x, powers)

        dPda = grad_P_a(x, powers)
        ddPda = np.matmul(dgrad_Pda.T, v)
        dd2Pda = np.zeros(n_coeff)
        for i in range(n_coeff):
            dd2Pda[i] = np.dot(np.matmul(dhess_Pda[:,:,i], v), v) + np.dot(dgrad_Pda[:,i], a)
        
        grad = np.zeros(n_coeff + 2)
        grad[:n_coeff] = 2 * n * n * dPda * P - 2 * (n-1) * ddPda * dP + n * dPda * d2P + n * P * dd2Pda
        return(grad[1:])

    return cons, grad_cons

def test_cons():
    X, V, A = points(10, 10, 10, 10)
    powers = get_param_polyN_mini(4)
    coeff = np.zeros(10)
    for x, v, a in zip(X, V, A):
        cons, grad_cons = constraint_convex(x, v, a, powers)
        print(cons(coeff))
        print(grad_cons(coeff))

def remove_linear_part(df, type_test):
    """
        Remove the elastic part from a experimental test (NT6 of SH only) to only have
        elastic part.

        Input :
            - df : dataFrame, test csv
            - type_test : string in ["NT6", "SH"]

        Output :
            - df_trimmed : dataFrame, test csv with elastic part removed
    """
    if type_test == "NT6":
        col_x = "Displacement[mm]"
    elif type_test == "SH":
        col_x = "Displacement longi[mm]"
    col_y = "Force[kN]"

    window_size = 10  
    a_list = []
    for i in range(len(df) - window_size):
        model = LinearRegression()
        window = df.iloc[i:i + window_size]
        X, Y = window[col_x], window[col_y]
        X = X.values.reshape(-1, 1)
        model.fit(X, Y)
        a = model.coef_
        a_list.append(a)
    
    if type_test == "SH":
        p = 4.5
    if type_test == "NT6":
        p = 2

    a_list = np.array(a_list)
    dap_list = np.float_power(np.abs(a_list[:-1] - a_list[1:]), p)

    m = np.mean(dap_list)
    linear_end_index = np.where(dap_list > m)[0][-1]

    df_trimmed = df.iloc[linear_end_index + window_size:]

    return(df_trimmed)

def slope_study(df, type_test, material):

    if type_test == "NT6":
        col_x = "Displacement[mm]"
    elif type_test == "SH":
        col_x = "Displacement longi[mm]"
    col_y = "Force[kN]"

    # Use a rolling window to find the linear part
    window_size = 10  # Adjust the window size as needed
    a_list = []
    for i in range(len(df) - window_size):
        model = LinearRegression()
        window = df.iloc[i:i + window_size]
        X, Y = window[col_x], window[col_y]
        X = X.values.reshape(-1, 1)
        model.fit(X, Y)
        a = model.coef_
        a_list.append(a)
    
    a_list = np.array(a_list)
    m = np.mean(np.square(a_list[:-1] - a_list[1:]))

    x = df[col_x].iloc[window_size:][:-1]
    y = df[col_y].iloc[window_size:][:-1]

    fig,ax1 = plt.subplots()
    ax1.plot(x, y, color = "blue")
    ax1.set_xlabel(r"$\delta u \ [mm]$")
    ax1.set_ylabel(r"$F \ [kN]$")

    ax2 = ax1.twinx()
    ax2.plot(x, np.square(a_list[:-1] - a_list[1:]), color = "red", linewidth=0.2)
    ax2.set_ylabel(r"$\delta a^2 \ [kN/mm]$")

    
    plt.grid(1)

    plt.suptitle(f"{material} : {type_test}")
    plt.title(f"window size : {window_size}")

    plt.show()

def get_yield_stress_NT6(sigma0, material):
    """
        Returns a dict ys_ratio_nt6 containing every yield stress ratio for the available orientation of the material given.

        Input :
            - sigma0 : float, value of the yield stress for UT_0
            - material : string, the material given

        Output :
            - ys_ratio_nt6 : dict, {"ori": ys_ratio}, yield stress ratio for the given orientation
    """
    mat_tests = analyze_exp_data(material)
    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material
    ys_ratio_nt6 = {}

    for type_test in mat_tests:

        if type_test == "NT6":
            for ori in mat_tests[type_test]:

                    test = type_test + "_" + ori
                    ys_exp = np.array([])
                    for i in range(mat_tests[type_test][ori]):

                        filename = test + "_" + str(i+1) + ".csv"
                        filepath = results_exp_dir + sep + filename

                        df_NT6 = pd.read_csv(filepath)

                        thick = df_NT6.loc[0, "Thickness[mm]"]
                        innerwidth = df_NT6.loc[0, "InnerWidth[mm]"]
                        l0 = 30

                        df_NT6 = remove_linear_part(df_NT6, type_test)
                        df_NT6["PlasticStrain"] = df_NT6["AxStrain_1"] - df_NT6["AxStrain_1"].iloc[0]
                        df_NT6["PlasticStress"] = df_NT6["Force[kN]"] * 1000 / (thick * innerwidth)
                        yield_stress = df_NT6[df_NT6["PlasticStrain"] > 0.002]["PlasticStress"].iloc[0]
                        ys_exp = np.append(ys_exp,yield_stress)

                    ys = np.mean(ys_exp)
                    ys_ratio = ys / sigma0
                    ys_ratio_nt6[ori] = ys_ratio

    return(ys_ratio_nt6)

def get_yield_stress_SH(sigma0, material):
    """
        Returns a dict ys_ratio_sh containing every yield stress ratio for the available orientation of the material given.

        Input :
            - sigma0 : float, value of the yield stress for UT_0
            - material : string, the material given

        Output :
            - ys_ratio_sh : dict, {"ori": ys_ratio}, yield stress ratio for the given orientation
    """
    mat_tests = analyze_exp_data(material)
    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material
    ys_ratio_SH = {}

    for type_test in mat_tests:

        if type_test == "SH":
            for ori in mat_tests[type_test]:
                    
                    test = type_test + "_" + ori
                    ys_exp = np.array([])
                    for i in range(mat_tests[type_test][ori]):

                        filename = test + "_" + str(i+1) + ".csv"
                        filepath = results_exp_dir + sep + filename
                        df_SH = pd.read_csv(filepath)

                        thick = df_SH.loc[0, "Thickness[mm]"]
                        width = 3.15
                        l0 = 11

                        df_SH = remove_linear_part(df_SH, type_test)

                        df_SH["PlasticStress"] = df_SH["Force[kN]"] * 1000 / (thick  * width)
                        df_SH["PlasticStrain"] = (df_SH["Displacement longi[mm]"] - df_SH["Displacement longi[mm]"].iloc[0]) / l0
                        yield_stress = df_SH[df_SH["PlasticStrain"] > 0.002]["PlasticStress"].iloc[0]
                        ys_exp = np.append(ys_exp,yield_stress)

                    ys = np.mean(ys_exp)
                    ys_ratio = ys / sigma0
                    ys_ratio_SH[ori] = ys_ratio

    return(ys_ratio_SH)

def get_dir_pst(coeff, powers):
    """
        For given coefficients of the yield function minimalistic, returns the direction of the PST
        points for orientation [0°, 45°, 90°]
        
        Input :
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function minimalistic
            - powers : ndarray of shape (nmon, 3), powers[i,j] = powers of the variable j in the monomial i

        Output :
            - dir_pst : dict, {"ori": (phi, theta)} the angles being defined for spherical coordonnates on wikipedia,
            the direction of the pst points
    """
    def gradg(S):
        return(grad_f_min_squared(S, coeff, powers))
    
    n_phi = 100
    phis = np.linspace(0, np.pi/2, n_phi)
    S = np.zeros((n_phi, 6))

    for i in range(n_phi):
        phi = phis[i]
        S[i] = np.array([np.cos(phi), np.sin(phi), 0, 0, 0, 0]) 

    grad00 = np.abs(gradg(S)[:,1])
    grad90 = np.abs(gradg(S)[:,0])

    phi_00 = phis[np.argmin(grad00)]
    theta_00 = np.pi/2
    phi_90 = phis[np.argmin(grad90)]
    theta_90 = np.pi/2

    #PART 45 DEGREES
    n_theta = 100
    thetas = np.linspace(0, np.pi/2)
    score = np.zeros((n_phi, n_theta))
    v = np.array([1/2, 1/2, 1])
    
    def f(phi, theta):
        u = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), 0, np.cos(theta), 0, 0]) 
        gradf = gradg(u)[0,[0, 1, 3]]
        w = np.cross(gradf, v)
        norm = np.linalg.norm(w)
        return(norm)
    
    X, Y = np.meshgrid(phis, thetas, indexing='ij')

    score = np.vectorize(f)(X, Y)
    index_colinear = np.argmin(score)
    min_index_2d = np.unravel_index(index_colinear, score.shape)
    phi_45 = phis[min_index_2d[0]]
    theta_45 = thetas[min_index_2d[1]]

    dir_pst = {"00" : (phi_00, theta_00), "45": (phi_45, theta_45), "90" : (phi_90, theta_90)}
    return(dir_pst)

def add_nt6(old_df, ys_ratio_nt6, dir_pst):
    """
        Add PST data to the old_df dataFrame used to calibrate the yield function.

        Input :
            - old_df : dataFrame, must contain the columns ["s11", "s22", "s33", "s12", "s13", "s23"]
            - ys_ratio_nt6 : dict, {"ori": ys_ratio}, ys_ratio for the given orientation
            - dir_pst : dict, {"ori": (phi, theta)} the angles being defined for spherical coordonnates on wikipedia,
            the direction of the pst points

        Output :
            - new_df : dataFrame, old_df with added rows corresponding to the PST points
    """
    new_df = old_df.copy(deep=True)

    for ori in ys_ratio_nt6:
        phi, theta = dir_pst[ori]
        ys = ys_ratio_nt6[ori]

        psi = float(ori) / 360 * 2 * np.pi
        dep = np.array([np.square(np.cos(psi)), np.square(np.sin(psi)), np.sin(2 * psi)])
        #only sigma axial
        ds = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])

        r = ys / np.dot(ds, dep)
        s = r * ds

        x = s[0]
        y = s[1]
        z = s[2]
        row = {"s11" : x, "s22":y , "s33":0, "s12":z, "s13":0, "s23":0, "Type":"e2", "Rval":0.0}

        for key in new_df:
            if not(key in row):
                row[key] = np.nan

        new_df.loc[len(new_df)] = row
    return(new_df)

def add_sh(old_df, ys_ratio_sh):
    """
        Add sh data to the old_df dataFrame used to calibrate the yield function.

        Input :
            - old_df : dataFrame, must contain the columns ["s11", "s22", "s33", "s12", "s13", "s23"]
            - ys_ratio_sh : dict, {"ori": ys_ratio}, ys_ratio for the given orientation

        Output :
            - new_df : dataFrame, old_df with added rows corresponding to the shear points
    """

    new_df = old_df.copy(deep=True)
    for ori in ys_ratio_sh:
        ys = ys_ratio_sh[ori]

        theta = float(ori[1:]) / 360 * 2 * np.pi
        if ori[0] == "m":
            theta = - theta

        #only shear comp
        ds = np.array([-np.sin(2 * theta), np.sin(2 * theta), np.cos(2 * theta)])
        s = ys * ds

        x = s[0]
        y = s[1]
        z = s[2]

        row = {"s11" : x, "s22":y , "s33":0, "s12":z, "s13":0, "s23":0, "Type":"e2", "Rval":0.0 }
        for key in new_df:
            if not(key in row):
                row[key] = np.nan
        new_df.loc[len(new_df)] = row
    return(new_df)



def leading_principal_minors(lmatrix):
    """
        Returns the leading principal minors of a matrix or a list of matrixes
        Input :
            - lmatrix : ndarray of shape (k, n, n) or (n, n)
        Output :
            - lead_minors : ndarray of shape (k, n)
    """
    if lmatrix.ndim == 2 :
        lmatrix = np.expand_dims(lmatrix, axis = 0)
    k = lmatrix.shape[0]
    n = lmatrix.shape[1]

    lead_minors = np.zeros((k, n))
    
    for i in range(n):
        lead_minors[:, i] = np.linalg.det(lmatrix[:, :i, :i])

    return lead_minors

def check_convexity(coeff, powers):
    """
        Returns True if a polyN minimalistic yield function is convex or not

        Input :
            - coeff : ndarray of shape (nmon + 2,), coefficients or the yield function
            - powers : ndarray of shape (nmon, 3), powers of polyN2D
    """
    coeff_mini = coeff[:-2]
    k1 = coeff[-2]
    k2 = coeff[-1]

    Nh = 100
    Mh = 30
    Nv = 100
    Mv = 30
    eps = 10e-7
    X,_,_ = points(Nh, Mh, Nv, Mv)

    coeff_grad, powers_grad = jac_polyN_2d_param(coeff_mini, powers)
    coeff_hess, powers_hess = hessian_polyN_2d_param(coeff_grad, powers_grad)

    hess = hessian_polyN_2d_opti(X, coeff_hess, powers_hess)
    L = leading_principal_minors(hess)
    m = min(L.flatten())

    if m > - eps and k1 > 0 and k2 > 0:
        return(True)
    
    return(False)


def optiCoeff_polyN_mini(df, degree, weight_ut, weight_e2, weight_vir, init_guess):
    """
        Returns the optimized coefficients of the yield surface minimaslistic on experimental data (UTs) and virtual data from a protomodel.

        Input :
            - df : pd.Dataframe, must contain columns ["d11", "d22", "s12", "s13", "s23"] of yield stresses points. And if available ["Rval"] with ["LoadAngle"].
            - degree : integer, degree of polyN
            - weight_ut : float, weight for experimental stress data
            - weight_e2 : float, proportion exp_e2 / exp_ut
            - weight_vir : float, proportion exp_vir / exp_ut
            - init_guess : ndarray of shape (nmon,), initial guess for optimization

        Output :
            - coeff : ndarray of shape (nmon,), coeff of yield surface model
    """

    data = df[["s11", "s22", "s33", "s12", "s13", "s23"]].values

    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X_stress = polyN.fit_transform(t_linear(data)[:,[0,1,2]])

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
    n_data = len(data)

    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    ndata_rval = np.count_nonzero(index_rval)
    r_vals = df["Rval"].iloc[index_rval].values
    thetas = df["LoadAngle"].iloc[index_rval].values

    index_s = np.where(df["Type"] == "e", True ,False)
    index_ns = np.where(df["Type"] == "e2", True ,False)
    index_v = ~index_s & ~index_ns

    n_data_ut = np.sum(index_s)
    n_data_ns = np.sum(index_ns)
    n_data_v = n_data - n_data_ut - n_data_ns
    weight_s = np.zeros(n_data)

    w = weight_ut / (n_data_ut + weight_e2 * n_data_ns + weight_vir * n_data_v)

    w_ut = w
    w_e2 = w * weight_e2
    w_v = w * weight_vir

    weight_s[index_s] = w_ut
    weight_s[index_ns] = w_e2
    weight_s[index_v] = w_v

    weight_r = np.ones(ndata_rval)
    weight_r = (1 - weight_ut) *  weight_r / np.sum(weight_r)

    print(weight_s, np.sum(weight_s), np.sum(weight_r))

    def J(a):
        b = np.ones(len(a) + 1)
        b[1:] = a

        g = f_min_squared(data, b, powers)
        f = np.float_power(g, 0.5)

        eq_ys = f - 1
        eq_rval = eq_rval_ut(thetas, r_vals, b[:-2], powers)

        c = np.sum(weight_s * np.square(eq_ys))
        d = np.sum(weight_r * np.square(eq_rval))

        res = 0.5 * (c + d)
        print(res)
        return(res)

    def grad_J(a):
        b = np.ones(len(a) + 1)
        b[1:] = a
        n_coeff = len(a)
        grad = np.zeros(n_coeff)

        p = polyN_2d(data, b[:-2], powers)
        g = f_min_squared(data, b, powers)

        f = np.float_power(g, 0.5)
        eq_ys = f - 1
        eq_rval = eq_rval_ut(thetas, r_vals, b[:-2], powers)
        grad_rval = grad_rval_ut(thetas, r_vals, b[:-2], powers)

        for i in range(len(a) - 2):
            grad_ys = 0.5 * (2/degree) * X_stress[:, i+1] * np.float_power(p,(2/degree - 1)) / f
            grad_c = np.sum(weight_s * eq_ys * grad_ys)
            grad_d = np.sum(weight_r * eq_rval * grad_rval[:, i+1])
            grad[i] = grad_c + grad_d

        grad_ys = 0.5 * np.square(data[:,4]) / f
        grad[-2] = np.sum(weight_s * eq_ys * grad_ys)

        grad_ys = 0.5 * np.square(data[:,5]) / f
        grad[-1] = np.sum(weight_s * eq_ys * grad_ys)

        return(grad)
    
    constraints = []

    X, V, A = points(10, 15, 10, 15)
    cons_pos = constraint_pos(X, powers)
    constraints.append(cons_pos)

    for x,v,a in zip(X, V, A):
        cons, grad_cons = constraint_convex(x, v, a, powers)
        cons_convex = scipy.optimize.NonlinearConstraint(cons, lb=0, ub=np.inf, jac=grad_cons, keep_feasible=True)
        constraints.append(cons_convex)

    print(f"Number of coefficients for degree {degree}:", nmon + 2)

    b = np.zeros(len(init_guess))
    b[-2] = 3
    b[-1] = 3

    res = 0
    for i in range(degree // 2, -1, -1):
        n = 2 * i + 1
        val = np.max(np.abs(init_guess[res:res + n]))
        b[res:res + n] = 5 * val
        res = res + n

    a0 = init_guess[1:]
    b = b[1:]
    I = np.eye(len(a0))
    cons_hypercube = scipy.optimize.LinearConstraint(I, lb = a0 - b, ub= a0 + b, keep_feasible=True)
    #constraints.append(cons_hypercube)

    options = {"verbose" : 3, "maxiter" : 1000}

    opt = scipy.optimize.minimize(J, x0=a0, jac=grad_J, constraints=constraints, tol=10e-25, options=options)

    while (not opt.success and opt.nit == options['maxiter']):
        new_a0 = a0 + 0.1 * np.random.uniform(low= - np.ones(len(a0)), high=np.ones(len(a0)))
        try :
            opt = scipy.optimize.minimize(J, x0=new_a0, jac=grad_J, method="trust-constr", constraints=constraints, options=options)
        except ValueError:
            pass

    coeff_temp = opt.x
    coeff = np.ones(len(coeff_temp) + 1)
    coeff[1:] = coeff_temp

    return(coeff)

def optiCoeff_pflow_mini(law, coeff, material, powers):
    """
        Returns the a, b, c coefficients of the hardening law defined by the user and the young modulus of the material.

        Input :
            - law : string, law in ["swift", "voce"]
            - coeff : ndarray of shape (nmon + 2,), yield surface minimalistic coefficients
            - material : string
            - powers : ndarray of shape (nmon, 3), powers[i,j] power of the variable j in the monomial j

        Output :
            - coeff_law : ndarray of shape (3,) or (7,), hardening law coefficients
            - ymod : float, Young modulus
    """
    
    def f(S):
        return np.power(f_min_squared(S, coeff, powers), 1/2)

    foldername = polyN_dir + sep + "calibration_data" + sep + material
    filename_out = f"data_plasticlaw_{material}.csv"
    filepath = foldername + sep + filename_out

    df_law = pd.read_csv(filepath)

    ymod = df_law.iloc[0, df_law.columns.get_loc("YoungModulus")]
    ep = df_law["PlasticStrain"].values

    sigma = df_law["PlasticStress"].values
    S = np.expand_dims(sigma, axis = 1)
    S = np.concatenate((S, np.zeros((sigma.shape[0],5))), axis = 1)

    sig0 = f(S[0])[0]
    if law == "swift" :

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
        return(np.array([a, b, c]), ymod)
    
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

        print(a,b,c)
        return(np.array([a, b, c]), ymod)
    
    if law == "swiftvoce":

        print("Looking for swiftvoce coeff in '_YLD2000_pre.csv'")
        filepath = polyN_dir + sep + "results_exp" + sep + material + sep + "_YLD2000_pre.csv"
        try :
            df_law = pd.read_csv(filepath)
            coeff_law = df_law.loc[5].values[:7]
            return(coeff_law, ymod)
        except :
            print("No file named '_YLD2000_pre.csv' containing the coefficients of swiftvoce")
    

def write_coeff_abq_mini(coeff, coeff_law, ymod, enu, protomodel, degree, material, law, density, powers, p=0, m=0):
    """
        Write the user material file in running directory.

        Input :
            - coeff : ndarray of shape (nmon + 2,), yield surface minimalistic coefficients
            - coeff_law :
            - ymod : float, Young modulus
            - enu : float, Poisson ratio
            - protomodel : string
            - degree : int
            - material : string
            - law : string
            - density : float
            - powers : ndarray of shape (nmon, 3), powers[i,j] power of the variable j in the monomial j
            - p : list
            - m : int, to save file
    
    """
    nmon = len(powers)

    if str(p)=="0" and m==0:
        filename = "{}_abq_deg{}mini_{}_{}.inp".format(material, degree, law, protomodel)
    else:
        filename = "{}_abq_deg{}mini_{}_{}_{}_{}.inp".format(material, degree, law, protomodel, p, m)

    foldername = polyN_dir + sep + "running"
    filepath = foldername + sep + filename

    if law == "swiftvoce":
        a0 = 11
        a = coeff_law[0]
        b = coeff_law[1]
        c = coeff_law[2]
        d = coeff_law[3]
        e = coeff_law[4]
        f = coeff_law[5]
        alpha = coeff_law[6]
    else:
        a0 = 7
        a = coeff_law[0]
        b = coeff_law[1]
        c = coeff_law[2]

    with open(filepath, "w") as file:
        file.write("*USER MATERIAL, constants={}\n".format(a0 + 2 + nmon))

        if law == "swiftvoce":
            file.write("{}, {}, {}, {}, {}, {}, {}, {}, \n".format(ymod, enu, a, b, c, d, e, f))
            file.write("{}, {}, {}, ".format(alpha, degree, nmon + 2))
        else:
            file.write("{}, {}, {}, {}, {}, {}, {}, ".format(ymod, enu, a, b, c, degree, nmon + 2))

        n0 = 0
        while n0 < nmon + 2:
            for k in range(0, degree + 1):
                for j in range(0, degree + 1 - k):
                    i = degree - k - j
                    if n0 < nmon:
                        i0, j0, k0 = powers[n0]
                        if (i==i0) and (j==j0) and (k==k0):
                            file.write("{}, ".format(coeff[n0]))
                            n0 = n0 + 1
                            if (n0 + a0) % 8 == 0:
                                file.write("\n")
                    elif n0 < nmon + 2:
                        file.write("{}, ".format(coeff[n0]))
                        n0 = n0 + 1
                        if (n0 + a0) % 8 == 0:
                            file.write("\n")
                            
        file.write("\n")
        file.write("*DENSITY\n")
        file.write("{}".format(density))

"""-----------------------------------------------TESTS--------------------------------------------------------"""

def test_gradf1(degree):
    powers = get_param_polyN_mini(degree)
    ncoeff = len(powers)
    coeff = np.ones(ncoeff)
    S = np.array([[1,0,0,0,0,0], [1,0,0,0,0,0]])
    print(gradf1(S, coeff, powers))

def test_hessf1(degree):
    powers = get_param_polyN_mini(degree)
    ncoeff = len(powers)
    coeff = np.ones(ncoeff)
    S = np.array([[1,0,0,0,0,0], [1,0,0,0,0,0]])
    print(hessf1(S, coeff, powers))

def test_gradf2(degree):
    powers = get_param_polyN_mini(degree)
    ncoeff = len(powers)
    coeff = 7 * np.ones(2)
    S = np.array([[1,0,0,0,1,1],[1,0,0,0,1,1]] )
    print(gradf2(S, coeff))

def test_hessf2(degree):
    powers = get_param_polyN_mini(degree)
    ncoeff = len(powers)
    coeff = 7 * np.ones(2)
    S = np.array([[1,0,0,0,1,1],[1,0,0,0,1,1]] )
    print(hessf2(S, coeff))

def test_sh():
    p = read_param()
    material = p["material"]
    protomodel = p["protomodel"]
    
    get_hardening_law(material)
    df = readData_2d(material, protomodel)
    sigma0 = df["YieldStress"].iloc[0]
    
    print(get_yield_stress_SH(sigma0, material))

def test_remove():
    p = read_param()
    material = p["material"]

    mat_tests = analyze_exp_data(material)
    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material

    for type_test in mat_tests:

        if type_test == "SH":
            for ori in mat_tests[type_test]:
                    
                    test = type_test + "_" + ori
                    ys_exp = np.array([])
                    for i in range(mat_tests[type_test][ori]):

                        filename = test + "_" + str(i+1) + ".csv"
                        filepath = results_exp_dir + sep + filename
                        df_SH = pd.read_csv(filepath)
                        plt.plot(df_SH["Displacement longi[mm]"], df_SH["Force[kN]"])
                        df = remove_linear_part(df_SH, type_test)
                        plt.plot(df["Displacement longi[mm]"], df["Force[kN]"])
                        plt.show()

        if type_test == "NT6":
            for ori in mat_tests[type_test]:

                    test = type_test + "_" + ori
                    ys_exp = np.array([])
                    for i in range(mat_tests[type_test][ori]):

                        filename = test + "_" + str(i+1) + ".csv"
                        filepath = results_exp_dir + sep + filename

                        df_NT6 = pd.read_csv(filepath)
                        plt.plot(df_NT6["Displacement[mm]"], df_NT6["Force[kN]"])
                        df = remove_linear_part(df_NT6, type_test)
                        plt.plot(df["Displacement[mm]"], df["Force[kN]"])
                        plt.show()

def test_slope_study():
    p = read_param()
    material = p["material"]

    mat_tests = analyze_exp_data(material)
    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material

    for type_test in mat_tests:

        if type_test == "SH":
            for ori in mat_tests[type_test]:
                    
                    test = type_test + "_" + ori
                    ys_exp = np.array([])
                    for i in range(mat_tests[type_test][ori]):

                        filename = test + "_" + str(i+1) + ".csv"
                        filepath = results_exp_dir + sep + filename
                        df_SH = pd.read_csv(filepath)
                        slope_study(df_SH, type_test, material)

        if type_test == "NT6":
            for ori in mat_tests[type_test]:

                    test = type_test + "_" + ori
                    ys_exp = np.array([])
                    for i in range(mat_tests[type_test][ori]):

                        filename = test + "_" + str(i+1) + ".csv"
                        filepath = results_exp_dir + sep + filename

                        df_NT6 = pd.read_csv(filepath)
                        slope_study(df_NT6, type_test, material)


def test_points(Nh, Mh, Nv, Mv):
    X, V, A = points(Nh, Mh, Nv, Mv)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x, y, z = X[:,0], X[:,1], X[:,2]
    # Plot the points
    ax.scatter(x, y, z, c='b', marker='o')

    # Set labels
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show() 

def firstopti_mini():
    p = read_param()

    material = p["material"]
    gseed = int(p["gseed"])
    enu = float(p["enu"])
    density = float(p["density"])
    nb_virtual_pt = int(p["nb_virtual_pt"])
    degree = int(p["degree"])
    sh = int(p["sh"])
    nt6 = int(p["nt6"])
    weight_ut = float(p["weight_ut"])
    weight_vir = float(p["weight_vir"])
    weight_e2 = float(p["weight_e2"])
    protomodel = p["protomodel"]
    law = p["law"]

    gen_v_data = int(p["gen_v_data"])
    gen_e_data = int(p["gen_e_data"])

    opti = int(p["opti"])
    loadcoeff = int(p["loadcoeff"])

    export_coeff_abq = int(p["export_coeff_abq"])
    savecoeff = int(p["savecoeff"])

    print(material)

    if gen_v_data :
        print("Generation of virtual data")
        export_virtual_data(protomodel, material, nb_virtual_pt)
        print("Generation ended")

    if gen_e_data:
        print("Processing experimental data")
        export_exp_data(material)
        print("Processing ended")

    get_hardening_law(material)
    df = readData_2d(material, protomodel)
    sigma0 = df["YieldStress"].iloc[0]

    if sh:
        ys_ratio_sh = get_yield_stress_SH(sigma0, material)
        df = add_sh(df, ys_ratio_sh)

    powers = get_param_polyN_mini(degree)
    nmon = len(powers)


    if degree==2:
        C = np.array([1, -1, 1, 3, 3, 3])
    if degree==4:
        C = np.array([1, -2, 3, -2, 1, 6, -6, 6, 9, 3, 3])
    if degree==6:
        C = np.array([1, -3, 6, -7, 6, -3, 1, 12, -24, 30, -24, 12, 14, -2, 14, 30, 3, 3])
    if degree==8:
        C = np.array([1, -4, 14, -30, 40, -30, 14, -4, 1, 11, -21, 17, -1, 17, -21, 11, 40, -4, 18, -4, 40, 42, -6, 42, 120, 3, 3])

    if opti:
        new_df = df.copy(deep=True)

        if nt6:
            ys_ratio_nt6 = get_yield_stress_NT6(sigma0, material)
            for i in range(3):
                coeff = optiCoeff_polyN_mini(new_df, degree,  weight_ut, weight_e2, weight_vir, C)
                dir_pst = get_dir_pst(coeff, powers)
                new_df = add_nt6(df, ys_ratio_nt6, dir_pst)
        else : 
            coeff = optiCoeff_polyN_mini(new_df, degree,  weight_ut, weight_e2, weight_vir, C)

    elif loadcoeff:
        coeff = np.load(polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy")
    else :
        coeff = C

    coeff_law, ymod = optiCoeff_pflow_mini(law, coeff, material, powers)

    if savecoeff:
        np.save(polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy", coeff)
        np.save(polyN_dir + sep + material + f"_{law}_mini_coeff.npy", np.append(coeff_law, ymod))

    if export_coeff_abq:
        print(nmon, coeff)
        write_coeff_abq_mini(coeff, coeff_law, ymod, enu, protomodel, degree, material, law, density, powers, p=0, m=0)

if __name__ == "__main__":
    firstopti_mini() 
    pass