
import pandas as pd
import sklearn
import sklearn.preprocessing
import matplotlib.pyplot as plt
import numpy as np
import time
import multiprocessing
import os
import sys

file_dir = os.path.dirname(os.path.abspath(__file__))
polyN_dir = os.path.dirname(file_dir)
sep = os.sep
sys.path.append(polyN_dir)

from read import read_param

p = read_param()

degree = int(p["degree"])
material = p["material"]
nb_dir = 1

tol_kg = 1e-7
tol_yf = 0.001
tol_minor = 1e-9

kg_max = tol_kg
kg_min = -5e-4
minor_min = - tol_minor
minor_max = tol_minor

np.random.seed(6)

if degree==2:
    M = np.array([3, 3, 3, 3, 3, 3])
if degree==4:
    M = np.array([9, 18, 27, 18, 9, 18, 18, 18, 9, 18, 18, 18, 18, 9, 0, 0, 18, 18, 18, 18, 18, 9])
if degree==6:
    M = np.array([27, 81, 162, 189, 162, 81, 27, 81, 162, 243, 162, 81, 81, 81, 81, 27, 81, 162, 243, 162, 81, 162, 162, 162, 81, 81, 81, 81, 81, 27, 0, 0, 0, 0, 81, 162, 243, 162, 81, 162, 162, 162, 81, 162, 162, 162, 162, 81, 81, 81, 81, 81, 81, 27,])

nb_coeff = len(M)

def get_param_polyN(degree):
    """
        Returns the parameters of polyN according to the degree
        Input :
            - degree : integer, degree of the polyN function
        Output :
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
    nmon = len(powers_hessian[0][0])

    for i in range(5):
        for j in range(5):
            for k in range(nmon):
                p = powers_hessian[i][j][k]
                jac_grad_f[:,i,j] = jac_grad_f[:,i,j] + coeff_hessian[i][j][k] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev()
    hessian_polyN = np.transpose(np.dot(jac_d.T,np.dot(jac_grad_f[:], jac_d)), (1, 2, 0))
    return(hessian_polyN)

"""-----------------------------------------------------------------------------------------"""
def plot_implicit_coeff(coeff, powers, bbox=(-1.5,1.5)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''

    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig, ax = plt.subplots(nrows = 1, ncols = 1, subplot_kw={"projection":"3d"})
    A = np.linspace(xmin, xmax, 50) # resolution of the contour
    B = np.linspace(xmin, xmax, 20) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    def f(S):
        return polyN(S, coeff, powers)

    for i in range(1):

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
            ax.contour(X, Y, Z+z, [z], zdir='z')
            # [z] defines the only level to plot for this contour for this value of z
        print("Fin du plot sur XY")
        print("Début du plot sur XZ")
        for y in B: # plot contours in the XZ plane
            X,Z = A1,A2
            Y = f_plane(X,y,Z) - 1
            ax.contour(X, Y+y, Z, [y], zdir='y')

        print("Fin du plot sur XZ")
        print("Début du plot sur YZ")
        for x in B: # plot contours in the YZ plane
            Y,Z = A1,A2
            X = f_plane(x,Y,Z) - 1
            ax.contour(X+x, Y, Z, [x], zdir='x')
        print("Fin du plot sur YZ")
        # must set plot limits because the contour will likely extend
        # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
        # to encompass all values in the contour.
        ax.set_zlim3d(zmin,zmax)
        ax.set_xlim3d(xmin,xmax)
        ax.set_ylim3d(ymin,ymax)

   
    plt.show()

def generate_dir(n, dim):
    """
        Returns n random directions on the unit sphere of dimension dim
        Input :
            - n : integer, number of directions
            - dim : integer, number of dimensions
        Output :
            - dirs : ndarray of shape (n, dim)
    """
    dirs = np.random.normal(0, 1, (n, dim))
    norms = np.linalg.norm(dirs, axis=1)
    dirs =  0.1 * (dirs.T/norms).T
    return(dirs)

def specific_dir(n, dim, theta, eps):
    u = 0.1 * np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
    dirs = np.random.normal(0, 1, (n, dim))
    norms = np.linalg.norm(dirs, axis=1)
    dirs = u + eps * (dirs.T/norms).T
    return(dirs)

def cofactor_matrix(lmatrix):
    """
        Returns the cofactor matrix of a matrix or list of matrixes
        Input :
            - lmatrix : ndarray of shape (k, n, n) or (n, n)
        Output :
            - cofactors : ndarray of shape (k, n, n)
    """
    if lmatrix.ndim == 2 :
        lmatrix = np.expand_dims(lmatrix, axis = 0)
    k = lmatrix.shape[0]
    n = lmatrix.shape[1]

    cofactors = np.zeros((k, n, n))
    
    for i in range(n):
        for j in range(n):
            minor = np.delete(np.delete(lmatrix, i, axis=1), j, axis=2)
            cofactors[:, i, j] = (-1) ** (i + j) * np.linalg.det(minor)

    return cofactors

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

def data_yf_sphere(f, itermax, nb_pt):
    """
        Returns random points such as f(x) = 1.
        Input :
            - f : (6,) or (n, 6) -> float
            - itermax : integer, maximum iterations of the error ad trial method
            - nb_pt : integer, number of points
        Output :
            - data : ndarray of shape (nb_pt, 6)
    """
    data = np.zeros((nb_pt, 6))
    j = 0

    while j < nb_pt:

        E = np.zeros(6)
        m = -1
        u = generate_dir(1, 6)[0]
        it = 0
        lamb = 1

        #Find a point such as f(x) > 1
        while m < tol_yf and it < itermax :

            E = E + lamb * u
            m = f(E) - 1
            lamb = lamb * 2
            it = it + 1

        # If f(x) = 1
        if it != itermax :
            S = np.zeros(6)
            res = f(E) - 1
            it = 0

            #I suppose it always ends cause of the nature of f (here always polyN)
            while abs(res) > tol_yf:
                I = (S+E)/2
                res = f(I) - 1
                if res > 0 :
                    E = I
                else:
                    S = I
                it = it + 1
            
            data[j] = I
            j = j + 1
    
    return(data)

def data_yf_sphere_spec(f, itermax, n_pt_angle):
    """
        Returns random points such as f(x) = 1.
        Input :
            - f : (6,) or (n, 6) -> float
            - itermax : integer, maximum iterations of the error ad trial method
            - nb_pt : integer, number of points
        Output :
            - data : ndarray of shape (nb_pt, 6)
    """
    data = np.zeros((n_pt_angle * 6, 6))
    thetas = np.linspace(0, np.pi/2, 6)
    eps = 0.01
    for i in range (len(thetas)):
        j = 0
        theta = thetas[i]
        while j < n_pt_angle:

            E = np.zeros(6)
            m = -1
            u = specific_dir(1, 6, theta, eps)[0]
            it = 0
            lamb = 1

            #Find a point such as f(x) > 1
            while m < tol_yf and it < itermax :

                E = E + lamb * u
                m = f(E) - 1
                lamb = lamb * 2
                it = it + 1

            # If f(x) = 1
            if it != itermax :
                S = np.zeros(6)
                res = f(E) - 1
                it = 0

                #I suppose it always ends cause of the nature of f (here always polyN)
                while abs(res) > tol_yf:
                    I = (S+E)/2
                    res = f(I) - 1
                    if res > 0 :
                        E = I
                    else:
                        S = I
                    it = it + 1
                
                data[i * n_pt_angle + j] = I
                j = j + 1
    
    return(data)

def check_convexity_lpm(coeff, powers, itermax, nb_pt):
    """
        Returns the minimum leading principal minor among nb_pt points of the polyN definded by the coefficients coeff
        Input :
            - coeff : ndarray of shape (22,)
            - itermax : integer, maximum iterations of the error and trial method for data_yf_sphere
            - nb_pt : integer, number of points
    """
    coeff_grad, powers_grad = jac_polyN_param(coeff, powers)
    coeff_hessian, powers_hessian = hessian_polyN_param(coeff_grad, powers_grad)

    def f(S):
        return(polyN(S, coeff, powers))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))

    print("Checking convexity on the whole surface")
    data = data_yf_sphere(f, itermax, nb_pt)
    hessian_data = hessian_f(data)
    L = leading_principal_minors(hessian_data)
    m = min(L.flatten())
    return(m)

def check_convexity_lpm_spec(coeff, powers, itermax, n_pt_angle):
    """
        Returns the minimum leading principal minor among nb_pt points of the polyN definded by the coefficients coeff
        Input :
            - coeff : ndarray of shape (22,)
            - itermax : integer, maximum iterations of the error and trial method for data_yf_sphere
            - nb_pt : integer, number of points
    """
    coeff_grad, powers_grad = jac_polyN_param(coeff, powers)
    coeff_hessian, powers_hessian = hessian_polyN_param(coeff_grad, powers_grad)

    def f(S):
        return(polyN(S, coeff, powers))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))

    print("Checking convexity on specified angles")
    data = data_yf_sphere_spec(f, itermax, n_pt_angle)
    hessian_data = hessian_f(data)
    L = leading_principal_minors(hessian_data)
    m = min(L.flatten())
    return(m)

coeff = np.load(polyN_dir + sep + material + "_polyN_coeff.npy")
powers = get_param_polyN(degree)

coeff[0] = coeff[0]
lpm_spec = check_convexity_lpm_spec(coeff, powers, 100, 100)
lpm = check_convexity_lpm(coeff, powers, 100, 1000)

print(lpm_spec, lpm)
if lpm_spec < - tol_minor or lpm < - tol_minor :
    plot_implicit_coeff(coeff, powers)