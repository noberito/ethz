import pandas as pd
import sklearn
import sklearn.preprocessing
import matplotlib.pyplot as plt
import numpy as np
import time
import multiprocessing
import os

degree = 6
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



"""-----------------------------------PARAMETERS FOR POLYN------------------------------------------"""
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
    nmon_abq = len(powers)

    selected_indices = []
    for i in range(nmon_abq):
        k,l,m,n,p = powers[i]
        if ((m==n) and (n==p)) or ((m%2==0) and (n%2==0) and (p%2==0)):
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
    powers = powers.T[sorted_indices]
    X = X[:,sorted_indices]

    ndata, nmon = X.shape

    return(powers, nmon, nmon_abq)

powers, nmon, nmon_abq = get_param_polyN(degree)

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

"""----------------------------------------------GRADIENT AND HESSIAN OF POLYN DEFINITION-----------------------------------------------------------"""

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
            - coeff_grad (float ndarray of shape (5, nmon)) : coeff_grad[i] contains the 
            coefficients of each monomial derived with respect to dev[i]
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

    grad = np.zeros(6)
    grad_f = np.zeros((len(X),5))

    for i in range(5):
        for j in range(nmon):
            p = powers_grad[i][j]
            grad_f[:,i] = grad_f[:,i] + coeff_grad[i][j] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev()

    grad= np.dot(grad_f, jac_d)
    return(grad)

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

    hessian = np.zeros((6,6))
    jac_grad_f = np.zeros((len(X), 5, 5))

    for i in range(5):
        for j in range(5):
            for k in range(nmon):
                p = powers_hessian[i][j][k]
                jac_grad_f[:,i,j] = jac_grad_f[:,i,j] + coeff_hessian[i][j][k] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev()
    hessian = np.transpose(np.dot(jac_d.T,np.dot(jac_grad_f[:], jac_d)), (1, 2, 0))
    return(hessian)
"""-----------------------------------------------------------------------------------------"""

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
    dirs = (dirs.T/norms).T
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


def gauss_curv(f, grad_f, hessian_f, nb_pt_check):
    """
        Returns the gauss curvature. Not used cause of bad performance
    """
    dirs = generate_dir(nb_pt_check, 6)
    itermax = 100
    data = np.zeros((nb_pt_check, 6))
    j = 0

    for i in range(nb_pt_check):
        E = np.zeros(6)
        m = -1
        u = dirs[i]
        it = 0
        lamb = 1
        while m < tol_yf and it < itermax :
            E = E + lamb * u
            m = f(E) - 1
            lamb = lamb * 1.01
            it = it + 1

        if it != itermax :
            
            S = np.zeros(6)
            res = f(E) - 1
            it = 0
            while abs(res) > tol_yf and it < itermax:
                I = (S+E)/2
                res = f(I) - 1
                if res > 0 :
                    E = I
                else:
                    S = I
                it = it + 1
            data[j] = I
            j = j + 1
    
    data = data[:j]
    
    K = - np.ones(j)
    grad_data = grad_f(data)
    norm_grad_data = np.linalg.norm(grad_data, axis = 1)
    hessian_data = hessian_f(data)

    cofactor_hessian_data = cofactor_matrix(hessian_data)
    L = leading_principal_minors(hessian_data)
    pos_L = np.empty(j)

    for i in range(j):
        #print(i, grad_data[i], norm_grad_data[i], cofactor_hessian_data[i])
        K[i] = K[i] * np.dot(grad_data[i], np.dot(grad_data[i], cofactor_hessian_data[i])) / (norm_grad_data[i] ** 7)
        
        if np.count_nonzero(L[i] > - tol_minor) < 6:
            pos_L[i] = True
        else :
            pos_L[i] = False
    
    bad_points_K = np.count_nonzero(K > tol_kg)
    bad_points_L = np.count_nonzero(pos_L)
    print("The proportion of convex spots with K is : {}".format(1 - bad_points_K/len(K)))
    print("The proportion of convex spots with L is : {}".format(1 - bad_points_L/len(L)))

    return(K)

def check_convexity_kg(coeff, nb_pt_check):
    """
        Checking of gauss curvature. Not used cause of bad perf
    """
    coeff_grad, powers_grad = jac_polyN_param(coeff, powers)
    coeff_hessian, powers_hessian = hessian_polyN_param(coeff_grad, powers_grad)

    def f(S):
        return(polyN(S, coeff))
                
    def grad_f(S):
        return(grad_polyN(S, coeff_grad, powers_grad))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))

    K = gauss_curv(f, grad_f, hessian_f, nb_pt_check)
    M, m = np.max(K), np.min(K)
    moy = np.mean(K)
    n_pos = np.count_nonzero(K > 0)
    n_neg = np.count_nonzero(K < 0)


    return(M, m, moy, n_pos, n_neg)

def gauss_curv_coeff(coeff, nb_pt_check, val_test=np.array([2, 1, 0, 1, 1, 1])):
    """
        Returns the gauss curvature. Not used cause of bad performance
    """
    coeff_grad, powers_grad = jac_polyN_param(coeff, powers)
    coeff_hessian, powers_hessian = hessian_polyN_param(coeff_grad, powers_grad)

    def f(S):
        return(polyN(S, coeff))
                
    def grad_f(S):
        return(grad_polyN(S, coeff_grad, powers_grad))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))
    
    return(gauss_curv(f, grad_f, hessian_f, nb_pt_check))

def check_convexity_lpm(coeff, itermax, nb_pt):
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
        return(polyN(S, coeff))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))

    print("gen data")
    data = data_yf_sphere(f, itermax, nb_pt)
    hessian_data = hessian_f(data)
    L = leading_principal_minors(hessian_data)
    m = min(L.flatten())
    return(m)

def plot_implicit(yf, bbox=(-5,5)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 100) # resolution of the contour
    B = np.linspace(xmin, xmax, 20) # number of slices
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

def dc_kg(M, nb_dir):
    """
        Not used
    """
    nb_coeff = len(M)
    dirs = generate_dir(nb_dir, nb_coeff)
    nb_pt_check_1 = 10
    nb_pt_check_2 = 1000

    pt_set = np.empty((0, nb_coeff))
    kg_set = np.empty(0)

    for i in range(nb_dir):
        for j in [1, -1]:

            u = dirs[i] * j
            lamb = 1
            kg = -1
            itermax = 1
            it = 0

            """Find coefficients where polyN not convex by trial and error"""
            while kg < kg_max :
                coeff = M + lamb * u
                kg, kgm, kgmoy, n_pos, n_neg = check_convexity_kg(coeff, nb_pt_check_1)
                print(i, j, it, lamb)
                print("kg max", kg)
                print("kg min", kgm)
                print("kg moy", kgmoy)
                print("pos", n_pos)
                print("neg", n_neg)
                lamb = lamb * 1.5
                it = it + 1

            S = M
            E = coeff
            pt_set = np.concatenate((pt_set, np.atleast_2d(E)))
            kg_set = np.append(kg_set, kg)
            it = 0

            """Find coefficients such as polyN in on the border of the convexity domain"""
            while (kg > kg_max or kg < kg_min) and it < itermax:
                
                I = (S+E)/2
                kg, kgM, p, d, u = check_convexity_kg(I, nb_pt_check_2)
                print(kg, it)
                pt_set = np.concatenate((pt_set, np.atleast_2d(I)), axis=0)
                kg_set = np.append(kg_set, kg)
                if kg < 0 :
                    S = I
                else:
                    E = I
                it = it + 1

        return(pt_set, kg_set)
    

def dc_lpm(direction, M=M):
    """
        Returns a small dataset with coefficients and 0 if not convex, 1 if convex. 
        Comment : Designed for multiprocessing otherwise give only the number of directions, not the direction itself
        Input :
            - direction : ndarray of shape (nmon), direction on the unit sphere of dimension nmon
        Output :
            - pt_set : float ndarray of shape (n, 22), coefficients
            - convex_set : integer ndarray of shape (n, ), 0 if not convex, 1 if convex
    """
    nb_coeff = len(M)
    itermax = 5
    nb_pt_check_1 = 100   #For trial and error, just one pt is enough to prove non convexity
    nb_pt_check_2 = 1000    #To find a convex function, more points to test in the bisection method

    pt_set = np.empty((0, nb_coeff))
    #0 if not-convex, 1 if convex
    convex_set = np.empty(0, dtype=int)

    for j in [1, -1]:

        u = direction * j
        lamb = 1
        lpm = 1
        it = 0

        #Find coefficients where polyN not convex by trial and error
        while lpm > - tol_minor :
            coeff = M + lamb * u
            lpm = check_convexity_lpm(coeff, itermax, nb_pt_check_1)
            #print(i, j, it, lamb)
            print("lpm", lpm)
            lamb = lamb * 3
            it = it + 1

        #Find coefficients such as polyN in on the border of the convexity domain using bisection
        #with M the center of the convexity domain (Mises)
        S = M
        E = coeff
        pt_set = np.concatenate((pt_set, np.atleast_2d(E)))
        convex_set = np.append(convex_set, 0)
        it = 0
        
        #Find the point close to the border
        while (lpm < minor_min or lpm > minor_max):
                
            I = (S+E)/2
            print(I)

            lpm = check_convexity_lpm(I, itermax, nb_pt_check_2)
            pt_set = np.concatenate((pt_set, np.atleast_2d(I)), axis=0)
                
            if lpm > 0 :
                convex_set = np.append(convex_set, 1)
                S = I
            else:
                convex_set = np.append(convex_set, 0)
                E = I
            it = it + 1

        #Due to numerical round error, theres always negative lpm. So the last 
        #coefficients to be found corresponds to a polyN convex
        convex_set[-1] = 1

    return(pt_set, convex_set)


if __name__ == "__main__":

    file_dir = os.path.dirname(os.path.abspath(__file__))
    sep = os.sep

    t = time.time()
    dirs = generate_dir(nb_dir, nb_coeff)

    pool = multiprocessing.Pool()

    results = pool.map(dc_lpm, dirs)

    pool.close()
    pool.join()

    X = np.empty((0, nb_coeff))
    Y = np.empty(0)

    for i in range(nb_dir):
        X = np.concatenate((X, results[i][0]), axis = 0)
        Y = np.concatenate((Y, results[i][1]), axis = 0)
    
    np.save(file_dir + sep + "X_convex_{}.npy".format(degree), X)
    np.save(file_dir + sep + "Y_convex_{}.npy".format(degree), Y)

    print(time.time() - t, ":", len(X),"points")
    


"""print(max(K))
print(min(K))
plt.boxplot(K)
plt.show()

def g(S):
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    return(np.sum(S**2, axis = 1))

def grad_g(S):
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    return(2 * S)

def hessian_g(S):
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    A = 2 * np.eye(6,6)
    res = np.tile(A, (len(S), 1, 1))
    return(res)

K = gauss_curv(g, grad_g, hessian_g, 10000)
g_wrap = np.vectorize(lambda x, y, z : g(np.array([x, y, z, 0, 0, 0])))
print(max(K))
print(min(K))

X, Y = dc(M, nb_dir)
np.save("X_convex.npy", X)
np.save("Y_convex.npy", Y)"""