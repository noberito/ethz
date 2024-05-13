from tabnanny import check
import numpy as np
import sklearn
import sklearn.preprocessing
import matplotlib.pyplot as plt

degree = 2
nb_dir = 1
kg_min = 1e-5
kg_max = 5e-4

if degree==4:
    M = np.array([9, 18, 27, 18, 9, 18, 18, 18, 9, 18, 18, 18, 18, 9, 0, 0, 18, 18, 18, 18, 18, 9])

"""-----------------------------------PARAMETERS FOR POLYN------------------------------------------"""
def get_param_polyN(degree):
    data = np.zeros((2,5))
    polyN = sklearn.preprocessing.PolynomialFeatures((degree, degree), include_bias=False)
    X = polyN.fit_transform(data)
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
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = S.copy()
    trace_dev = S[:,0] + S[:,1] + S[:,2]
    D[:,0] = S[:,0] - (1/3) * trace_dev
    D[:,1] = S[:,1] - (1/3) * trace_dev
    D[:,2] = S[:,2] - (1/3) * trace_dev
    return(D)

def polyN(S, coeff):
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

def jac_dev(S):
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
            coefficients of each monomial derived with respect to dev[i]¨
            - powers_grad (float ndarray of shape (5, nmon, 5)) : powers_grad[i][j] contains
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
            - S (float ndarray of shape : len(data) * 6) : Stress
            - coeff_grad : coeff_grad (float ndarray of shape (5, nmon)) : coeff_grad[i] contains the 
            coefficients of each monomial derived with respect to dev[i]¨
            - powers_grad (float ndarray of shape (5, nmon, 5)) : powers_grad[i][j] contains
            the monomial j derived with respect to dev[i]
    """
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = dev(S)
    X = D[:,[0, 1, 3, 4, 5]]

    grad_polyN = np.zeros(6)
    grad_f = np.zeros((len(X),5))

    for i in range(5):
        for j in range(nmon):
            p = powers_grad[i][j]
            grad_f[:,i] = grad_f[:,i] + coeff_grad[i][j] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev(S)

    grad_polyN = np.dot(grad_f, jac_d)
    return(grad_polyN)

def hessian_polyN_param(coeff_grad, powers_grad):
    """
    Compute the different parameters and coefficients to compute the Hessian of polyN
        Input :
            - coeff_grad (float ndarray of shape (5, nmon)) : coeff_grad[i] contains the 
            coefficients of each monomial derived with respect to dev[i]
            - powers_grad (float ndarray of shape (5, nmon, 5)) : powers_grad[i][j] contains
            the powers of the monomial j derived with respect to dev[i]
        
        Output :
            - coeff_hessian (float ndarray of shape (5, 5, nmon)) : coeff_hessian[i][j][k] contains
            the coefficients of the monomial k of polyN derived with respect to dev[i] and then to dev[j]
            - powers_hessian (float ndarray of shape (5, 5, nmon, 5)) : powers_hessian[i][j][k] contains
            the powers of the monomial k of polyN derived with respect to dev[i] and then to dev[j]

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
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = dev(S)
    X = D[:,[0, 1, 3, 4, 5]]

    hessian_polyN = np.zeros((6,6))
    jac_grad_f = np.zeros((len(X), 5, 5))

    for i in range(5):
        for j in range(5):
            for k in range(nmon):
                p = powers_hessian[i][j][k]
                jac_grad_f[:,i,j] = jac_grad_f[:,i,j] + coeff_hessian[i][j][k] * np.prod(X ** p, axis=1)
    
    jac_d = jac_dev(S)
    hessian_polyN = np.transpose(np.dot(jac_d.T,np.dot(jac_grad_f[:], jac_d)), (1, 2, 0))
    return(hessian_polyN)
"""-----------------------------------------------------------------------------------------"""

def generate_dir(n, dim):
    np.random.seed(6)
    dirs = np.random.normal(0, 1, (n, dim))
    norms = np.linalg.norm(dirs, axis=1)
    dirs = (dirs.T/norms).T
    return(dirs)

def cofactor_matrix(lmatrix):
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

def gauss_curv(f, grad_f, hessian_f, nb_pt_check):
    dirs = generate_dir(nb_pt_check, 6)
    tol = 0.01
    itermax = 20
    data = np.zeros((nb_pt_check, 6))


    for i in range(nb_pt_check):
        E = np.zeros(6)
        m = -1
        u = dirs[i]
        it = 0
        lamb = 1

        while m < tol and it < itermax :
            E = E + lamb * u
            m = f(E) - 1
            lamb = lamb * 1.001
            it = it + 1

        S = np.zeros(6)
        res = f(E) - 1
        it = 0
        while abs(res) > tol and it < itermax:
            I = (S+E)/2
            res = f(I) - 1
            if res > 0 :
                E = I
            else:
                S = I
            it = it + 1
        data[i] = I
    
    K = - np.ones(nb_pt_check)
    grad_data = grad_f(data)
    norm_grad_data = np.linalg.norm(grad_data, axis = 1)
    hessian_data = hessian_f(data)

    cofactor_hessian_data = cofactor_matrix(hessian_data)

    for i in range(nb_pt_check):
        K[i] = K[i] // (norm_grad_data[i] ** 7) * np.dot(grad_data[i], np.dot(grad_data[i], cofactor_hessian_data[i]))

    return(K)

def check_convexity(coeff, nb_pt_check):
    coeff_grad, powers_grad = jac_polyN_param(coeff, powers)
    coeff_hessian, powers_hessian = hessian_polyN_param(coeff_grad, powers_grad)

    def f(S):
        return(polyN(S, coeff))
                
    def grad_f(S):
        return(grad_polyN(S, coeff_grad, powers_grad))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))

    m = np.min(gauss_curv(f, grad_f, hessian_f, nb_pt_check))

    return(m)

def gauss_curv_coeff(coeff, nb_pt_check):
    coeff_grad, powers_grad = jac_polyN_param(coeff, powers)
    coeff_hessian, powers_hessian = hessian_polyN_param(coeff_grad, powers_grad)

    def f(S):
        return(polyN(S, coeff))
                
    def grad_f(S):
        return(grad_polyN(S, coeff_grad, powers_grad))
                
    def hessian_f(S):
        return(hessian_polyN(S, coeff_hessian, powers_hessian))
    
    return(gauss_curv(f, grad_f, hessian_f, nb_pt_check))


def plot_implicit2(yf, bbox=(-1.5,1.5)):
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

def dc(M, nb_dir, powers):
    nb_coeff = len(M)
    dirs = generate_dir(nb_dir, nb_coeff)
    tol = 10e-4

    nb_pt_check_1 = 10
    nb_pt_check_2 = 10000

    pt_set = np.empty((0, nb_coeff))
    kg_set = np.empty(0)

    for i in range(nb_dir):
        for j in [1, -1]:

            u = dirs[i] * j
            lamb = 1
            kg = 1
            itermax = 10000
            it = 0

            while kg > - tol :
                coeff = M + lamb * u
                kg = check_convexity(coeff, nb_pt_check_1)

                print(i, j, it, lamb)
                print("kg", kg)

                lamb = lamb * 2
                it = it + 1

            S = M
            E = coeff

            pt_set = np.concatenate((pt_set, np.atleast_2d(E)))
            kg_set = np.append(kg_set, kg)

            it = 0

            while (kg > kg_max or kg < kg_min) and it < itermax:

                print(kg, it)

                I = (S+E)/2
                
                kg = check_convexity(I, nb_pt_check_2)

                pt_set = np.concatenate((pt_set, np.atleast_2d(I)), axis=0)
                kg_set = np.append(kg_set, kg)

                if kg > 0 :
                    S = I
                else:
                    E = I

                it = it + 1

            """def f(S):
                return(polyN(S, coeff))
            f_plane_wrap = np.vectorize(lambda x, y, z : f(np.array([x, y, 0, z, 0, 0])))
            plot_implicit(f_plane_wrap)"""

        
        return(pt_set, kg_set)


"""coeff = np.ones(22)
coeff[1] = -1
coeff[2] = -1
coeff[]
K = gauss_curv_coeff(coeff, 100000)
def f(S):
        return(polyN(S, coeff))
f_wrap = np.vectorize(lambda x, y, z : f(np.array([x, y, 0, z, 0, 0])))

print(max(K))
print(min(K))
plt.boxplot(K)
plt.show()
plot_implicit(f_wrap)"""

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
plot_implicit(g_wrap)