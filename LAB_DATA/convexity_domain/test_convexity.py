import numpy as np
import pandas as pd
import sklearn 
import sklearn.preprocessing
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

degree = 4
bbox=(-5,5)
"""----------------------------------------------DEF POLYN-----------------------------------------------------------------------------"""

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

X_coeff = np.load("X_convex.npy")
Y_coeff = np.load("Y_convex.npy")

print(X_coeff)
n_samples = len(X_coeff)

j_max = n_samples // 3 + 1
k_max = 3

fig, ax = plt.subplots(j_max, k_max,subplot_kw={'projection':'3d'})
ax = np.atleast_2d(ax)
print(ax.shape)
for i in range(n_samples):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''

    j, k = i // 3, i % 3
    coeff = X_coeff[i]
    print(coeff, j, k)

    def f(S):
        return polyN(S, coeff)
    
    f_plane = lambda x, y, z : f(np.array([x, y, 0, z, 0, 0]))
    f_plane = np.vectorize(f_plane)

    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    A = np.linspace(xmin, xmax, 100) # resolution of the contour
    B = np.linspace(xmin, xmax, 20) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    print("Début du plot sur XY")
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = f_plane(X,Y,z) - 1
        cset = ax[j, k].contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for this value of z
    print("Fin du plot sur XY")
    print("Début du plot sur XZ")
    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = f_plane(X,y,Z) - 1
        cset = ax[j, k].contour(X, Y+y, Z, [y], zdir='y')

    print("Fin du plot sur XZ")
    print("Début du plot sur YZ")
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = f_plane(x,Y,Z) - 1
        cset = ax[j, k].contour(X+x, Y, Z, [x], zdir='x')
    print("Fin du plot sur YZ")
    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax[j, k].set_zlim3d(zmin,zmax)
    ax[j, k].set_xlim3d(xmin,xmax)
    ax[j, k].set_ylim3d(ymin,ymax)

    val = Y_coeff[i]
    if val: 
        ax[j, k].set_title("CONVEX")
    else:
        ax[j, k].set_title("NOT CONVEX")

plt.show()