# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
# ---

import numpy as np
import pandas as pd
import sklearn 
import sklearn.preprocessing
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import multiprocessing

degree = 4
bbox=(-1.5,1.5)

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

specific = True

start = 200
end = 300

X_coeff = np.load("X_convex_{}.npy".format(degree))
Y_coeff = np.load("Y_convex_{}.npy".format(degree))

if specific :
    indexes = Y_coeff > 0.5
else:
    indexes = range(start, end)

print(indexes)


X_coeff = X_coeff[indexes]
Y_coeff = Y_coeff[indexes]

n_samples = len(X_coeff)

j_max = (n_samples - 1) // 3 + 1
k_max = 3

n_plots = 3

# fig, ax = plt.subplots(j_max, k_max,subplot_kw={'projection':'3d'})
# ax = np.atleast_2d(ax)

def plot_subplot(i):

    fig, ax = plt.subplots(nrows=1, ncols=n_plots, subplot_kw={'projection':'3d'})
    coeff = X_coeff[i]
    convex_bool = Y_coeff[i]

    if convex_bool: 
        plt.title("{} : CONVEX".format(indexes[i]))
    else:
        plt.title("{} : NOT CONVEX".format(indexes[i]))

    def f(S):
        return polyN(S, coeff)
    
    for j in range(n_plots):
    
        if j == 0:
            f_plane = lambda x, y, z : f(np.array([x, y, 0, z, 0, 0]))
        elif j == 1:
            f_plane = lambda x, y, z : f(np.array([x, y, 0, 0, z, 0]))
        else :
            f_plane = lambda x, y, z : f(np.array([x, y, 0, 0, 0, z]))

        f_plane = np.vectorize(f_plane)

        xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
        A = np.linspace(xmin, xmax, 50) # resolution of the contour
        B = np.linspace(xmin, xmax, 15) # number of slices
        A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

        print("Début du plot sur XY du plan {}".format(j))
        for z in B: 
            X,Y = A1,A2
            Z = f_plane(X,Y,z) - 1
            ax[j].contour(X, Y, Z+z, [z], zdir='z')
        print("Fin du plot sur XY du plan {}".format(j))

        print("Début du plot sur XZ du plan {}".format(j))
        for y in B: 
            X,Z = A1,A2
            Y = f_plane(X,y,Z) - 1
            ax[j].contour(X, Y+y, Z, [y], zdir='y')
        print("Fin du plot sur XZ du plan {}".format(j))

        print("Début du plot sur YZ du plan {}".format(j))
        for x in B: 
            Y,Z = A1,A2
            X = f_plane(x,Y,Z) - 1
            ax[j].contour(X+x, Y, Z, [x], zdir='x')
        print("Fin du plot sur YZ du plan {}".format(j))

        ax[j].set_zlim3d(zmin,zmax)
        ax[j].set_xlim3d(xmin,xmax)
        ax[j].set_ylim3d(ymin,ymax)

    plt.show()

if __name__ == "__main__":

    # Create a pool of processes
    pool = multiprocessing.Pool()

    results = pool.map(plot_subplot, range(n_samples))
    
    pool.close()
    pool.join()
