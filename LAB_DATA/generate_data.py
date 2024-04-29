import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

material = "DP780"
nb_virtual_pt = 100000
protomodel = "mises"
current_dir = "./"  # Assuming current directory
dir = "/"
if os.name == "nt":
    current_dir = ".\\"
    dir = "\\"

def generate_dir(nb_virtual_pt):
    u = np.random.normal(0, 1, (nb_virtual_pt,6))
    return(u)


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
    return res - 1

def mises_plane(s11, s22, s12):
    """
    Vectorized Mises function
    """

    res = np.sqrt(0.5 * (s11 - s22)**2 + 3 * s12**2)
    return res - 1

# Find points where the function is close to zero
def plot_implicit(yf, bbox=(-3,3)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''
    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    A = np.linspace(xmin, xmax, 100) # resolution of the contour
    B = np.linspace(xmin, xmax, 10) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = yf(X,Y,z)
        cset = ax.contour(X, Y, Z+z, [z], zdir='z')
        # [z] defines the only level to plot for this contour for this value of z

    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = yf(X,y,Z)
        cset = ax.contour(X, Y+y, Z, [y], zdir='y')

    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = yf(x,Y,Z)
        cset = ax.contour(X+x, Y, Z, [x], zdir='x')

    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)

    plt.show()

#plot_implicit(mises_plane)
#TODO ca peut planter etre faux ici si 0 n'est pas compris dans l'intervalle des valeurs prises
def generate_data(nb_virtual_pt):
    np.random.seed(98)
    data = np.zeros((nb_virtual_pt, 6))
    us = generate_dir(nb_virtual_pt)
    alpha = np.linspace(0, 5, 100000)
    alpha = alpha[np.newaxis, :]
    for i in range(nb_virtual_pt):
        u = np.expand_dims(us[i], axis=1)
        sigmas = np.dot(u, alpha).T
        yss=mises(sigmas)
        #m = np.min(ys)
        #M = np.max(ys)
        #if not (m < 0 < M):
            #print("hey")
        k = np.argmin(np.abs(yss))
        data[i] = sigmas[k]
    return(data)


data = generate_data(nb_virtual_pt)


db = pd.DataFrame(data, columns=["s11", "s22", "s33", "s12", "s13", "s23"])
db["Type"] = ["v"] * len(data)

stress = mises(db[["s11", "s22", "s33", "s12", "s13", "s23"]])
print("Le max est : ", max(stress), ", Le min est : ", min(stress))

tol = 1e-3
print("Le nombre de points où l'écart à 0 est supérieur à ", tol," est : ", db[mises(db[["s11", "s22", "s33", "s12", "s13", "s23"]]) > 0.001].size)



filename = "data_virtual_" + material + "_" + protomodel + ".csv"
filepath = current_dir + material + "_results" + dir + "DATA" + dir + filename

db.to_csv(filepath, index=False)