import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import os
import sklearn
import sklearn.preprocessing 

material = "DP600"
degree = 4
weigth_exp = 0.9
protomodel = "mises"

current_dir = "./"  # Assuming current directory
dir = "/"
if os.name == "nt":
    current_dir = ".\\"
    dir = "\\"

def readData(material, protomodel):
    filename_v = "data_virtual_" + material + "_" + protomodel + ".csv"
    filepath_v = current_dir + material + "_results" + dir + "DATA" + dir + filename_v

    db_v = pd.read_csv(filepath_v)

    filename_e = "data_cal_" + material + ".csv"
    filepath_e = current_dir + material + "_results" + dir + "DATA" + dir + filename_e

    db_e = pd.read_csv(filepath_e)

    sigma0 = db_e["YieldStress"].iloc[0]

    db_e["s11"] = db_e["YieldStress"] / sigma0 * np.square(np.cos(db_e["LoadAngle"])) + db_e["q"] * np.square(np.sin(db_e["LoadAngle"]))
    db_e["s22"] = db_e["YieldStress"] / sigma0 * np.square(np.sin(db_e["LoadAngle"])) + db_e["q"] * np.square(np.cos(db_e["LoadAngle"]))
    db_e["s33"] = db_e["YieldStress"] * 0
    db_e["s12"] = db_e["YieldStress"] / sigma0 * (1 - db_e["q"]) * np.sin(db_e["LoadAngle"]) * np.cos(db_e["LoadAngle"]) / sigma0
    db_e["s13"] = db_e["YieldStress"] * 0
    db_e["s23"] = db_e["YieldStress"] * 0

    db = db_e.loc[:, ["s11", "s22", "s33", "s12", "s13", "s23", "Type"]]
    db = pd.concat([db, db_v])

    db["d11"] = (2/3) * db["s11"] - (1/3) * db["s22"] - (1/3) * db["s33"]
    db["d22"] = - (1/3) * db["s11"] + (2/3) * db["s22"] - (1/3) * db["s33"]

    return(db)

db = readData(material, protomodel)
nb_virtual_pt = len(db[db["Type"] == "v"])
print(nb_virtual_pt)
data = db[["d11", "d22", "s12", "s13", "s23"]].values

def dev(S):
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = S.copy()
    trace_dev = S[:,0] + S[:,1] + S[:,2]
    D[:,0] = S[:,0] - (1/3) * trace_dev
    D[:,1] = S[:,1] - (1/3) * trace_dev
    D[:,2] = S[:,2] - (1/3) * trace_dev
    return(D)

def optiCoeff(data, degree, weigth_exp):
    polyN = sklearn.preprocessing.PolynomialFeatures((degree, degree), include_bias=False)
    X = polyN.fit_transform(data[:, :5],)
    powers = polyN.powers_

    selected_indices = []
    for i in range(len(powers)):
        k,l,m,n,p = powers[i]
        if ((m==n) and (n==p)) or ((m%2==0) and (n%2==0) and (p%2==0)):
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
    powers = powers.T[sorted_indices]
    X = X[:,sorted_indices]

    weigth = np.where(db["Type"] == "e", weigth_exp, 1-weigth_exp)
    ndata, nmon = X.shape
    M = np.zeros((nmon, nmon))

    for i in range(ndata):
        Xi = np.expand_dims(X[i], axis=0)
        M = M + weigth[i] * np.dot(Xi.T, Xi)
    M = M / 2

    V = np.sum(np.expand_dims(weigth, axis=1) * X, axis = 0)

    D = np.sum(weigth, axis = 0)
    print(M.shape, V.shape, D.shape)

    def J(a):
        return np.dot(np.dot(a, M), a) - np.dot(V, a) + D

    def Grad_J(a):
        return 2 * np.dot(M, a) - V

    a = np.zeros(nmon)
    a[0] = 1

    print("Optimisation en cours")
    opt = scipy.optimize.minimize(J, a, method='BFGS', jac=Grad_J)
    print("Optimisation terminée")
    coeff = opt.x

    return(coeff, powers, nmon)

coeff, powers, nmon = optiCoeff(data, degree, weigth_exp)

def polyN(S):
    if S.ndim==1:
        S = np.expand_dims(S, axis=0)
    D = dev(S)
    X = D[:,[0, 1, 3, 4, 5]]
    res = np.zeros(len(X))
    i = 0
    for i in range(nmon) :
        p = powers[i]
        res = res + coeff[i] * np.prod(X ** p, axis=1)
        
    return(res)

def polyN_plane(s11, s22, s12):
    S_plane = [s11, s22, s12]
    S = np.zeros(6)
    S[0] = S_plane[0]
    S[1] = S_plane[1]
    S[3] = S_plane[2]
    return(polyN(S))

polyN_plane = np.vectorize(polyN_plane)


def plot_implicit(yf, bbox=(-1.5,1.5)):
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

def export_coeff(coeff, protomodel, degree, material, nb_virtual_pt):
    filename = "deg{}_{}.txt".format(degree, protomodel)
    foldername = current_dir + material + "_results" + dir + "COEFF" + dir
    filepath = foldername + filename

    n = len(coeff)
    with open(filepath, "w") as file:
        file.write("#Coefficients poly{} for {} based on {} points from the {} protomodel\n".format(degree, material, nb_virtual_pt, protomodel))
        file.write("Number of coefficients : {}\n".format(n))
        for i in range(n):
            file.write("{} : {}\n".format(i + 1, coeff[i]))


data = db[["s11", "s22", "s33", "s12", "s13", "s23"]].values
S = np.zeros(6)
S[0] = 1
plot_implicit(polyN_plane)
export_coeff(coeff, protomodel, degree, material, nb_virtual_pt)