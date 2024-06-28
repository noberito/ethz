import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

file_dir = os.path.dirname(os.path.abspath(__file__))
polyN_dir = os.path.dirname(file_dir)
sep = os.sep
sys.path.append(polyN_dir)

from optimize_polyN import get_param_polyN, polyN, jac_polyN_param, grad_polyN, readData, read_param, mises



def plot_check(df, coeff_polyN, powers, material, weight_exp, weight_rval, nb_virtual_pt, degree, protomodel):
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

    foldername_out = polyN_dir + sep + "plots" + sep + material
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    filename = f"ysrval_{material}_poly{degree}_{weight_exp}_{weight_rval}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)

def plot_planestress(material, coeff, powers, weight_exp, weight_rval, nb_virtual_pt, degree, protomodel):
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

    foldername_out = polyN_dir + sep + "plots" + sep + material
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    filename = f"planexy_{material}_poly{degree}_{weight_exp}_{weight_rval}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)

def plot_implicit_coeff(material, coeff, powers, weight_exp, weight_rval, nb_virtual_pt, degree, protomodel, bbox=(-1.5,1.5)):
    ''' create a plot of an implicit function
    fn  ...implicit function (plot where fn==0)
    bbox ..the x,y,and z limits of plotted interval'''

    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig, ax = plt.subplots(nrows = 1, ncols = 1, subplot_kw={"projection":"3d"})
    A = np.linspace(xmin, xmax, 100) # resolution of the contour
    B = np.linspace(xmin, xmax, 50) # number of slices
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

    foldername_out = polyN_dir + sep + "plots" + sep + material
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    filename = f"3dcuts_{material}_poly{degree}_{weight_exp}_{weight_rval}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)
    plt.show()

def main():
    p = read_param()
    material = p["material"]
    law = p["law"]
    degree = int(p["degree"])
    protomodel = p["protomodel"]
    input_type = p["input_type"]
    n_opti = int(p["n_opti"])
    density = p["density"]
    weight_exp = float(p["weight_exp"])
    weight_rval = float(p["weight_rval"])
    nb_virtual_pt = int(p["nb_virtual_pt"])
    powers = get_param_polyN(degree)
    df = readData(material, protomodel)

    coeff_polyN = np.load(polyN_dir + sep + "polyN_coeff.npy")

    plot_check(df, coeff_polyN, powers, material, weight_exp, weight_rval, nb_virtual_pt, degree, protomodel)
    plot_planestress(material, coeff_polyN, powers, weight_exp, weight_rval, nb_virtual_pt, degree, protomodel)
    plot_implicit_coeff(material, coeff_polyN, powers, weight_exp, weight_rval, nb_virtual_pt, degree, protomodel)

main()