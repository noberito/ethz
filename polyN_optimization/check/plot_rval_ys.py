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
from optimize_polyN_mini import get_param_polyN_mini, f_min_squared, jac_polyN_2d_param, check_convexity
from optimize_polyN_mini import grad_polyN_2d, polyN_2d, get_dir_pst, get_yield_stress_SH, get_yield_stress_NT6, add_sh, add_nt6

def plot_check(df, coeff_polyN, powers, material, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel):
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
    filename = f"ysrval_{material}_poly{degree}_{weight_ut}_{weight_e2}_{weight_vir}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)

def plot_planestress(material, coeff, powers, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel):
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
    filename = f"planexy_{material}_poly{degree}_{weight_ut}_{weight_e2}_{weight_vir}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)


def ys_ut_mini(thetas, coeff, powers):
    """
        For given loading angles, returns the stress value for a uniaxial tensile test in the direction theta
        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff : ndarray of shape (nmon + 2,), coefficients of the yield function
            - powers : ndarray of shape (nmon, 3), powers of the polyN function
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

def rval_ut_mini(thetas, coeff_mini, powers):
    """
        For given loading angles, returns the rvalue for a uniaxial tensile test in the direction theta
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

def plot_rval_ut_mini(df, coeff, powers):
    """
        Plot the rvalues for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["Rval"] and ["LoadAngle"] for UTs
            - coeff : ndarray of shape (nmon + 2,), coeff of the polyN mini function
            - powers : ndarray of shape (nmon, 3), powers of the polyN mini function
    """
    
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi/2, n_thetas)
    r_vals_theo = rval_ut_mini(thetas_theo, coeff[:-2], powers)
    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    r_vals_exp = df["Rval"].iloc[index_rval].values
    thetas_exp = df["LoadAngle"].iloc[index_rval].values
    plt.plot(thetas_theo, r_vals_theo, c="red", label="R_val model")
    plt.scatter(thetas_exp, r_vals_exp, c="red", marker="x", label="R_val exp")
    plt.show()

def plot_yield_stresses_ut_mini(df, coeff, powers):
    """
        Plot the ys for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["YieldStress"], ["Rval"] and ["LoadAngle"] for UTs
            - coeff : ndarray of shape (nmon + 2,), coeff of the polyN mini function
            - powers : ndarray of shape (nmon, 3), powers of the polyN mini function
    """

    def f(S):
        return(f_min_squared(S, coeff, powers))
    
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
    
    ntheta=100
    thetas=np.linspace(0, np.pi / 2, ntheta)
    ys_theta = np.vectorize(ys_theta)
    ys_model = ys_theta(thetas)

    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    fig, ax1 = plt.subplots()

    ys_exp = df_exp["YieldStress"]
    theta_exp = df_exp["LoadAngle"].values

    plt.scatter(theta_exp, ys_exp/sigma0, color="blue", label="YS exp", marker="x")
    plt.xlabel(r"$\theta$[rad]", size=12)
    plt.ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)
    plt.plot(thetas, ys_model, color="blue", label="YS model")
    plt.show()

def plot_check_mini(df, coeff, powers, material, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel):
    """
        Plot the rvalues and ys for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["YieldStress"], ["Rval"] and ["LoadAngle"] for UTs
            - coeff : ndarray of shape (nmon + 2,), coeff of the polyN mini function
            - powers : ndarray of shape (nmon, 3), powers of the polyN mini function
            - rest of the usual parameters for plot's title
    """
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi / 2, n_thetas)

    ys_model = ys_ut_mini(thetas_theo, coeff, powers)
    r_vals_model = rval_ut_mini(thetas_theo, coeff[:-2], powers)

    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    ys_exp = df_exp["YieldStress"]

    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    thetas_exp = df["LoadAngle"].iloc[index_rval].values
    r_vals_exp = df["Rval"].iloc[index_rval].values

    plt.plot(thetas_theo, r_vals_model, c="red", label="R_val model")
    plt.plot(thetas_theo, ys_model, color="blue", label="YS model")
    plt.scatter(thetas_exp, r_vals_exp, c="red", marker="x", label="R_val exp")
    plt.scatter(thetas_exp, ys_exp/sigma0, color="blue", label="YS exp", marker="x")
    plt.title("Check poly")
    plt.xlabel(r"$\theta$[rad]", size=12)
    plt.ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.grid(1)
    plt.title(f"Poly{degree}_mini model on {material} : Yield Stresses and R-values from UT tests", size=12)

    foldername_out = polyN_dir + sep + "plots" + sep + material
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    filename = f"ysrval_{material}_poly{degree}mini_{weight_ut}_{weight_e2}_{weight_vir}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)
    plt.show()

def plot_planestress_mini(material, coeff, powers, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel):
    """
        Plot the yield surface in the plane sx sy for different sxy
    """
    zs = [0.0, 0.3, 0.4, 0.5, 0.6, 0.63]
    mises_plot = False
    fig, ax = plt.subplots()


    sx = np.linspace(-1.5, 1.5, 100)
    sy = np.linspace(-1.5, 1.5, 100)
    sx, sy = np.meshgrid(sx,sy)

    for z in zs:
        def mises_plane(x,y):
            return(mises(np.array([x,y,0,z,0,0])))

        def f(x,y):
            return(f_min_squared(np.array([x,y,0,z,0,0]), coeff, powers))

        f = np.vectorize(f)
        mises_plane = np.vectorize(mises_plane)

        ys_polyN = f(sx, sy)
        ys_mises = mises_plane(sx, sy)

        # Create the contour plot
        cs1 = ax.contour(sx, sy, ys_polyN, levels=[1], linewidths=1)
        ax.clabel(cs1, inline=1, fontsize=10, fmt={1: z})
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
    plt.legend(nm1, ["polyN_mini"])
    if mises_plot:
        plt.legend(nm2, ["Mises"])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    foldername_out = polyN_dir + sep + "plots" + sep + material
    if not os.path.exists(foldername_out):
        os.makedirs(foldername_out)
    filename = f"planexy_{material}_poly{degree}mini_{weight_ut}_{weight_e2}_{weight_vir}_{nb_virtual_pt}_{protomodel}.png"
    filepath = foldername_out + sep + filename
    plt.savefig(filepath)

def check_pst_points(df, material, coeff, powers):
    sigma0 = df["YieldStress"].iloc[0]
    ys_ratio_nt6 = get_yield_stress_NT6(sigma0, material)
    dir_pst = get_dir_pst(coeff, powers)
    new_df = add_nt6(df, ys_ratio_nt6, dir_pst)

    zs = [0]
    mises_plot = False
    fig, ax = plt.subplots()

    try:
        ebt = df[df["q"] == 1][["s11","s22"]].values[0]
        plt.scatter(ebt[0], ebt[1], marker="x", label = "UT_EBT")
    except IndexError:
        print("EBT non disponible")
        
    ut0 = df[df["LoadAngle"] == 0.0][["s11","s22"]].values[0]
    plt.scatter(ut0[0], ut0[1], marker="x", label = "UT_00")
    ut90 = df[df["LoadAngle"] > np.pi/2.1][["s11","s22"]].values[0]
    plt.scatter(ut90[0], ut90[1], marker="x", label = "UT_90")

    sx = np.linspace(-1.5, 1.5, 100)
    sy = np.linspace(-1.5, 1.5, 100)
    sx, sy = np.meshgrid(sx,sy)

    for z in zs:
        def mises_plane(x,y):
            return(mises(np.array([x,y,0,z,0,0])))

        def f(x,y):
            return(f_min_squared(np.array([x,y,0,z,0,0]), coeff, powers))

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
    nt6 = new_df[new_df["Type"] == "e2"][["s11", "s22"]].values

    for i in range(len(nt6)):
        
        plt.scatter(nt6[i,0], nt6[i,1], marker="x", label = i)

    plt.legend(nm1, ["polyN_mini"])
    plt.title(f"{material} : PST points")
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

def check_sh_points(df, material, coeff, powers):
    sigma0 = df["YieldStress"].iloc[0]
    ys_ratio_sh = get_yield_stress_SH(sigma0, material)
    new_df = add_sh(df, ys_ratio_sh)

    fig, ax = plt.subplots()

    sx = np.linspace(-1.5, 1.5, 100)
    sy = np.linspace(-1.5, 1.5, 100)
    sx, sy = np.meshgrid(sx,sy)

    def f(x,y):
        return(f_min_squared(np.array([x * 1 / np.sqrt(2),x * 1 / np.sqrt(2),0,y,0,0]), coeff, powers))

    f = np.vectorize(f)

    ys_polyN = f(sx, sy)

    # Create the contour plot
    cs1 = ax.contour(sx, sy, ys_polyN, levels=[1], colors='blue', linewidths=1)

    # Set labels
    ax.set_xlabel(r"$ \frac{1}{\sqrt{2}}(\sigma_{xx} + \sigma_{yy})/\sigma_0$[-]", size=12)
    ax.set_ylabel(r"$\sigma_{xy}/\sigma_0$[-]", size=12)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(1)
    
    nm1, labels = cs1.legend_elements()

    sh = new_df[new_df["Type"] == "e2"][["s11", "s22", "s12"]]
    sh["sbis"] = 1/ np.sqrt(2) * (sh["s11"] + sh["s22"])
    sh = sh[["sbis", "s12"]].values
    
    for i in range(len(sh)):
        
        plt.scatter(sh[i,0], sh[i,1], marker="x", label = i)

    plt.legend(nm1, ["polyN_mini"])
    plt.title(f"{material} : SH_000 and SH_090 points")
    plt.legend()
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.show()

def check_all_pt(df, material, coeff, powers, bbox=(-1.5,1.5)):

    xmin, xmax, ymin, ymax, zmin, zmax = bbox*3
    fig, ax = plt.subplots(nrows = 1, ncols = 1, subplot_kw={"projection":"3d"})
    A = np.linspace(xmin, xmax, 20) # resolution of the contour
    B = np.linspace(xmin, xmax, 20) # number of slices
    A1,A2 = np.meshgrid(A,A) # grid on which the contour is plotted

    sigma0 = df["YieldStress"].iloc[0]
    print("Getting NT6 points and SH points")
    ys_ratio_nt6 = get_yield_stress_NT6(sigma0, material)
    dir_pst = get_dir_pst(coeff, powers)
    new_df = add_nt6(df, ys_ratio_nt6, dir_pst)
    ys_ratio_sh = get_yield_stress_SH(sigma0, material)
    new_df = add_sh(new_df, ys_ratio_sh)

    pts = new_df[new_df["Type"]!= "v"][["s11", "s22", "s12"]].values

    X = pts[:,0]
    Y = pts[:,1]
    Z = pts[:,2]

    ax.scatter(X,Y,Z, marker="x", color="red", s=100)

    def f(S):
        return f_min_squared(S, coeff, powers)

    f_plane = lambda x, y, z : f(np.array([x, y, 0, z, 0, 0]))
        
    f_plane = np.vectorize(f_plane)

    print("Début du plot sur XY")
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = f_plane(X,Y,z) - 1
        ax.contour(X, Y, Z+z, [z], zdir='z', alpha=0.5)
        # [z] defines the only level to plot for this contour for this value of z
    print("Fin du plot sur XY")
    print("Début du plot sur XZ")
    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = f_plane(X,y,Z) - 1
        ax.contour(X, Y+y, Z, [y], zdir='y', alpha=0.5)

    print("Fin du plot sur XZ")
    print("Début du plot sur YZ")
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = f_plane(x,Y,Z) - 1
        ax.contour(X+x, Y, Z, [x], zdir='x', alpha=0.5)
    print("Fin du plot sur YZ")
    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)
    
    
    plt.show()



def main():
    p = read_param()
    func = p["func"]
    material = p["material"]
    degree = int(p["degree"])
    protomodel = p["protomodel"]
    weight_ut = float(p["weight_ut"])
    weight_e2 = float(p["weight_e2"])
    weight_vir = float(p["weight_vir"])
    nb_virtual_pt = int(p["nb_virtual_pt"])
    var_optim = p["var_optim"]

    df = readData(material, protomodel)

    if func == "polyN":
        coeff_polyN = np.load(polyN_dir + sep + "poly_coeff.npy")
        powers = get_param_polyN(degree)
        plot_check(df, coeff_polyN, powers, material, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel)
        plot_planestress(material, coeff_polyN, powers, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel)
        

    if func == "polyN_mini":

        if 1:
            coeff_file_initial = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy"
            coeff_polyN_mini = np.load(coeff_file_initial)

            if 1:
                coeff_file_final = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff_final_scipy_" + "[8]" + ".npy"
                coeff_polyN_mini = np.load(coeff_file_final)

            powers = get_param_polyN_mini(degree)
            #check_all_pt(df, material, coeff_polyN_mini, powers)
            check_sh_points(df, material, coeff_polyN_mini, powers)
            check_pst_points(df, material, coeff_polyN_mini, powers)
            plot_check_mini(df, coeff_polyN_mini, powers, material, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel)
            #plot_planestress_mini(material, coeff_polyN_mini, powers, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel)

        else :
            coeff_file_initial = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy"
            coeff_polyN_mini = np.load(coeff_file_initial)

            if 1:
                coeff_file_final = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff_final_scipy_" + str(var_optim) + ".npy"
                coeff_polyN_mini = np.load(coeff_file_final)

            powers = get_param_polyN_mini(degree)


main()