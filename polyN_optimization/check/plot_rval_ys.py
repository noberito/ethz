from email.utils import collapse_rfc2231_value
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys
import os

file_dir = os.path.dirname(os.path.abspath(__file__))
polyN_dir = os.path.dirname(file_dir)
sep = os.sep
sys.path.append(polyN_dir)

from optimize_polyN import get_param_polyN, polyN, jac_polyN_param, grad_polyN, readData, read_param, mises
from optimize_polyN_mini import get_param_polyN_mini, f_min_squared, jac_polyN_2d_param
from optimize_polyN_mini import grad_polyN_2d, polyN_2d, get_dir_pst, get_yield_stress_SH, get_yield_stress_NT6, add_sh, add_nt6
from optimize_Hill48 import hill48, grad_hill48, load_coeff_hill48
from optimize_YLD2000 import yld2000, grad_yld2000, load_coeff_yld2000

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
    coeff = coeff_mini[:-2]
    n_theta = len(thetas)
    S = np.zeros((n_theta, 6))
    v2 = np.zeros((n_theta, 3))
    r_val = np.zeros(n_theta)

    for i in range(n_theta):
        theta = thetas[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), 0, np.cos(theta) * np.sin(theta), 0, 0])
        v2[i] = np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])

    degree = np.sum(powers[0])
    coeff_grad, powers_grad = jac_polyN_2d_param(coeff, powers)
    grad_P = grad_polyN_2d(S, coeff_grad, powers_grad)[:,[0,1,3]]
    grad_f_plane = np.zeros((n_theta,3))

    for i in range(n_theta):
        for j in range(3):
            grad_f_plane[i,j] = (2 / degree) * grad_P[i,j] * np.float_power(polyN_2d(S[i], coeff_mini, powers),(2/degree - 1))[0]

    for i in range(n_theta):
        r_val[i] = - np.dot(grad_f_plane[i], v2[i]) / (grad_f_plane[i,0] + grad_f_plane[i,1])
    return(r_val)

def plot_rval_ut_mini(df, coeff_mini, powers):
    """
        Plot the rvalues for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["Rval"] and ["LoadAngle"] for UTs
            - coeff_mini : ndarray of shape (nmon + 2,), coeff of the polyN mini function
            - powers : ndarray of shape (nmon, 3), powers of the polyN mini function
    """
    
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi/2, n_thetas)
    r_vals_theo = rval_ut_mini(thetas_theo, coeff_mini, powers)
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

def plot_check_mini(df, coeff_mini, powers, material, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel):
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

    ys_model = ys_ut_mini(thetas_theo, coeff_mini, powers)
    r_vals_model = rval_ut_mini(thetas_theo, coeff_mini, powers)

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

def check_pst_points_mini(df, material, coeff, powers):
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
        print("EBT non available")
        
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

def check_sh_points_mini(df, material, coeff, powers):
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

def check_3D_mini(df, material, coeff, powers, bbox=(-1.5,1.5)):

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

    print("Start plotting XY")
    for z in B: # plot contours in the XY plane
        X,Y = A1,A2
        Z = f_plane(X,Y,z) - 1
        ax.contour(X, Y, Z+z, [z], zdir='z', alpha=0.5)
        # [z] defines the only level to plot for this contour for this value of z
    print("End plotting XY")
    print("Start plotting XZ")
    for y in B: # plot contours in the XZ plane
        X,Z = A1,A2
        Y = f_plane(X,y,Z) - 1
        ax.contour(X, Y+y, Z, [y], zdir='y', alpha=0.5)

    print("End plotting XZ")
    print("Start plotting YZ")
    for x in B: # plot contours in the YZ plane
        Y,Z = A1,A2
        X = f_plane(x,Y,Z) - 1
        ax.contour(X+x, Y, Z, [x], zdir='x', alpha=0.5)
    print("End plotting YZ")
    # must set plot limits because the contour will likely extend
    # way beyond the displayed level.  Otherwise matplotlib extends the plot limits
    # to encompass all values in the contour.
    ax.set_zlim3d(zmin,zmax)
    ax.set_xlim3d(xmin,xmax)
    ax.set_ylim3d(ymin,ymax)
    
    
    plt.show()

def ys_ut_hill48(thetas, coeff):
    """
        For given loading angles, returns the stress value for a uniaxial tensile test in the direction theta
        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff : ndarray of shape (nmon + 2,), coefficients of YLD2000
        Output :
            - cs : ndarray of shape (n_theta,)
    """
    n_theta = len(thetas)
    us = np.zeros((n_theta, 3))
    cs = []
    for i in range(n_theta):
        theta = thetas[i]
        us[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.cos(theta) * np.sin(theta)])
    
    def f(S):
        return(hill48(S, coeff))

    for u in us:
        ys = 0.1
        lamb = 1.1
        eps = 1e-7
        
        while f(u * ys) < 1:
            ys = ys * lamb
        s = 0.1
        e = ys
        c = (s + e) / 2
        res = f(u * c) - 1

        while (np.abs(res) > eps):
            if res > 0:
                e = c
            else:
                s = c
            c = (s + e) / 2
            res = f(u * c) - 1
        cs.append(c)
    return(np.array(cs))

def rval_ut_hill48(thetas, coeff_hill):
    """
        For given loading angles, returns the rvalue for a uniaxial tensile test in the direction theta
        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff_hill : ndarray of shape (nmon + 2,), coefficients of the Hill48 function
        Output :
            - rval : ndarray of shape (n_theta,)
    """
    n_theta = len(thetas)
    S = np.zeros((n_theta, 3))
    v2 = np.zeros((n_theta, 3))
    r_val = np.zeros(n_theta)
    
    for i in range(n_theta):
        theta = thetas[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.cos(theta) * np.sin(theta)])
        v2[i] = np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])

    grad_f_plane = grad_hill48(S, coeff_hill)

    for i in range(n_theta):
        r_val[i] = - np.dot(grad_f_plane[i], v2[i]) / (grad_f_plane[i,0] + grad_f_plane[i,1])

    return(r_val)

def plot_s_hill48(df, coeff, material):
    """
        Plot the rvalues and ys for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["YieldStress"], ["Rval"] and ["LoadAngle"] for UTs
            - coeff : ndarray of shape (nmon + 2,), coeff of the yld2000
            - material : string
            - m : degree
    """
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi / 2, n_thetas)

    ys_model = ys_ut_hill48(thetas_theo, coeff)

    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    ys_exp = df_exp["YieldStress"]

    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    thetas_exp = df["LoadAngle"].iloc[index_rval].values
    #r_vals_exp = df["Rval"].iloc[index_rval].values

    #plt.plot(thetas_theo, r_vals_model, c="red", label="R_val model")
    plt.plot(thetas_theo, ys_model, color="blue", label="YS model")
    #plt.scatter(thetas_exp, r_vals_exp, c="red", marker="x", label="R_val exp")
    plt.scatter(thetas_exp, ys_exp/sigma0, color="blue", label="YS exp", marker="x")
    plt.title("Check poly")
    plt.xlabel(r"$\theta$[rad]", size=12)
    plt.ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.grid(1)
    plt.title(f"YLD2000 model on {material} : Yield Stresses and R-values from UT tests", size=12)
    plt.show()

def ys_ut_yld2000(thetas, coeff, m):
    """
        For given loading angles, returns the stress value for a uniaxial tensile test in the direction theta
        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff : ndarray of shape (nmon + 2,), coefficients of YLD2000
            - m : degree
        Output :
            - cs : ndarray of shape (n_theta,)
    """
    n_theta = len(thetas)
    us = np.zeros((n_theta, 3))
    cs = []
    for i in range(n_theta):
        theta = thetas[i]
        us[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.cos(theta) * np.sin(theta)])
    
    def f(S):
        return(yld2000(S, coeff, m))

    for u in us:
        ys = 0.1
        lamb = 1.1
        eps = 1e-7
        
        while f(u * ys) < 1:
            ys = ys * lamb
        s = 0.1
        e = ys
        c = (s + e) / 2
        res = f(u * c) - 1

        while (np.abs(res) > eps):
            if res > 0:
                e = c
            else:
                s = c
            c = (s + e) / 2
            res = f(u * c) - 1
        cs.append(c)
    return(np.array(cs))

def rval_ut_yld2000(thetas, coeff_yld2000, m):
    """
        For given loading angles, returns the rvalue for a uniaxial tensile test in the direction theta
        Input :
            - thetas : ndarray of shape (n_theta,), angles to check
            - coeff_yld2000 : ndarray of shape (nmon + 2,), coefficients of the yld2000 function
        Output :
            - rval : ndarray of shape (n_theta,)
    """
    n_theta = len(thetas)
    S = np.zeros((n_theta, 3))
    v2 = np.zeros((n_theta, 3))
    r_val = np.zeros(n_theta)
    
    for i in range(n_theta):
        theta = thetas[i]
        S[i] = np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.cos(theta) * np.sin(theta)])
        v2[i] = np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.cos(theta) * np.sin(theta)])

    grad_f_plane = grad_yld2000(S, coeff_yld2000, m)

    for i in range(n_theta):
        r_val[i] = - np.dot(grad_f_plane[i], v2[i]) / (grad_f_plane[i,0] + grad_f_plane[i,1])

    return(r_val)

def plot_s_yld2000(df, coeff, material, m):
    """
        Plot the rvalues and ys for UTs for polyN mini
        Input :
            - df : Dataframe, must contains ["YieldStress"], ["Rval"] and ["LoadAngle"] for UTs
            - coeff : ndarray of shape (nmon + 2,), coeff of the yld2000
            - material : string
            - m : degree
    """
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi / 2, n_thetas)

    ys_model = ys_ut_yld2000(thetas_theo, coeff, m)

    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    ys_exp = df_exp["YieldStress"]

    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    thetas_exp = df["LoadAngle"].iloc[index_rval].values
    #r_vals_exp = df["Rval"].iloc[index_rval].values

    #plt.plot(thetas_theo, r_vals_model, c="red", label="R_val model")
    plt.plot(thetas_theo, ys_model, color="blue", label="YS model")
    #plt.scatter(thetas_exp, r_vals_exp, c="red", marker="x", label="R_val exp")
    plt.scatter(thetas_exp, ys_exp/sigma0, color="blue", label="YS exp", marker="x")
    plt.title("Check poly")
    plt.xlabel(r"$\theta$[rad]", size=12)
    plt.ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.grid(1)
    plt.title(f"YLD2000 model on {material} : Yield Stresses and R-values from UT tests", size=12)
    plt.show()

def get_data_prefile(material):
    """
        Load the ys ratios and r-values from the yld2000 prefile if available
    """
    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material 
    filename = "_YLD2000_pre.csv"
    filepath = results_exp_dir + sep + filename

    try:
        df = pd.read_csv(filepath)
        ys = df.loc[1].values[:-1]
        r_val = df.loc[2].values[:-1]
        return(ys, r_val)

    except:
        print(f"No prefile for {material}")
        return([], [])


def plot_ys_all(df, material, coeff_hill, coeff_yld, coeff_mini, m, powers, plot_hill, plot_yld, plot_mini):
    """
        Plot the ys ratios and r-values from the hill'48 and/or yld2000 and/or polyN_mini model.

        Input :

            - df : dataFrame, data extracted from the experimental file with readData2D
            - material : string
            - coeff_hill : ndarray of shape(6,), hill'48 model
            - coeff_yld2000 : ndarray of shape(8,), yld2000 model
            - coeff_mini : ndarray of shape(nmon + 2,), polyN minimalistic model
            - m : int, degree of yld2000 (usually 8)
            - powers : ndarray of shape(nmon,), powers of the polyN2D used in polyN minimalistic
            - plot_hill : 1 to plot hill if not 0
            - plot_yld2000 : 1 to plot yld2000 if not 0
            - plot_mini : 1 to plot polyN if not 0
    """
    degree = np.sum(powers[0])
    n_thetas = 100
    thetas_theo = np.linspace(0, np.pi / 2, n_thetas)

    ys_hill = ys_ut_hill48(thetas_theo, coeff_hill)
    ys_yld = ys_ut_yld2000(thetas_theo, coeff_yld, m)
    ys_polyN = ys_ut_mini(thetas_theo, coeff_mini, powers)

    r_vals_hill = rval_ut_hill48(thetas_theo, coeff_hill)
    r_vals_yld2000 = rval_ut_yld2000(thetas_theo, coeff_yld, m)
    r_vals_polyN = rval_ut_mini(thetas_theo, coeff_mini, powers)

    df_exp = df[df["Rval"] > 0.001]
    sigma0 = df_exp["YieldStress"].iloc[0]
    ys_exp = df_exp["YieldStress"]

    index_rval = np.where(df["Rval"]< 0.00001, False, True)
    thetas_exp = df["LoadAngle"].iloc[index_rval].values
    r_vals_exp = df["Rval"].iloc[index_rval].values

    ys_ratio_exp2, r_vals_exp2 = get_data_prefile(material)

    #plt.plot(thetas_theo, r_vals_model, c="red", label="R_val model")
    if plot_hill:
        plt.plot(thetas_theo, ys_hill, color="blue", linewidth=1, label="Hill'48")
        plt.plot(thetas_theo, r_vals_hill, color="blue", linestyle="dashed", linewidth=1)
    if plot_yld:
        plt.plot(thetas_theo, ys_yld, color="red", linewidth=1, label="Yld2000")
        plt.plot(thetas_theo, r_vals_yld2000, color="red", linestyle="dashed", linewidth=1)
    if plot_mini:
        plt.plot(thetas_theo, ys_polyN, color="black", linewidth=1, label=f"Poly{degree}")
        plt.plot(thetas_theo, r_vals_polyN, color="black", linestyle="dashed", linewidth=1)

    #plt.scatter(thetas_exp, r_vals_exp, c="red", marker="x", label="R_val exp")
    plt.scatter(thetas_exp, ys_exp/sigma0, color="black", label="Exp.", linewidths=1, marker="x")
    plt.scatter(thetas_exp, r_vals_exp, color="black", linewidths=1, marker="x")
    plt.scatter(thetas_exp, ys_ratio_exp2, color="red", label="Exp. 2", linewidths=1, marker="x")
    plt.scatter(thetas_exp, r_vals_exp2, color="red", linewidths=1, marker="x")
    plt.title("Check poly")
    plt.xlabel(r"$\theta$[rad]", size=12)
    plt.ylabel(r'$\sigma$ / $\sigma_0$[-]', size=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.grid(1)
    plt.title(f"{material} : Yield Stresses and R-values from UT tests", size=12)
    plt.show()  

def plot_planestress_all(material, coeff_hill, coeff_yld, coeff_polyN, m, powers, plot_mises, plot_hill, plot_yld, plot_mini):
    """
        Plot yield surfaces (Mises, Hill48, Yld2000, PolyN_mini) in the plane sx sy for sxy = 0

        Input :
            - material : string
            - coeff_hill : ndarray of shape(6,), hill'48 model
            - coeff_yld2000 : ndarray of shape(8,), yld2000 model
            - coeff_mini : ndarray of shape(nmon + 2,), polyN minimalistic model
            - m : int, degree of yld2000 (usually 8)
            - powers : ndarray of shape(nmon,), powers of the polyN2D used in polyN minimalistic
            - plot_mises : 1 to plot mises if not 0
            - plot_hill : 1 to plot hill if not 0
            - plot_yld2000 : 1 to plot yld2000 if not 0
            - plot_mini : 1 to plot polyN if not 0
    """

    degree = np.sum(powers[0])
    zs = [0.0]
    fig, ax = plt.subplots()

    sx = np.linspace(-1.5, 1.5, 100)
    sy = np.linspace(-1.5, 1.5, 100)
    sx, sy = np.meshgrid(sx, sy)

    for z in zs:
        def mises_plane(x, y):
            return mises(np.array([x, y, 0, z, 0, 0]))

        def hill48_plane(x, y):
            return hill48(np.array([x, y, z]), coeff_hill)

        def yld2000_plane(x, y):
            return yld2000(np.array([x, y, z]), coeff_yld, m)

        def polyN_mini_plane(x, y):
            return f_min_squared(np.array([x, y, 0, z, 0, 0]), coeff_polyN, powers)

        mises_plane = np.vectorize(mises_plane)
        hill48_plane = np.vectorize(hill48_plane)
        yld2000_plane = np.vectorize(yld2000_plane)
        polyN_mini_plane = np.vectorize(polyN_mini_plane)

        ys_mises = mises_plane(sx, sy)
        ys_hill = hill48_plane(sx, sy)
        ys_yld2000 = yld2000_plane(sx, sy)
        ys_polyN_mini = polyN_mini_plane(sx, sy)

        handles = []
        labels = []

        if plot_mises:
            cs1 = ax.contour(sx, sy, ys_mises, levels=[1], linewidths=1, colors="green")
            handles.append(Line2D([0], [0], color="green", lw=1))
            labels.append("Mises")

        if plot_hill:
            cs2 = ax.contour(sx, sy, ys_hill, levels=[1], linewidths=1, colors="red")
            handles.append(Line2D([0], [0], color="red", lw=1))
            labels.append("Hill48")

        if plot_yld:
            cs3 = ax.contour(sx, sy, ys_yld2000, levels=[1], linewidths=1, colors="blue")
            handles.append(Line2D([0], [0], color="blue", lw=1))
            labels.append("Yld2000")

        if plot_mini:
            cs4 = ax.contour(sx, sy, ys_polyN_mini, levels=[1], linewidths=1, colors="black")
            handles.append(Line2D([0], [0], color="black", lw=1))
            labels.append(f"Poly{degree}")

    ax.legend(handles, labels)

    # Set labels
    ax.set_xlabel(r"$\sigma_{xx}/\sigma_0$ [-]", size=12)
    ax.set_ylabel(r"$\sigma_{yy}/\sigma_0$ [-]", size=12)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True)
    ax.set_title(rf'{material} Yield surface in the $\sigma_{{xx}},\sigma_{{yy}}$ plane', size=12)

    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
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
        if 0:
            coeff_file_initial = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy"
            coeff_polyN_mini = np.load(coeff_file_initial)

            if 1:
                coeff_file_final = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff_final_scipy_" +  + ".npy"
                coeff_polyN_mini = np.load(coeff_file_final)

            powers = get_param_polyN_mini(degree)
            #check_all_pt(df, material, coeff_polyN_mini, powers)
            check_sh_points_mini(df, material, coeff_polyN_mini, powers)
            check_pst_points_mini(df, material, coeff_polyN_mini, powers)
            plot_check_mini(df, coeff_polyN_mini, powers, material, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel)
            #plot_planestress_mini(material, coeff_polyN_mini, powers, weight_ut, weight_e2, weight_vir, nb_virtual_pt, degree, protomodel)

        else :
            
            powers = get_param_polyN_mini(degree)
            coeff_file_initial = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy"
            coeff_polyN_mini = np.load(coeff_file_initial)

            if 0:
                coeff_file_final = polyN_dir + sep + material + "_poly" + str(degree) + "_mini_coeff_final_scipy_" + "[5]" + ".npy"
                coeff_polyN_mini = np.load(coeff_file_final)

            coeff_hill = load_coeff_hill48(material)
            coeff_yld = load_coeff_yld2000(material)

            plot_ys_all(df, material, coeff_hill, coeff_yld, coeff_polyN_mini, 8, powers, 1, 1, 1)
            plot_planestress_all(material, coeff_hill, coeff_yld, coeff_polyN_mini, 8, powers, 1, 1, 1, 1)


main()