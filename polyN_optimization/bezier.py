

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import binom
import random

from read import readdata_exp

polyN_dir = os.path.dirname(os.path.abspath(__file__))
sep = os.sep


n_pt_seg_dir = 10
"""
    For any question regarding the creation of the Bezier surface, please refer to the shyQP article written by
    Soare in 2023.
"""

#Functions and coefficients to calculate Bezier segment of order 5
c50 = binom(5,0)
c51 = binom(5,1)
c52 = binom(5,2)
c53 = binom(5,3)
c54 = binom(5,4)
c55 = binom(5,5)

def phi0(t):
    return c50 * np.power(1 - t,5)

def phi1(t):
    return c51 * np.power(t, 1) * np.power(1 - t,4)

def phi2(t):
    return c52 * np.power(t, 2) * np.power(1 - t,3)

def phi3(t):
    return c53 * np.power(t, 3) * np.power(1 - t,2)

def phi4(t):
    return c54 * np.power(t, 4) * np.power(1 - t,1)

def phi5(t):
    return c55 * np.power(t, 5)
    
def bezier_seg(bs, be, us, ue, ls, le, n):
    """
        Returns a Bezier segment of order 5 between bs and be. us (resp. ue) being the tangent
        at bs (resp. be) and ls (resp. le) a parameter on bs (resp. be)

        Input :
            - bs : ndarray of shape(m,), start segment
            - be : ndarray of shape(m,), end segment
            - us : ndarray of shape(m,), tangent at bs
            - ue : ndarray of shape(m,), tangent at be
            - ls : float, parameter at bs
            - le : float, parameter at be
            - n : int, number of points from the segment
        
        Output :
            - Y : ndarray of shape (n, m), segment
    """
    b0 = bs
    b1 = bs + ls * us
    b2 = bs + 2 * ls * us
    b3 = be - 2 * le * ue
    b4 = be - le * ue
    b5 = be

    def y(t):
        p = phi0(t) * b0 + phi1(t) * b1 + phi2(t) * b2 + phi3(t) * b3 + phi4(t) * b4 + phi5(t) * b5
        return(p)

    T = np.linspace(0, 1, n)
    Y = np.array([y(t) for t in T])

    return(Y)

def param_data_dir(material):
    """
        Returns all the parameters needed to create the segment containing
        the directionnal data.

        WARNING : To adapt if tension compression asymetry. I would suggest 
        to adapt readData_2D too.

        Input :
            - material : string
        
        Output :
            - YS : ndarray of shape((n_data, 2)), directionnal points in the space (theta, ys_ratio)
            - mus_ys : ndarray of shape((n_seg, 2)), starting tangents of every part of the yield stress segment
            - mue_ys : ndarray of shape((n_seg, 2)), ending tangents of every part of the yield stress segment
            - lmax_ys : int, lambda to have convex parts of the yield stress segment (hypothesis : lambda = cst for all the parts)
            - R : ndarray of shape((n_data, 2)), directionnal points in the space (theta, r)
            - mus_r : ndarray of shape((n_seg, 2)), starting tangents of every part of the r-value segment
            - mue_ys : ndarray of shape((n_seg, 2)), ending tangents of every part of the r-value segment
            - lmax_r : int, lambda to have convex parts of the r-value segment (hypothesis : lambda = cst for all the parts)
            - thetas : ndarray of shape (n_data), experimental thetas
            - ys_ratio : ndarray of shape (n_data), experimental yield stresses
            - r_val : ndarray of shape (n_data), experimental r-values
    """

    df = readdata_exp(material)
    data_dir = df[df["Type"] == "e"][df["q"] == 0]
    
    thetas = data_dir["LoadAngle"].values
    ys = data_dir["YieldStress"].values
    ys_ratio = ys / ys[0]

    r_val = data_dir["Rval"].values

    n_data = len(thetas)
    n_seg = n_data - 1

    YS = np.zeros((n_data, 2))
    YS[:,0] = thetas
    YS[:,1] = ys_ratio

    R = np.zeros((n_data, 2))
    R[:,0] = thetas
    R[:,1] = r_val

    mu = 1/2
    mus_ys = np.zeros((n_seg , 2))
    mus_ys[0] = np.array([1, 0])
    mus_r = np.zeros((n_seg, 2))
    mus_r[0] = np.array([1, 0])

    ys1, ys2, ys3 = YS[-1], YS[0], YS[1]
    r1, r2, r3 = R[-1], R[0], R[1]
    
    for i in range(1, n_seg):
        ys1 = ys2 
        ys2 = ys3 
        ys3 = YS[i + 1]

        r1 = r2 
        r2 = r3 
        r3 = R[i + 1]

        v_ys = (1 - mu) * (ys2 - ys1) + mu * (ys3 - ys2)
        t_ys = v_ys/np.linalg.norm(v_ys)

        mus_ys[i] = t_ys

        v_r = (1 - mu) * (r2 - r1) + mu * (r3 - r2)
        t_r = v_r/np.linalg.norm(v_r)

        mus_r[i] = t_r
    
    mue_ys = np.roll(mus_ys, -1, axis=0)
    mue_r = np.roll(mus_r, -1, axis=0)

    dtheta = thetas[1] - thetas[0]

    lmax_ys = np.min(dtheta/(2 * (mus_ys[:-1, 0] + mus_ys[1:, 0])))
    lmax_r = np.min(dtheta/(2 * (mus_r[:-1, 0] + mus_r[1:, 0])))

    return(YS, mus_ys, mue_ys, lmax_ys, R, mus_r, mue_r, lmax_r, thetas, ys_ratio, r_val)

def seg_data_dir(material, alpha):
    """
        Returns the directionnal segments (yield stress ratios, r-values)

        Input :
            - material : string
            - alpha : float between 0 and 1, 0 : Lowest curvature, 1: Highest curvature

        Output :
            - B_ys : ndarray of shape ((n * n_seg, 2)), yield stress directionnal segment in space (theta, ys_ratio)
            - B_r : ndarray of shape ((n * n_seg, 2)), r-values directionnal segment in space (theta, r-values)
    """

    YS, mus_ys, mue_ys, lmax_ys, R, mus_r, mue_r, lmax_r, _, _, _ = param_data_dir(material)
    l_ys = alpha * lmax_ys
    l_r = alpha * lmax_r
    n_seg = len(mus_ys)

    B_ys = np.zeros((n_pt_seg_dir * n_seg, 2))
    B_r = np.zeros((n_pt_seg_dir * n_seg, 2))

    for i in range(n_seg):

        B_ys[(i * n_pt_seg_dir):((i + 1) * n_pt_seg_dir)] = bezier_seg(YS[i], YS[i + 1], mus_ys[i], mue_ys[i], l_ys, l_ys, n_pt_seg_dir)
        B_r[(i * n_pt_seg_dir):((i + 1) * n_pt_seg_dir)] = bezier_seg(R[i], R[i + 1], mus_r[i], mue_r[i], l_r, l_r, n_pt_seg_dir)

    return(B_ys, B_r)



def bound_lambda(bs, be, us, ue):
    """
        Returns the maximum value for the starting (t/2) and ending lambda (tau/2) to have
        a convex curve in a plane PI(theta)

        Input :
            - bs : ndarray of shape(3,), start segment
            - be : ndarray of shape(3,), end segment
            - us : ndarray of shape(3,), tangent at bs
            - ue : ndarray of shape(3,), tangent at be

        Output :
            - halft : float, max value for ls
            - halftau : float, max value for le
    """
    v = np.dot(us, ue)
    u = (be - bs) /(1 - v * v) 

    t = np.dot(u, us - v * ue)
    tau = np.dot(u, ue - v * us)

    halft = t/2
    halftau = tau/2

    return(halft, halftau)


def param_bezier_section_plane(material, theta, ys, ys_s, r, r_s, r_tb=1, r_cb=1, tc_sym=1):
    df = readdata_exp(material)

    #Intersection of surface Pi(theta) and the yield surface
    theta_s = np.pi/2 - theta
    v = np.array([- np.sin(2 * theta), np.sin(2 * theta), 2 * np.cos(2 * theta)])

    p_t = ys * np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.sin(theta) * np.cos(theta)])
    p_ts = ys_s * np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.sin(theta) * np.cos(theta)])

    p_c = ys * np.array([- np.square(np.cos(theta)), - np.square(np.sin(theta)), - np.sin(theta) * np.cos(theta)])
    p_cs = ys_s * np.array([- np.square(np.sin(theta)), - np.square(np.cos(theta)), np.sin(theta) * np.cos(theta)])

    try:
        b_t = df[df["q"] == 1][["s11", "s22", "s12"]].values[0]
    except :
        #print("EBT non available, we use the approximation : sEBT = (s0 + 2s45 + s90) / 4")
        s0 = df[np.abs(df["LoadAngle"] - 0) < 10e-2].iloc[0]["YieldStress"]
        s45 = df[np.abs(df["LoadAngle"] - 45 * 2 * np.pi / 360) < 10e-2].iloc[0]["YieldStress"]
        s90 = df[np.abs(df["LoadAngle"] - 90 * 2 * np.pi / 360) < 10e-2].iloc[0]["YieldStress"]
        sEBT = (s0 + 2 * s45 + s90) / 4
        b_t = sEBT * np.array([1/np.sqrt(2), 1/np.sqrt(2), 0]) / s0

    b_c = - b_t

    P = [b_t, p_ts, p_c, b_c, p_cs, p_t]
    P_name = ["b_t", "p_ts", "p_c", "b_c", "p_cs", "p_t"] 

    mus = np.zeros((6, 3))
    mue = np.zeros((6, 3))

    wm = np.array([r + np.square(np.sin(theta)), r + np.square(np.cos(theta)), - np.sin(theta) * np.cos(theta)])
    wp = np.array([r_s + np.square(np.sin(theta_s)), r_s + np.square(np.cos(theta_s)), np.sin(theta_s) * np.cos(theta_s)])

    #Tangent calculation
    for i in range(len(P)):
        p = P_name[i]
        if str(p) == "b_t":
            n = 1/np.sqrt(1 + np.square(r_tb)) * np.array([1, r_tb, 0])
        elif str(p) == "b_c":
            n = 1/np.sqrt(1 + np.square(r_cb)) * np.array([-1, -r_cb, 0])
        else :
            if str(p) == "p_t":
                dg = np.array([- np.sin(2 * theta), np.sin(2 * theta), np.cos(2 * theta)])
                w = wm
            elif str(p) == "p_ts":
                dg = np.array([np.sin(2 * theta), - np.sin(2 * theta), - np.cos(2 * theta)])
                w = wp
            elif str(p) == "p_c":
                dg = np.array([np.sin(2 * theta), - np.sin(2 * theta), - np.cos(2 * theta)])
                w = wm
            elif str(p) == "p_cs":
                dg = np.array([- np.sin(2 * theta), np.sin(2 * theta), np.cos(2 * theta)])
                w = wp
            
            u = np.cross(w, dg)
            n = u/np.linalg.norm(u)

        u = np.cross(v, n)
        mus[i] = u/np.linalg.norm(u)
    
    mue = np.roll(mus, -1, axis=0)

    #Maximum lambdas calculations
    l_max_s = np.zeros(6)
    l_max_e = np.zeros(6)
    l_max = np.zeros(6)

    for i in range(6):
        halft, halftau = bound_lambda(P[i - 1], P[i], mus[i - 1], mue[i - 1])
        l_max_s[i - 1] = halft
        l_max_e[i - 1] = halftau
    
    for i in range(6):
        l_max[i] = min(l_max_s[i], l_max_e[i-1])
    
    l_max = np.expand_dims(l_max, axis=0)
    return(P, mus, mue, l_max)

def bezier(material, n_pt_curve):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    B_ys, B_r = seg_data_dir(material, 0.6)

    thetas = B_ys[:,0]
    ys = B_ys[:,1]
    rval = B_r[:,1]

    n_curves = len(B_ys)

    P = np.zeros((n_curves, 6, 3))
    mus = np.zeros((n_curves, 6, 3))
    mue = np.zeros((n_curves, 6, 3))
    l_max = np.zeros((n_curves, 6))

    for i in range(n_curves):
        P_new, mus_new, mue_new, l_max_new = param_bezier_section_plane(material, thetas[i], ys[i], ys[n_curves - i - 1], rval[i], rval[n_curves - i - 1])

        P[i] = P_new
        mus[i] = mus_new
        mue[i] = mue_new
        l_max[i] = l_max_new
    
    s = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
    l = np.min(l_max, axis=0) * s
    B = np.zeros((n_curves * n_pt_curve * 6, 3))

    for i in range(n_curves):
        for j in range(6):

            bs = P[i, j - 1]
            be = P[i, j]
            us = mus[i, j - 1]
            ue = mue[i, j - 1]
            ls = l[j - 1]
            le = l[j]
            B_new = bezier_seg(bs, be, us, ue, ls, le, n_pt_curve)
            B[i * 6 * n_pt_curve + j * n_pt_curve : i * 6 * n_pt_curve + (j + 1) * n_pt_curve] = B_new

    return(B)

def bezier_3D(material, n_pt_total):
    n_pt_curve = 200
    B = bezier(material, n_pt_curve)
    
    m = len(B)

    indexes = random.sample(range(m), n_pt_total)
    B_3D = np.zeros((len(B), 6))
    B_3D[:,0] = B[:,0]
    B_3D[:,1] = B[:,1]
    B_3D[:,3] = B[:,2]

    B_3D = B_3D[indexes]

    print(len(B_3D))
    return(B_3D)

"""---------------------------------------------PLOT FUNCTIONS--------------------------------------------------------------------"""
def plot_seg_data_dir(material):
    """
        Plot the directionnal segments.

        Input :
            - material : string
    """
    _, _, _, _, _, _, _, _, thetas, ys_ratio, r_val = param_data_dir(material)
    alphas = [0,0.6,1]
    linestyles = ["solid", "dotted", "dashed"]
    fig, ax = plt.subplots(2)

    for j in range(len(alphas)):
        alpha = alphas[j]
        B_ys, B_r = seg_data_dir(material, alpha)
        
        X = B_ys[:, 0] / (2 * np.pi) * 360
        Y = B_ys[:, 1]
        ax[0].plot(X, Y, color="blue", linewidth=1,linestyle=linestyles[j], label=str(alpha))

        X = B_r[:, 0] / (2 * np.pi) * 360
        Y = B_r[:, 1]
        ax[1].plot(X, Y, color="red", linewidth=1, linestyle=linestyles[j], label=str(alpha))

    legend_lines = []
    for i in range(len(alphas)):
        legend_lines.append(Line2D([0], [0], color='black', linestyle=linestyles[i], label=alphas[i])) 

    ax[0].legend(handles=legend_lines)
    thetas = thetas / (2 * np.pi) * 360
    ax[0].scatter(thetas, ys_ratio, marker="x", linewidths=1, color="blue")
    ax[1].scatter(thetas, r_val, marker="x", linewidths=1, color="red")

    ax[0].set_ylabel(r"$\sigma/\sigma_0[-]$")
    ax[1].set_ylabel(r"$\text{r-value}[-]$")
    ax[1].set_xlabel(r"$\theta[\text{°}]$")

    ax[0].grid(1)
    ax[1].grid(1)

    specific_xtick_labels = ['Zero', 'Two', 'Four']  # Labels for these ticks

    ax[0].set_xticks(thetas, thetas)
    ax[1].set_xticks(thetas, thetas)
    plt.suptitle(f"Bezier model {material} : Yield stresses and Lankford ratios")
    plt.show()

def plot_bezier(material):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    B_ys, B_r = seg_data_dir(material, 0.6)

    thetas = B_ys[:,0]
    ys = B_ys[:,1]
    rval = B_r[:,1]

    n_curves = len(B_ys)
    n_pt_curve = 50

    P = np.zeros((n_curves, 6, 3))
    mus = np.zeros((n_curves, 6, 3))
    mue = np.zeros((n_curves, 6, 3))
    l_max = np.zeros((n_curves, 6))

    for i in range(n_curves):
        P_new, mus_new, mue_new, l_max_new = param_bezier_section_plane(material, thetas[i], ys[i], ys[n_curves - i - 1], rval[i], rval[n_curves - i - 1])

        P[i] = P_new
        mus[i] = mus_new
        mue[i] = mue_new
        l_max[i] = l_max_new

    s = 0
    l = np.min(l_max, axis=0) * s
    B = np.zeros((n_curves * n_pt_curve * 6, 3))

    for i in range(n_curves):
        for j in range(6):

            bs = P[i, j - 1]
            be = P[i, j]
            us = mus[i, j - 1]
            ue = mue[i, j - 1]
            ls = l[j - 1]
            le = l[j]
            B_new = bezier_seg(bs, be, us, ue, ls, le, n_pt_curve)
            B[i * 6 * n_pt_curve + j * n_pt_curve : i * 6 * n_pt_curve + (j + 1) * n_pt_curve] = B_new

    X = B[:,0]
    Y = B[:,1]
    Z = B[:,2]

    ax.scatter(X, Y, Z, linewidth=0.1, s=1)
    ax.set_xlabel(r"$\sigma_{xx}/\sigma_0[-]$")
    ax.set_ylabel(r"$\sigma_{yy}/\sigma_0[-]$")
    ax.set_zlabel(r"$\sigma_{xy}/\sigma_0[-]$")
    plt.suptitle(f"Bezier model {material}, s = {s}")

    plt.show()

def bezier_section_plane_alloutput(material, theta, ys, ys_s, n_pt_curve, r, r_s, r_tb=1, r_cb=1, tc_sym=1):
    df = readdata_exp(material)

    #Intersection of surface Pi(theta) and the yield surface
    theta_s = np.pi/2 - theta
    v = np.array([- np.sin(2 * theta), np.sin(2 * theta), 2 * np.cos(2 * theta)])

    p_t = ys * np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.sin(theta) * np.cos(theta)])
    p_ts = ys_s * np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.sin(theta) * np.cos(theta)])

    p_c = ys * np.array([- np.square(np.cos(theta)), - np.square(np.sin(theta)), - np.sin(theta) * np.cos(theta)])
    p_cs = ys_s * np.array([- np.square(np.sin(theta)), - np.square(np.cos(theta)), np.sin(theta) * np.cos(theta)])

    try:
        b_t = df[df["q"] == 1][["s11", "s22", "s12"]].values[0]
    except :
        #print("EBT non available, we use the approximation : sEBT = (s0 + 2s45 + s90) / 4")
        s0 = df[np.abs(df["LoadAngle"] - 0) < 10e-2].iloc[0]["YieldStress"]
        s45 = df[np.abs(df["LoadAngle"] - 45 * 2 * np.pi / 360) < 10e-2].iloc[0]["YieldStress"]
        s90 = df[np.abs(df["LoadAngle"] - 90 * 2 * np.pi / 360) < 10e-2].iloc[0]["YieldStress"]
        sEBT = (s0 + 2 * s45 + s90) / 4
        b_t = sEBT * np.array([1/np.sqrt(2), 1/np.sqrt(2), 0]) / s0

    b_c = - b_t

    P = [b_t, p_ts, p_c, b_c, p_cs, p_t]
    P_name = ["b_t", "p_ts", "p_c", "b_c", "p_cs", "p_t"]   
    B = np.zeros((n_pt_curve * 6, 3))

    mus = np.zeros((6, 3))
    mue = np.zeros((6, 3))
    ns = np.zeros((6, 3))
    dgs = np.zeros((6, 3))
    ws = np.zeros((6, 3)) 

    wm = np.array([r + np.square(np.sin(theta)), r + np.square(np.cos(theta)), - np.sin(theta) * np.cos(theta)])
    wp = np.array([r_s + np.square(np.sin(theta_s)), r_s + np.square(np.cos(theta_s)), np.sin(theta_s) * np.cos(theta_s)])

    #Tangent calculation
    for i in range(len(P)):
        p = P_name[i]
        if str(p) == "b_t":
            n = 1/np.sqrt(1 + np.square(r_tb)) * np.array([1, r_tb, 0])
            dgs[i] = n
            ws[i] = n
        elif str(p) == "b_c":
            n = 1/np.sqrt(1 + np.square(r_cb)) * np.array([-1, -r_cb, 0])
            dgs[i] = n
            ws[i] = n
        else :
            if str(p) == "p_t":
                dg = np.array([- np.sin(2 * theta), np.sin(2 * theta), np.cos(2 * theta)])
                w = wm
            elif str(p) == "p_ts":
                dg = np.array([np.sin(2 * theta), - np.sin(2 * theta), - np.cos(2 * theta)])
                w = wp
            elif str(p) == "p_c":
                dg = np.array([np.sin(2 * theta), - np.sin(2 * theta), - np.cos(2 * theta)])
                w = wm
            elif str(p) == "p_cs":
                dg = np.array([- np.sin(2 * theta), np.sin(2 * theta), np.cos(2 * theta)])
                w = wp
            
            u = np.cross(w, dg)
            n = u/np.linalg.norm(u)
            dgs[i] = dg
            ws[i] = w

        ns[i] = n
        u = np.cross(v, n)
        mus[i] = u/np.linalg.norm(u)
    
    mue = np.roll(mus, -1, axis=0)

    #Maximum lambdas calculations
    l_max = np.zeros(6)
    for i in range(6):
        halft, halftau = bound_lambda(P[i - 1], P[i], mus[i - 1], mue[i - 1])
        l_max[i - 1] = min(halft, halftau)
    
    print(l_max)
    #Orthotrophy hypothesis
    L1 = l_max[0]
    L2 = min(l_max[1], l_max[5])
    L3 = min(l_max[2], l_max[4])
    L4 = l_max[3]

    if tc_sym:
        L2 = min(L2, L3)
        L3 = L2

    ls = np.array([L1, L2, L3, L4, L3, L2])

    for i in range(6):
        B[i * n_pt_curve:(i+1) * n_pt_curve] = bezier_seg(P[i - 1], P[i], mus[i - 1], mue[i - 1], ls[i - 1], ls[i], n_pt_curve)
    
    return(B, P, mus, mue, ns, dgs, ws, v)

def plot_bezier_alloutput(material):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    B_ys, B_r = seg_data_dir(material, 0.6)

    thetas = B_ys[:,0]
    ys = B_ys[:,1]
    rval = B_r[:,1]

    n_data = len(B_ys)
    n_pt_curve = 100

    B = np.empty((0,3))
    P = np.empty((0,3))
    mus = np.empty((0,3))
    mue = np.empty((0,3))
    ns = np.empty((0,3))
    dgs = np.empty((0,3))
    ws = np.empty((0,3))
    v = np.empty((0,3))

    for i in range(n_data):

        B_new, P_new, mus_new, mue_new, ns_new, dgs_new, ws_new, v_new = bezier_section_plane_alloutput(material, thetas[i], ys[i], ys[n_data - i - 1], n_pt_curve, rval[i], rval[n_data - i - 1])
        B = np.concatenate((B, B_new))
        P = np.concatenate((P, P_new))
        mus = np.concatenate((mus, mus_new))
        mue = np.concatenate((mue, mue_new))
        ns = np.concatenate((ns, ns_new))
        dgs = np.concatenate((dgs, dgs_new))
        ws = np.concatenate((ws, ws_new))
        v_new = np.expand_dims(v_new, axis=0)
        v = np.concatenate((v, v_new))

    X = B[:,0]
    Y = B[:,1]
    Z = B[:,2]

    mus = P + mus  
    ns = P + ns
    dgs = P + dgs
    ws = P + ws
    k = 1
    ax.scatter(X[k * 600:(k+1) * 600], Y[k * 600:(k+1) * 600], Z[k * 600:(k+1) * 600], cmap='viridis', marker="o",s=1, linewidth=1)

    for i in range(k * 6, (k+1) * 6):
        x = [P[i, 0], mus[i, 0]]
        y = [P[i, 1], mus[i, 1]]
        z = [P[i, 2], mus[i, 2]]
        ax.plot(x, y, z, marker='o', markersize=1)

    for i in range(k * 6, (k+1) * 6):
        x = [P[i, 0], ns[i, 0]]
        y = [P[i, 1], ns[i, 1]]
        z = [P[i, 2], ns[i, 2]]
        #ax.plot(x, y, z, marker='o', color="green")

    for i in range(k * 6, (k+1) * 6):
        x = [P[i, 0], ws[i, 0]]
        y = [P[i, 1], ws[i, 1]]
        z = [P[i, 2], ws[i, 2]]
        #ax.plot(x, y, z, marker='o', color="red")

    for i in range(k * 6, (k+1) * 6):
        x = [P[i, 0], dgs[i, 0]]
        y = [P[i, 1], dgs[i, 1]]
        z = [P[i, 2], dgs[i, 2]]
        #ax.plot(x, y, z, marker='o', color="blue")

    Pe = np.roll(P, -1, axis=0)
    mue = Pe + mue

    for i in range(k * 6, (k+1) * 6):
        x = [Pe[i, 0], mue[i, 0]]
        y = [Pe[i, 1], mue[i, 1]]
        z = [Pe[i, 2], mue[i, 2]]
        #ax.plot(x, y, z, marker='o')
    normal_vector = v[k]
    point = np.array([0, 0, 0])

    # Create the figure and 3D axis

    # Define the plane equation ax + by + cz = d
    d = np.dot(normal_vector, point)

    # Create a grid of x, y values
    xx, yy = np.meshgrid(range(-2, 2), range(-2, 2))

    # Calculate the corresponding z values
    a, b, c = normal_vector[0], normal_vector[1], normal_vector[2]
    zz = (d - a * xx - b * yy) / c

    # Plot the plane
    ax.plot_surface(xx, yy, zz, alpha=0.5, rstride=100, cstride=100)

    plt.show()

if __name__ == "__main__":
    plot_seg_data_dir("DP600")