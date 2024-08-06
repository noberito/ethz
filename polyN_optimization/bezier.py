

from ctypes.wintypes import WPARAM
import pstats
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import binom

from optimize_polyN_mini import readData_2d
from read_param import read_param


c50 = binom(5,0)
c51 = binom(5,1)
c52 = binom(5,2)
c53 = binom(5,3)
c54 = binom(5,4)
c55 = binom(5,5)

print(c51)
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

    df = readData_2d(material, "bezier")
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
    YS, mus_ys, mue_ys, lmax_ys, R, mus_r, mue_r, lmax_r, thetas, ys_ratio, r_val = param_data_dir(material)
    l_ys = alpha * lmax_ys
    l_r = alpha * lmax_r
    n_seg = len(mus_ys)
    n = 10

    B_ys = np.zeros((n * n_seg, 2))
    B_r = np.zeros((n * n_seg, 2))

    for i in range(n_seg):

        B_ys[(i * n):((i + 1) * n)] = bezier_seg(YS[i], YS[i + 1], mus_ys[i], mue_ys[i], l_ys, l_ys, n)
        B_r[(i * n):((i + 1) * n)] = bezier_seg(R[i], R[i + 1], mus_r[i], mue_r[i], l_r, l_r, n)

    return(B_ys, B_r)

def plot_seg_data_dir(material):
    YS, mus_ys, mue_ys, lmax_ys, R, mus_r, mue_r, lmax_r, thetas, ys_ratio, r_val = param_data_dir(material)
    alphas = [0,0.6,1]
    linestyles = ["solid", "dotted", "dashed"]
    fig, ax = plt.subplots(2)

    for j in range(len(alphas)):
        alpha = alphas[j]
        _, B_ys, B_r = seg_data_dir(material, alpha)
        
        X = B_ys[:, 0]
        Y = B_ys[:, 1]
        ax[0].plot(X, Y, color="blue", linewidth=1,linestyle=linestyles[j], label=str(alpha))

        X = B_r[:, 0]
        Y = B_r[:, 1]
        ax[1].plot(X, Y, color="red", linewidth=1, linestyle=linestyles[j], label=str(alpha))

    legend_lines = []
    for i in range(len(alphas)):
        legend_lines.append(Line2D([0], [0], color='black', linestyle=linestyles[i], label=alphas[i])) 

    ax[0].legend(handles=legend_lines)

    ax[0].scatter(thetas, ys_ratio, marker="x", linewidths=1, color="blue")
    ax[1].scatter(thetas, r_val, marker="x", linewidths=1, color="red")
    ax[0].grid(1)
    ax[1].grid(1)
    plt.suptitle(f"Bezier model {material} : Yield stresses and Lankford ratios")
    plt.show()

def bound_lambda(bs, be, us, ue):
    v = np.dot(us, ue)
    u =(be - bs) /(1 - v) 

    t = np.dot(u, us - v * ue)
    tau = np.dot(u, ue - v * us)

    return(t/2, tau/2)


def bezier_section_plane(material, theta, ys, ys_s, n_pt, r, r_s, r_tb=1, r_cb=1, tc_sym=1):
    theta_s = np.pi/2 - theta
    df = readData_2d(material, "bezier")
    v = np.array([- np.sin(2 * theta), np.sin(2 * theta), np.cos(2 * theta)])

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

    p_t = ys * np.array([np.square(np.cos(theta)), np.square(np.sin(theta)), np.sin(theta) * np.cos(theta)])
    print(np.dot(p_t, v))
    p_ts = ys_s * np.array([np.square(np.sin(theta)), np.square(np.cos(theta)), - np.sin(theta) * np.cos(theta)])

    p_c = ys * np.array([- np.square(np.cos(theta)), - np.square(np.sin(theta)), - np.sin(theta) * np.cos(theta)])
    p_cs = ys_s * np.array([- np.square(np.sin(theta)), - np.square(np.cos(theta)), np.sin(theta) * np.cos(theta)])

    P = [b_t, p_ts, p_c, b_c, p_cs, p_t]
    P_name = ["b_t", "p_ts", "p_c", "b_c", "p_cs", "p_t"]   
    B = np.zeros((n_pt * 6, 3))

    mus = np.zeros((6, 3))
    mue = np.zeros((6, 3))
    ns = np.zeros((6, 3))
    dgs = np.zeros((6, 3))
    ws = np.zeros((6, 3)) 

    wm = np.array([r + np.square(np.sin(theta)), r + np.square(np.cos(theta)), - np.sin(theta) * np.cos(theta)])
    wp = np.array([r_s + np.square(np.sin(theta_s)), r_s + np.square(np.cos(theta_s)), np.sin(theta_s) * np.cos(theta_s)])

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
    
    print(np.dot(mus, v))
    mue = np.roll(mus, -1, axis=0)

    l_max = np.zeros(6)
    for i in range(6):
        halft, halftau = bound_lambda(P[i - 1], P[i], mus[i - 1], mue[i - 1])
        l_max[i - 1] = min(halft, halftau)
    
    L1 = l_max[0]
    L2 = min(l_max[1], l_max[5])
    L3 = min(l_max[2], l_max[4])
    L4 = l_max[3]

    if tc_sym:
        L2 = min(L2, L3)
        L3 = L2

    ls = np.array([L1, L2, L3, L4, L3, L2]) * 0.5

    for i in range(6):
        B[i * n_pt:(i+1) * n_pt] = bezier_seg(P[i - 1], P[i], mus[i - 1], mue[i - 1], ls[i - 1], ls[i], n_pt)
    
    return(B, P, mus, mue, ns, dgs, ws, v)


def plot_bezier(material):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    B_ys, B_r = seg_data_dir(material, 0.6)

    thetas = B_ys[:,0]
    ys = B_ys[:,1]
    rval = B_r[:,1]

    n_data = len(B_ys)
    n_pt = 100

    B = np.empty((0,3))

    for i in range(n_data):
        B_new, P_new, mus_new, mue_new, ns_new, dgs_new, ws_new = bezier_section_plane(material, thetas[i], ys[i], ys[n_data - i - 1], n_pt, rval[i], rval[n_data - i - 1])
        B = np.concatenate((B, B_new))

    X = B[:,0]
    Y = B[:,1]
    Z = B[:,2]

    ax.scatter(X, Y, Z, linewidth=0.01, s=1)

    plt.show()

def plot_bezier_tangent(material):
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    B_ys, B_r = seg_data_dir(material, 0.6)

    thetas = B_ys[:,0]
    ys = B_ys[:,1]
    rval = B_r[:,1]

    n_data = len(B_ys)
    n_pt = 100

    B = np.empty((0,3))
    P = np.empty((0,3))
    mus = np.empty((0,3))
    mue = np.empty((0,3))
    ns = np.empty((0,3))
    dgs = np.empty((0,3))
    ws = np.empty((0,3))
    v = np.empty((0,3))

    for i in range(n_data):

        B_new, P_new, mus_new, mue_new, ns_new, dgs_new, ws_new, v_new = bezier_section_plane(material, thetas[i], ys[i], ys[n_data - i - 1], n_pt, rval[i], rval[n_data - i - 1])
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

material = "DP780"
plot_bezier_tangent(material)