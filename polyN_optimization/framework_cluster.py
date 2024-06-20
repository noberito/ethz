import subprocess
import os
import glob
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np
import time
from scipy.interpolate import interp1d
from scipy.integrate import simpson
from optimize_polyN_cluster import first_opti, write_coeff_abq, get_param_polyN
from get_calibration_data import analyze_exp_data

sep = os.sep
polyN_cali_dir = os.path.dirname(os.path.abspath(__file__))
exec_dir = polyN_cali_dir + sep + "running"


from data_visualizing.compare_sim_exp import compare_large_strain
from read_param import read_param
from running.run_large_sim_cluster import run

p = read_param()
material = p["material"]
law = p["law"]
degree = int(p["degree"])
protomodel = p["protomodel"]
input_type = p["input_type"]
n_opti = int(p["n_opti"])
density = p["density"]
enu = p["enu"]
var_optim = int(p["var_optim"])
results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

mat_exp = analyze_exp_data(material)
tests = []
for type_test in mat_exp.keys():
    for ori in mat_exp[type_test]:
        if type_test != "UT":
            test = type_test + "_" + ori
            tests.append(test)

def numerical_gradient(f, x):
    h = 1e-5
    grad = np.zeros_like(x)
    for i in range(x.size):
        temp_val = x[i]
        x[i] = temp_val + h
        fxh1 = f(x)
        
        x[i] = temp_val - h 
        fxh2 = f(x)  
        
        grad[i] = (fxh1 - fxh2) / (2 * h)
        x[i] = temp_val  
        
    return grad

def gradient_descent(f, x0, learning_rate=0.1, n_opti=n_opti):
    x = x0
    for i in range(n_opti):
        grad = numerical_gradient(f, x)
        x = x - learning_rate * grad
    return x

def mean_square_error(test):
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
    sim_res_path = results_sim_dir + sep + test + "_" + input_type + ".csv"
    df_sim = pd.read_csv(sim_res_path)
    x_sim = df_sim["U2"]
    y_sim = df_sim["RF2"]
    f_sim = interp1d(x_sim, y_sim, kind="linear", fill_value="extrapolate")

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    type_test, ori = test.split("_")
    n_exp = mat_exp[type_test][ori]
    dfs = []
    lens = []
    disps = np.empty(0)
    forces = np.empty(0)
    for m in range(1):
        exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{m+1}.csv"
        df_exp = pd.read_csv(exp_res_path)
        e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
        s = df_exp["Force[kN]"]
        lens.append(len(e))
        disps = np.concatenate((disps, e))
        forces = np.concatenate((forces, s))
    
    f_exp = interp1d(disps, forces, kind="linear", fill_value="extrapolate")
    x_common = np.linspace(0, min(max(x_sim), max(e)), num=500)
    y1_common = f_sim(x_common)
    y2_common = f_exp(x_common)
    y_diff = np.abs(y1_common - y2_common)
    area_between_curves = simpson(y=y_diff, x=x_common)
    return(area_between_curves)

def framework():

    t0 = time.time()
    coeff_polyN = first_opti()
    powers = get_param_polyN(degree)
    nmon = len(powers)

    run(tests, material, degree, law, protomodel, input_type)
    compare_large_strain(material, degree, input_type)
    coeff_law = np.load(polyN_cali_dir + sep + f"{law}_coeff.npy")
    
    a = coeff_law[0]
    b = coeff_law[1]
    c = coeff_law[2]
    ymod = coeff_law[3]

    dcoeffs = np.zeros((n_opti, 22))
    dcoeffs[:,8] = np.arange(-5, -5 + n_opti)

    for dcoeff in dcoeffs :
        new_coeff = coeff_polyN + dcoeff
        write_coeff_abq(new_coeff, a, b, c, ymod, enu, nmon, protomodel, degree, material, law, density, powers)

        run(tests, material, degree, law, protomodel, input_type)
        compare_large_strain(material, degree, input_type)

def framework_2():

    t0 = time.time()
    coeff_polyN = first_opti()
    
    coeff_law = np.load(polyN_cali_dir + sep + f"{law}_coeff.npy")
    a = coeff_law[0]
    b = coeff_law[1]
    c = coeff_law[2]
    ymod = coeff_law[3]
    powers = get_param_polyN(degree)
    nmon = len(powers)

    def f_cost(x):
        powers = get_param_polyN(degree)
        nmon = len(powers)

        new_coeff = coeff_polyN
        new_coeff[var_optim] = new_coeff[var_optim] + x

        write_coeff_abq(new_coeff, a, b, c, ymod, enu, nmon, protomodel, degree, material, law, density, powers)
        run(tests, material, degree, law, protomodel, input_type)

        err = 0

        for test in tests:
            err = err + mean_square_error(test)
        
        return(err)

    result = gradient_descent(f_cost, [0])

    coeff_polyN[var_optim] = coeff_polyN[var_optim] + result
    coeff_file = polyN_cali_dir + sep + "coeff.npy"
    np.save(coeff_polyN, coeff_file)
    print(coeff_polyN)