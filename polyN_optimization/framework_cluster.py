import subprocess
import os
import glob
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np
import time
import scipy.optimize
from scipy.integrate import simps
from scipy.interpolate import interp1d


from optimize_polyN import first_opti, get_param_polyN, write_coeff_abq
from optimize_polyN_mini import firstopti_mini, write_coeff_abq_mini, get_param_polyN_mini, write_coeff_abq_mini
from get_calibration_data import analyze_exp_data
from tests_parameters import ut_tests_ext
from check.compare_sim_exp import compare_large_strain, compare_ut_s_2, compare_ut_fd
from read_param import read_param
from run_large_sim_cluster import launch_run, create_csv

sep = os.sep
polyN_cali_dir = os.path.dirname(os.path.abspath(__file__))
exec_dir = polyN_cali_dir + sep + "running"


p = read_param()
func = p["func"]
material = p["material"]
law = p["law"]
degree = int(p["degree"])
protomodel = p["protomodel"]
input_type = p["input_type"]
n_opti = int(p["n_opti"])
density = p["density"]
var_optim = p["var_optim"].split(",")
var_optim = np.array([int(var_optim[i]) for i in range(len(var_optim))])
enu = p["enu"]
results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

mat_exp = analyze_exp_data(material)
tests = []
for type_test in mat_exp.keys():
    for ori in mat_exp[type_test]:
        if type_test != "UT":
            test = type_test + "_" + ori
            tests.append(test)

def numerical_gradient(f, x, m):
    h = 0.05
    grad = np.zeros_like(x)
    for i in range(x.size):

        print(m)
        temp_val = x[i]
        x[i] = temp_val + h
        fxh1 = f(x, 100 + 2 * m)
        
        x[i] = temp_val - h 
        fxh2 = f(x, 100 + 2 * m + 1)  
        
        grad[i] = (fxh1 - fxh2) / (2 * h)
        x[i] = temp_val  
        
    return grad

def gradient_descent(f, x0, learning_rate=0.1, n_opti=n_opti):
    x = x0
    for i in range(n_opti):
        grad = numerical_gradient(f, x, i)
        print("grad, lr", grad, learning_rate)
        x = x - learning_rate * grad
        f(x, m=i)
    return x

def mean_square_error_fd(test, p=0, m=0):
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
    sim_res_path = results_sim_dir + sep + test + "_" + input_type + "_" + str(var_optim) + "_" + str(m) + ".csv"
    df_sim = pd.read_csv(sim_res_path)
    x_sim = df_sim["U2"]
    y_sim = df_sim["RF2"]
    f_sim = interp1d(x_sim, y_sim, kind="linear", fill_value="extrapolate")

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    type_test, ori = test.split("_")
    n_exp = mat_exp[type_test][ori]
    lens = []
    disps = np.empty(0)
    forces = np.empty(0)

    area_between_curves = 0

    for m in range(n_exp):
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
        y1_common[0] = 0
        y2_common[0] = 0
        y_diff = np.abs(y1_common - y2_common)
        area_between_curves = area_between_curves + simps(y=y_diff, x=x_common)

    return(area_between_curves)

def mean_square_error_str(test, p=0, m=0):
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
    sim_res_path = results_sim_dir + sep + test + "_" + input_type + "_" + str(var_optim) + "_" + str(m) + ".csv"
    df_sim = pd.read_csv(sim_res_path)

    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    type_test, ori = test.split("_")
    n_exp = mat_exp[type_test][ori]
    lens = []
    disps = np.empty(0)
    strains = np.empty(0)

    area_between_curves = 0

    if "Strain_ext" in df_sim.columns:
        x_sim = df_sim["U2"]
        y_sim = df_sim["Strain_ext"]
        f_sim = interp1d(x_sim, y_sim, kind="linear", fill_value="extrapolate")

        for m in range(n_exp):
            exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{m+1}.csv"
            df_exp = pd.read_csv(exp_res_path)
            e = df_exp["Displacement longi[mm]"] if type_test == "SH" else df_exp["Displacement[mm]"]
            s = df_exp["AxStrain_1"]
            lens.append(len(e))
            disps = np.concatenate((disps, e))
            strains = np.concatenate((strains, s))
            
            f_exp = interp1d(disps, strains, kind="linear", fill_value="extrapolate")
            x_common = np.linspace(0, min(max(x_sim), max(e)), num=500)
            y1_common = f_sim(x_common)
            y2_common = f_exp(x_common)
            y1_common[0] = 0
            y2_common[0] = 0
            y_diff = np.abs(y1_common - y2_common)
            area_between_curves = area_between_curves + simps(y=y_diff, x=x_common)
        
    return(area_between_curves)

def framework(var_optim):
    """
        Framework to calibrate the  yield function and
        large strain optimization on variables in var_optim using using gradient 
        descent coded above

        Input :
            - var_optim : int or list

    """
    t0 = time.time()
    if(1):
        coeff_polyN = np.load(polyN_cali_dir + sep + "polyN_coeff.npy")
    else:
        coeff_polyN = first_opti()
    coeff_law = np.load(polyN_cali_dir + sep + f"{law}_coeff.npy")
    a = coeff_law[0]
    b = coeff_law[1]
    c = coeff_law[2]
    ymod = coeff_law[3]

    def f_cost(x, m):
        powers = get_param_polyN(degree)
        nmon = len(powers)
        new_coeff = np.copy(coeff_polyN)
        new_coeff[var_optim - 1] = new_coeff[var_optim - 1] + x
        write_coeff_abq(new_coeff, a, b, c, ymod, enu, protomodel, degree, material, law, density, powers, m=m)
        time.sleep(5)
        launch_run(tests, func, material, degree, law, protomodel, input_type, m=m)
        create_csv(tests, material, input_type, m=m)
        compare_ut_fd(material, degree, input_type, m=m)

        err = 0

        for test in tests:
            err = err + mean_square_error_fd(test, m=m)
            print("err:", err)
        
        return(err)
    

    x0 = np.zeros(len(var_optim))
    result = gradient_descent(f_cost, x0)

    coeff_polyN[var_optim] = coeff_polyN[var_optim] + result
    coeff_file = polyN_cali_dir + sep + "polyN_coeff.npy"
    np.save(coeff_polyN, coeff_file)
    print(coeff_polyN)


def framework_mini_old(var_optim):
    """
        Framework to calibrate the polyN minimalistic yield function and
        large strain optimization on variables in var_optim using gradient 
        descent coded above

        Input :
            - var_optim : int or list

    """
    t0 = time.time()

    coeff_file = polyN_cali_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy"
    
    if os.path.exists(coeff_file):
        coeff_polyN_mini = np.load(coeff_file)
    else :
        coeff_polyN_mini = firstopti_mini()

    coeff_law = np.load(polyN_cali_dir + sep + f"{material}_{law}_mini_coeff.npy")
    ymod = coeff_law[-1]
    coeff_law = coeff_law[:-1]

    p = str(var_optim[0])

    for var in var_optim[1:]:
        p = "0" + str(var)

    def f_cost(x, m):

        powers = get_param_polyN_mini(degree)
        nmon = len(powers)

        new_coeff = np.copy(coeff_polyN_mini)
        new_coeff[var_optim] = new_coeff[var_optim] + x

        write_coeff_abq_mini(new_coeff, coeff_law, ymod, enu, protomodel, degree, material, law, density, powers, p=p, m=m)
        launch_run(tests, func, material, degree, law, protomodel, input_type, p=p, m=m)
        launch_run(ut_tests_ext, func, material, degree, law, protomodel, input_type, p=p, m=m)
        create_csv(tests, material, input_type,p=p, m=m)
        create_csv(ut_tests_ext, material, input_type,p=p, m=m)
        compare_large_strain(material, func, degree, input_type,p=p, m=m)
        compare_ut_s_2(material, func, degree, input_type, p=p, m=m)

        err = 0

        for test in tests:
            err = err + mean_square_error_fd(test,p=var_optim, m=m) + mean_square_error_str(test,p=var_optim, m=m)
            print("err:", err)
        
        return(err)
    
    x0 = np.zeros(len(var_optim))
    result = gradient_descent(f_cost, x0)

    coeff_polyN_mini[var_optim] = coeff_polyN_mini[var_optim] + result
    coeff_file = polyN_cali_dir + sep + material + "_poly" + str(degree) + "_mini_coeff_final_ " + str(var_optim) + ".npy"
    np.save(coeff_file, coeff_polyN_mini)
    print(coeff_polyN_mini)

def framework_mini(var_optim):
    """
        Framework to calibrate the polyN minimalistic yield function and
        large strain optimization on variables in var_optim using scipy

        Input :
            - var_optim : int or list

    """
    t0 = time.time()

    coeff_file = polyN_cali_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy"
    
    if os.path.exists(coeff_file):
        coeff_polyN_mini = np.load(coeff_file)
    else :
        coeff_polyN_mini = firstopti_mini()

    b = 0.1 * np.abs(coeff_polyN_mini)
    b = b[var_optim]

    coeff_law = np.load(polyN_cali_dir + sep + f"{material}_{law}_mini_coeff.npy")
    ymod = coeff_law[-1]
    coeff_law = coeff_law[:-1]
    p = str(var_optim[0])

    for var in var_optim[1:]:
        p = p + "0" + str(var)

    def f_cost(x, m=1000):

        powers = get_param_polyN_mini(degree)
        nmon = len(powers)

        new_coeff = np.copy(coeff_polyN_mini)
        new_coeff[var_optim] = new_coeff[var_optim] + x

        write_coeff_abq_mini(new_coeff, coeff_law, ymod, enu, protomodel, degree, material, law, density, powers, p=p, m=m)
        launch_run(tests, func, material, degree, law, protomodel, input_type, p=p, m=m)
        launch_run(ut_tests_ext, func, material, degree, law, protomodel, input_type, p=p, m=m)
        create_csv(tests, material, input_type,p=p, m=m)
        create_csv(ut_tests_ext, material, input_type,p=p, m=m)
        compare_large_strain(material, func, degree, input_type,p=p, m=m)
        compare_ut_s_2(material, func, degree, input_type, p=p, m=m)

        err = 0

        for test in tests:
            err = err + mean_square_error_fd(test,p=p, m=m) + mean_square_error_str(test,p=p, m=m)
            print("err:", err)
        
        return(err)
    
    x0 = np.zeros(len(var_optim))
    bounds = scipy.optimize.Bounds(lb = -b, ub = b, keep_feasible=True)
    result = scipy.optimize.minimize(f_cost, x0, method="SLSQP", jac="3-point", tol=10e-12, bounds=bounds)

    coeff_polyN_mini[var_optim] = coeff_polyN_mini[var_optim] + result.x
    coeff_file = polyN_cali_dir + sep + material + "_poly" + str(degree) + "_mini_coeff_final_scipy_" + str(var_optim) + ".npy"
    np.save(coeff_file, coeff_polyN_mini)
    print(coeff_polyN_mini)


if __name__ == "__main__":
    if func == "polyN":
        framework(var_optim)
    elif func == "polyN_mini":
        framework_mini(var_optim)