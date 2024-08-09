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


from optimize_polyN_mini import firstopti_mini, write_coeff_abq_mini, get_param_polyN_mini, write_coeff_abq_mini
from get_calibration_data import analyze_exp_data
from tests_parameters import ut_tests_ext
from check.compare_sim_exp import compare_large_strain, compare_ut_s_2
from read import read_param, get_coeff_mini, get_coeff_law
from run_large_sim_cluster import run, post_process

sep = os.sep
polyN_dir = os.path.dirname(os.path.abspath(__file__))
exec_dir = polyN_dir + sep + "running"

def mean_square_error_fd(material, test, input_type, p=0, m=0):
    mat_exp = analyze_exp_data(material)
    results_sim_dir = polyN_dir + sep + "results_sim" + sep + material
    sim_res_path = results_sim_dir + sep + test + "_" + input_type + "_" + str(p) + "_" + str(m) + ".csv"
    df_sim = pd.read_csv(sim_res_path)
    x_sim = df_sim["U2"]
    y_sim = df_sim["RF2"]
    f_sim = interp1d(x_sim, y_sim, kind="linear", fill_value="extrapolate")

    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material
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

def mean_square_error_str(material, test, input_type, p=0, m=0):
    mat_exp = analyze_exp_data(material)
    results_sim_dir = polyN_dir + sep + "results_sim" + sep + material
    sim_res_path = results_sim_dir + sep + test + "_" + input_type + "_" + str(p) + "_" + str(m) + ".csv"
    df_sim = pd.read_csv(sim_res_path)

    results_exp_dir = polyN_dir + sep + "results_exp" + sep + material
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

        for k in range(n_exp):
            exp_res_path = results_exp_dir + sep + type_test + "_" + ori + f"_{k+1}.csv"
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

def framework_mini(material, degree, law, enu, protomodel, input_type, density, tests, var_optim):
    """
        Framework to calibrate the polyN minimalistic yield function and
        large strain optimization on variables in var_optim using scipy

        Input :
            - var_optim : int or list

    """
    t0 = time.time()

    #Loading coefficients

    try :
        coeff_mini = get_coeff_mini(material, degree)  
    except:
        coeff_mini = firstopti_mini()
    coeff_law, ymod = get_coeff_law(material, law)


    p = str(var_optim[0])
    for var in var_optim[1:]:
        p = p + "0" + str(var)

    def f_cost(x, m=1000):
        powers = get_param_polyN_mini(degree)

        new_coeff = np.copy(coeff_mini)
        new_coeff[var_optim] = new_coeff[var_optim] + x

        write_coeff_abq_mini(new_coeff, coeff_law, ymod, enu, protomodel, degree, material, law, density, p=p, m=m)
        run(tests, material, degree, law, protomodel, input_type, p=p, m=m)
        run(ut_tests_ext, material, degree, law, protomodel, input_type, p=p, m=m)
        post_process(tests, material, input_type,p=p, m=m)
        post_process(ut_tests_ext, material, input_type,p=p, m=m)
        compare_large_strain(material, degree, input_type,p=p, m=m)
        compare_ut_s_2(material, degree, input_type, p=p, m=m)

        err = 0

        for test in tests:
            err = err + mean_square_error_fd(material, test, input_type, p=p, m=m) + mean_square_error_str(material, test, input_type, p=p, m=m)
            print("err:", err)
        
        return(err)
    
    x0 = np.zeros(len(var_optim))
    b = 0.1 * np.abs(coeff_mini)
    b = b[var_optim]
    bounds = scipy.optimize.Bounds(lb = -b, ub = b, keep_feasible=True)
    result = scipy.optimize.minimize(f_cost, x0, method="SLSQP", jac="3-point", tol=10e-12, bounds=bounds)

    coeff_mini[var_optim] = coeff_mini[var_optim] + result.x

    filedir = polyN_dir + sep + "coeff" + sep
    coeff_filename = material + "_poly" + str(degree) + "_mini_" + str(p) + ".npy"
    coeff_file = filedir + coeff_filename

    np.save(coeff_file, coeff_mini)
    print(coeff_mini)


if __name__ == "__main__":
    p = read_param()

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
    results_sim_dir = polyN_dir + sep + "results_sim" + sep + material

    mat_exp = analyze_exp_data(material)

    tests = []
    for type_test in mat_exp.keys():
        for ori in mat_exp[type_test]:
            if type_test != "UT":
                test = type_test + "_" + ori
                tests.append(test)
            
   
    framework_mini(material, degree, law, enu, protomodel, input_type, density, tests, var_optim)