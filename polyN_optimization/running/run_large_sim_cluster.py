import subprocess
import os
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np


# Get the directory where the Python exec is located

sep = os.sep


exec_dir = os.path.dirname(os.path.abspath(__file__))
polyN_cali_dir = os.path.dirname(exec_dir)

sys.path.append(polyN_cali_dir)

from read_param import read_param
from get_calibration_data import analyze_exp_data

p = read_param()

material = p["material"]
law = p["law"]
degree = int(p["degree"])
protomodel = p["protomodel"]
input_type = p["input_type"]

results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material

dt = 1

def change_paras(mat_tests):
    for type_test in mat_tests.keys():
        if type_test != "UT":
            for ori in mat_tests[type_test]:
                test = type_test + "_" + ori
                res_filename = results_exp_dir + sep + test + "_1.csv"
                df_exp = pd.read_csv(res_filename, index_col=False)
                dtime = np.max(df_exp.iloc[:,0])
                displ = np.max(df_exp.iloc[:,1])
                thickness = df_exp.loc[0, "Thickness[mm]"]
                width = df_exp.loc[0, "OuterWidth[mm]"]
                tramp1 = dtime / 100
                ramp1 = displ / 100
                dt = dtime / 10000

                para_filename = exec_dir + sep + type_test + "_" + ori + "_paras.inp"
                with open(para_filename, "w") as f:
                    f.write("*PARAMETER\n")
                    f.write(f"DISPL = {displ} * 10 \n")
                    f.write(f"DTIME = {dtime}\n")
                    f.write(f"THICKNESS = {thickness}\n")
                    f.write(f"WIDTH = {width}\n")
                    f.write(f"TRAMP1 = {tramp1}\n")
                    f.write(f"RAMP1 = {ramp1}\n")
                    f.write(f"DT = {dt}\n")
                    f.write("MINDT = DTIME * 1e-6\n")
                    f.write("MAXDT = DTIME * 1e-2")
                    f.close()

def change_usermat(usermatfile):
    """
        Change the input user material file in the all the abaqus files of test in tests
        Input :
            - usermatfile : string, filename of the user material file
    """
    for test in tests:

        filename = exec_dir + sep + test + f"_{input_type}.inp"
        with open(filename, "r") as f:
            content = f.readlines()
        i = 0
        line = content[i]
        vu = 0
        
        while not vu :
            if len(line) > 9 and line[:9] == "*MATERIAL":
                vu = 1
            i = i + 1
            line = content[i]
        content[i] = f"*INCLUDE, INPUT={usermatfile}\n"
        with open(filename, "w") as f:
            f.writelines(content)

def simulate(test):
    """
        Run the abaqus simulation of the test_polyN.inp and generate a csv file in results_sim folder
        Input :
            - test : string, name of the test (ex : UT_00)
    """
    print(f"Simulating {test}")
    job = f'{test}_{input_type}'
    subroutine = f"{input_type}_PolyN_3D_{law}.for"

    input_file = f"{job}.inp"
    cp_input_file = f'temp_{input_file}'
    odb = "{}.odb".format(job)

    copy_sim_cmd = f'cp {input_file} {cp_input_file} '
    subprocess.call(copy_sim_cmd, shell=True, cwd=exec_dir)

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    job_command = f'sbatch -n 24 -t 0-2 --mem-per-cpu=300 --tmp=300 --wrap "abaqus job=temp_{job} input={input_file} double interactive user={subroutine} cpus=24 scratch=\$TMPDIR"'
    
    subprocess.call(job_command, shell=True, cwd=exec_dir)

if __name__ == "__main__":
    mat_tests = analyze_exp_data(material)
    tests = []

    for type_test in mat_tests.keys():
        for ori in mat_tests[type_test]:
            if type_test != "UT":
                test = type_test + "_" + ori
                tests.append(test)
    
    if(0):
        tests = ["CH_00", "CH_45",
         "NT6_00", "NT6_45", "NT6_90",
         "NT20_00", "NT20_45",
         "SH_000", "SH_090", "SH_p45"]

    usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}.inp"
    change_paras(mat_tests)
    change_usermat(usermatfile)
    for test in tests:
        simulate(test)
