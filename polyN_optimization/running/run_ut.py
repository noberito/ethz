import subprocess
import os
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np


# Get the directory where the Python exec is located

sep = os.sep
run_dir = os.path.dirname(os.path.abspath(__file__))
polyN_cali_dir = os.path.dirname(run_dir)

sys.path.append(polyN_cali_dir)
from tests_parameters import load_points
from read_param import read_param

ut_tests = ["UT_00", "UT_15", "UT_30", "UT_45", "UT_60", "UT_75","UT_90", "UT_EBT"]

tests = ["UT_00", "UT_15", "UT_30", "UT_45", "UT_60", "UT_75","UT_90", "UT_EBT"]
""""CH_00", "CH_45",
"NT6_00", "NT6_45", "NT6_90",
"NT20_00", "NT20_45",
"SH_000", "SH_090", "SH_p45"]"""

dt = 1

def change_usermat(test, func, material, degree, law, protomodel, input_type, var_optim=0, n_try=0):
    """
        Change the input user material file in the all the abaqus files of test in tests
        Input :
            - usermatfile : string, filename of the user material file
    """
    if func=="polyN":
        if var_optim==0 and n_try==0:
            usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}.inp"
        else:
            usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}_{var_optim}_{n_try}.inp"
    if func=="polyN_mini":
        if var_optim==0 and n_try==0:
            usermatfile = f"{material}_abq_deg{degree}mini_{law}_{protomodel}.inp"
        else:
            usermatfile = f"{material}_abq_deg{degree}mini_{law}_{protomodel}_{var_optim}_{n_try}.inp"
    filename = run_dir + sep + test + f"_{input_type}.inp"
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

def change_paras(test, material):
    type_test = test.split("_")[0]
    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material
    facs_displ = {"UT": 20, "NT6" : 2., "NT20": 2., "CH": 2., "SH": 1.}
    facs_thick = {"UT": 20, "NT6" : 2., "NT20": 2., "CH": 1., "SH": 8/5}
    facs_width = {"UT": 20, "NT6" : 10, "NT20": 10, "CH": 20, "SH": 20}

    res_filename = results_exp_dir + sep + test + "_1.csv"
    df_exp = pd.read_csv(res_filename, index_col=False)
    dtime = np.max(df_exp.iloc[:,0])
    if type_test != "UT":
        displ = np.max(df_exp.iloc[:,1]) / facs_displ[type_test]
        thickness = df_exp.loc[0, "Thickness[mm]"] / facs_thick[type_test]
        width = df_exp.loc[0, "OuterWidth[mm]"] / facs_width[type_test]
    else:
        displ = 0.1
        thickness = 1
        width = 1
    tramp1 = dtime / 100
    ramp1 = displ / 1000
    dt = dtime / 20000

    para_filename = run_dir + sep + test + "_paras.inp"
    with open(para_filename, "w") as f:
        f.write("*PARAMETER\n")
        f.write(f"DISPL = {displ}\n")
        f.write(f"DTIME = {dtime}\n")
        f.write(f"THICKNESS = {thickness}\n")
        f.write(f"WIDTH = {width}\n")
        f.write(f"TRAMP1 = {tramp1}\n")
        f.write(f"RAMP1 = {ramp1}\n")
        f.write(f"DT = {dt}\n")
        f.write("MINDT = DTIME * 1e-6\n")
        f.write("MAXDT = DTIME * 1e-2")
        f.close()



def run(test, func, input_type, law, var_optim=0, n_try=0):
    """
        Run the abaqus simulation of the test_polyN.inp and generate a csv file in results_sim folder
        Input :
            - test : string, name of the test (ex : UT_00)
    """
    print(f"Simulating {test}")
    job = f'{test}_{input_type}'
    if func=="polyN":
        subroutine = f"{input_type}_PolyN_3D_{law}.for"
    elif func== "polyN_mini":
        subroutine = f"{input_type}_PolyNmini_3D_{law}.for"

    input_file = f"{job}.inp"
    cp_input_file = f'temp_{job}_{var_optim}_{n_try}.inp'
    odb = "{}.odb".format(job)

    copy_sim_cmd = f'copy {input_file} {cp_input_file} '
    subprocess.call(copy_sim_cmd, shell=True, cwd=run_dir)

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    job_command = f"abq2023 job={cp_input_file} double interactive user={subroutine} ask_delete=OFF && exit"
    sim_command = f"{abq} {job_command}"
    
    subprocess.call(sim_command, shell=True, cwd=run_dir)

def post_process(test, material, input_type, is_optim, var_optim=0, n_try=0):
    job = f'{test}_{input_type}_{var_optim}_{n_try}'
    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material

    cp_input_file = f'temp_{job}.odb'
    odb = "{}.odb".format(job)

    copy_odb = f'copy {cp_input_file} {odb}'
    subprocess.call(copy_odb, shell=True, cwd=run_dir)

    delete_cmd = f'del temp_{job}*'
    subprocess.call(delete_cmd, shell=True, cwd=run_dir)

    print("Simulation {} ended".format(job))

    #REPORTING
    print("Creating report {}".format(test))
    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    if is_optim:
        load_pt = load_points[test.split("_")[0]]
        histregion=f'"Node PART-1-1.{load_pt}"'
        history = "U2,RF2"
        print(job, odb)
        command = f"abq2023 odbreport job={job} odb={odb} histregion={histregion} history={history} components"
        sim_command = f"{abq} {command}"
        result = subprocess.call(sim_command, shell=True, cwd=run_dir)
        print("Report {} ended".format(job))

        #POST PROCESSING
        print("Post processing report {}".format(test))
        filepath = run_dir + sep + r"{}.rep".format(job) 

        t = time.time()

        with open(filepath, "r") as file:
            file_no_space = file.read().replace("\n", " ")
            for i in range(10):
                file_no_space = file_no_space.replace("  ", " ")
            content = file_no_space.split(" ")
        
        len_data = 0
        labels = ["Time"]
        vu = 0
        n_frame = 0
        
        i = 0
        while content[i] != "End":
            while not((content[i] == "History") and (content[i+1] == "Output")):
                i = i + 1
            i = i + 2
            label = content[i][1:-1]
            labels.append(label)

            while content[i] != "type":
                i = i + 1

            i = i + 2
            type_label = content[i]

            if type_label == "SCALAR":
                len_data = 1
            
            while content[i] != "Data":
                i = i + 1
            
            i = i + 1

            if vu == 0:
                ts = []
                Y = []

                while content[i] != "End" and content[i] != "History":
                    ts.append(float(content[i]))
                    i = i + 1
                    Y.append(content[i:i+len_data])
                    i = i + len_data

                vu = 1
                n_frame= len(ts)
                data = np.zeros((n_frame, 1 + len_data))
                data[:,0] = np.array(ts)
                try :
                    data[:,1:1+len_data] = np.array(Y, dtype=float)
                except ValueError as e:
                    data[:,1:1+len_data] = np.zeros(len(Y))

            else :
                X = content[i:i + (1 + len_data) * n_frame]

                # Initialize a list to store converted values
                converted_values = []

                for val in X:
                    try:
                        converted_values.append(float(val))
                    except ValueError:
                        converted_values.append(0.0)  # Set to zero if conversion fails

                # Convert the list to a numpy array and reshape
                Y = np.array(converted_values).reshape(-1, 1 + len_data)[:, range(1, len_data + 1)]

                i = i + (1 + len_data) * n_frame
                data = np.concatenate((data, Y), axis=1)
                
        df = pd.DataFrame(data, columns = labels)
        filename = "{}.csv".format(job)
        if not(os.path.exists(results_sim_dir)):
            os.makedirs(results_sim_dir)
        filepath = results_sim_dir + sep + filename
        print(filepath)
        print("Post processing {} ended".format(job))
        df.to_csv(filepath)

        delete_cmd = f'del {test}*{var_optim}_{n_try}*'
        subprocess.call(delete_cmd, shell=True, cwd=run_dir)

    else:
        field = "S,LE,SDV_EPBAR"
        command = f"abq2023 odbreport job={job} odb={odb} field={field} components"
        full_command = f"{abq} {command}"
        result = subprocess.run(full_command, shell=True, capture_output=True, text=True, cwd=run_dir)

        print("Standard Output:\n", result.stdout)
        print("Standard Error:\n", result.stderr)
        print("Report {} ended".format(job))

        #POST PROCESSING
        print("Post processing report {}".format(test))
        filepath = run_dir + sep + r"{}.rep".format(job) 

        t = time.time()

        with open(filepath, "r") as file:
            file_no_space = file.read().replace("\n", " ")
            for i in range(10):
                file_no_space = file_no_space.replace("  ", " ")
            content = file_no_space.split(" ")
        
        i = 0
        pos_start_file = 0
        len_frame = 0
        pos_time = 0
        n_field = 0
        pos_start_field = []
        n_comp_field = []
        n_comp_cumul = 0
        n_comp_field_cumul = []

        labels = ["Time"]
        n_labels = 0

        #DELETE HEADER
        while content[i] != "Frame":
            i = i + 1
        
        pos_start_file = i - 1
        content = content[pos_start_file:]

        #SIZE OF A FRAME
        i = 2 #---- AND FIRST FRAME IGNORED

        while content[i] != "Frame" :
            i = i + 1

        len_frame = i - 1
        i = 0

        #TIME POSITION IN THE STEP
        while content[i] != "Time" :
            i = i + 1
        
        pos_time = i + 2

        #NUMBER OF FIELD
        while i < len_frame :
            if content[i] == "Field" :
                n_field = n_field + 1
                pos_start_field.append(i)
            
            if content[i] == "type":
                i = i + 2
                if content[i] == "SCALAR":
                    n_comp_field.append(1)
                    n_comp_cumul = n_comp_cumul + 1
                else :
                    n_comp_field.append(6)
                    n_comp_cumul = n_comp_cumul + 6
                n_comp_field_cumul.append(n_comp_cumul)
                
            if content[i] == "IP":
                i = i + 1
                n_comp = n_comp_field[n_field - 1]
                if n_comp == 6:
                    for k in range(n_comp):
                        labels.append(content[i])
                        i = i + 1
                else :
                    pos_label = pos_start_field[n_field - 1] + 2
                    labels.append(content[pos_label][1:-1])
        
            i = i + 1
        
        pos_start_field.append(len_frame)
        
        n_labels = len(labels)
        n_data = len(content) // len_frame

        data = np.zeros((n_data, n_labels))
        data[:, 0] = np.array(content[pos_time::len_frame], dtype=float)
        
        for j in range(n_field):
            n_comp = n_comp_field[j]
            n_comp_tot = n_comp_field_cumul[j]
            for k in range(n_comp):
                data[:, n_comp_field_cumul[j] - k] = np.array(content[(pos_start_field[j + 1] - k - 1)::len_frame], dtype=float)
        
        df = pd.DataFrame(data, columns=labels)

        for i in range(1, 4):
            for j in range(i, 4):
                e = "E{}{}".format(i,j)
                le = "LE{}{}".format(i,j)
                df[e] = np.exp(df[le]) - 1
        
    
        filename = "{}.csv".format(job)
        if not(os.path.exists(results_sim_dir)):
            os.makedirs(results_sim_dir)
        filepath = results_sim_dir + sep + filename
        print("Post processing {} ended".format(job))
        df.to_csv(filepath)
        

def launch_run(tests, func, material, degree, law, protomodel, input_type, var_optim=0, n_try=0):
    for test in tests:
        change_usermat(test, func, material, degree, law, protomodel, input_type, var_optim, n_try)
        change_paras(test, material)

    time.sleep(5)
    for test in tests:
        run(test, func, input_type, law, var_optim, n_try)

def create_csv(tests, material, input_type, var_optim=0, n_try=0):
    for test in tests:
        #0 to ut report field, 1 to optim
        post_process(test, material, input_type, 0, var_optim, n_try)


p = read_param()

func = p["func"]
material = p["material"]
degree = int(p["degree"])
law = p["law"]
protomodel = p["protomodel"]
input_type = p["input_type"]

launch_run(ut_tests, func, material, degree, law, protomodel, input_type, var_optim=0, n_try=0)
create_csv(ut_tests, material, input_type, var_optim=0, n_try=0)