import subprocess
import os
import sys
import multiprocessing
import pandas as pd
import time
import numpy as np
import glob


# Get the directory where the Python exec is located

sep = os.sep
polyN_cali_dir = os.path.dirname(os.path.abspath(__file__))
run_dir = polyN_cali_dir + sep + "running"

sys.path.append(polyN_cali_dir)

from read_param import read_param
from get_calibration_data import analyze_exp_data
from tests_parameters import load_points, ext_points, ut_tests_ext
from check.compare_sim_exp import compare_large_strain, compare_ut_s_2
from optimize_polyN_mini import write_coeff_abq_mini, get_param_polyN_mini

facs_displ = {"UT": 1, "NT6" : 2., "NT20": 2., "CH": 2., "SH": 1}
facs_thick = {"UT": 1, "NT6" : 1., "NT20": 1., "CH": 0.5, "SH": 0.8}
facs_width = {"UT": 1, "NT6" : 20, "NT20": 20, "CH": 20, "SH": 10}
facs_force = {"UT": 1, "NT6" : 4, "NT20": 4, "CH": 4, "SH": 2}

dt = 1

def change_paras(test, material):
    type_test = test.split("_")[0]
    results_exp_dir = polyN_cali_dir + sep + "results_exp" + sep + material


    res_filename = results_exp_dir + sep + test + "_1.csv"
    df_exp = pd.read_csv(res_filename, index_col=False)
    dtime = np.max(df_exp.iloc[:,0])

    if type_test != "UT":
        displ = np.max(df_exp.iloc[:,1]) / facs_displ[type_test]
        thickness = df_exp.loc[0, "Thickness[mm]"] / (2 * facs_thick[type_test])
        width = df_exp.loc[0, "OuterWidth[mm]"] / facs_width[type_test]
    else:
        displ = 0.01
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

def change_usermat(test, func, material, degree, law, protomodel, input_type, p=0, m=0):
    """
        Change the input user material file in the all the abaqus files of test in tests
        Input :
            - usermatfile : string, filename of the user material file
    """
    if func=="polyN":
        if str(p)=="0" and m==0:
            usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}.inp"
        else:
            usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}_{p}_{m}.inp"
    if func=="polyN_mini":
        if str(p)=="0" and m==0:
            usermatfile = f"{material}_abq_deg{degree}mini_{law}_{protomodel}.inp"
        else:
            usermatfile = f"{material}_abq_deg{degree}mini_{law}_{protomodel}_{p}_{m}.inp"
            
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

def simulate(test, func, input_type, law,  p=0, m=0):
    """
        Run the abaqus simulation of the test_polyN.inp and generate a csv file in results_sim folder
        Input :
            - test : string, name of the test (ex : UT_00)
    """
    print(f"Simulating {test}")

    job = f'{test}_{input_type}'
    cpus = 24

    if func == "polyN":
        subroutine = f"{input_type}_PolyN_3D_{law}.for"
    elif func == "polyN_mini":
        subroutine = f"{input_type}_PolyNmini_3D_{law}.for"

    type_test = test.split("_")[0]

    input_file = f"{job}.inp"
    cp_input_file = f'temp_{job}_{p}_{m}.inp'

    print(input_file)

    copy_sim_cmd = f'cp {input_file} {cp_input_file} '
    subprocess.call(copy_sim_cmd, shell=True, cwd=run_dir)

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    if type_test == "UT":
        job_command = f'''sbatch -C ib -n 1 -t 0-2 --mem-per-cpu=300 --tmp=300 --wrap "bash -c 'abaqus job={cp_input_file} double interactive user={subroutine} cpus=1 scratch=\$TMPDIR'"'''
    else:
        job_command = f'''sbatch -C ib -n {cpus} --time=00:10:00 --mem-per-cpu=300 --tmp=300 --wrap "bash -c 'abaqus job={cp_input_file} double interactive user={subroutine} cpus={cpus} scratch=\$TMPDIR'"'''

    res = subprocess.run(job_command, shell=True, cwd=run_dir, capture_output=True)
    n_batch = res.stdout.decode('utf-8').strip("\n").split(" ")[-1]
    return(n_batch)

def post_process(test, material, input_type, p=0, m=0):
    """
        Run the abaqus simulation of the test_polyN.inp and generate a csv file in results_sim folder
        Input :
            - test : string, name of the test (ex : UT_00)
    """

    results_sim_dir = polyN_cali_dir + sep + "results_sim" + sep + material
    job = f'temp_{test}_{input_type}_{p}_{m}'
    odb = f'temp_{test}_{input_type}_{p}_{m}.odb'

    print("Simulation {} ended".format(job))

    #REPORTING
    print("Creating report {}".format(test))
    type_test = test.split("_")[0]

    if type_test == "UT":
        field = "S,LE,SDV_EPBAR"
        command = f"abaqus odbreport job={job} odb={odb} field={field} components"
        result = subprocess.run(command, shell=True, cwd=run_dir)
        time.sleep(10)
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
                frame_data = content[(pos_start_field[j + 1] - k - 1)::len_frame]
                try:
                    data[:, n_comp_field_cumul[j] - k] = np.array(frame_data, dtype=float)
                except ValueError as e:
                    data[:, n_comp_field_cumul[j] - k] = np.zeros(len(frame_data))

        df = pd.DataFrame(data, columns=labels)

        for i in range(1, 4):
            for j in range(i, 4):
                e = "E{}{}".format(i,j)
                le = "LE{}{}".format(i,j)
                df[e] = np.exp(df[le]) - 1
    
        filename = "{}_{}_{}_{}.csv".format(test, input_type, p, m)
        if not(os.path.exists(results_sim_dir)):
            os.makedirs(results_sim_dir)
        filepath = results_sim_dir + sep + filename
        print("Post processing {} ended".format(job))
        df.to_csv(filepath)
   

    else :
        load_pt = load_points[test.split("_")[0]]
        histregion=f'"Node PART-1-1.{load_pt}"'
        history = "U2,RF2"
        command = f"abaqus odbreport job={job} odb={odb} histregion={histregion} history={history} components"
        result = subprocess.run(command, shell=True, cwd=run_dir)
        time.sleep(10)
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


        if test.split("_")[0] in ext_points.keys():
            ext_pt = ext_points[test.split("_")[0]]
            histregion=f'"Node PART-1-1.{ext_pt}"'
            history = "U2"
            command = f"abaqus odbreport job={job} odb={odb} histregion={histregion} history={history} components"
            result = subprocess.run(command, shell=True, cwd=run_dir)
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
            labels.append("Time")
            vu = 0
            n_frame = 0
            
            i = 0
            while content[i] != "End":
                while not((content[i] == "History") and (content[i+1] == "Output")):
                    i = i + 1
                i = i + 2
                label_ext = content[i][1:-1] + "_ext"
                labels.append(label_ext)

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
                    data_ext = np.zeros((n_frame, 1 + len_data))
                    data_ext[:,0] = np.array(ts)
                    try :
                        data_ext[:,1:1+len_data] = np.array(Y, dtype=float)
                    except ValueError as e:
                        data_ext[:,1:1+len_data] = np.zeros(len(Y))

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
                    data_ext = np.concatenate((data_ext, Y), axis=1)

            data = np.concatenate((data, data_ext), axis = 1)

        df = pd.DataFrame(data, columns = labels)
        if test.split("_")[0] in ext_points.keys():
            df["Strain_eng"] = (df["U2_ext"] / 0.5)
            df["Strain_ext"] = np.log(1 + df["Strain_eng"])
        df["U2"] = df["U2"] * facs_displ[type_test]
        df["RF2"] = df["RF2"] / 1000 * facs_force[type_test]
        filename = "{}_{}_{}_{}.csv".format(test, input_type, p, m)
        if not(os.path.exists(results_sim_dir)):
            os.makedirs(results_sim_dir)
        filepath = results_sim_dir + sep + filename
        print("Post processing {} ended".format(job))
        df.to_csv(filepath)

def del_temp_file(test, input_type, p=0, m=0):
    pattern = f"temp_{test}_{input_type}_{p}_{m}*"
    os.chdir(run_dir)
    files_to_delete = glob.glob(pattern)

    for file in files_to_delete:
        if os.path.isdir(file):
            subprocess.run(['rm', '-r', file])
        else:
            subprocess.run(['rm', file])

def del_usermatfile(material, func, degree, law, protomodel, p=0, m=0):
    if func=="polyN":
        if str(p)=="0" and m==0:
            usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}.inp"
        else:
            usermatfile = f"{material}_abq_deg{degree}_{law}_{protomodel}_{p}_{m}.inp"
    if func=="polyN_mini":
        if str(p)=="0" and m==0:
            usermatfile = f"{material}_abq_deg{degree}mini_{law}_{protomodel}.inp"
        else:
            usermatfile = f"{material}_abq_deg{degree}mini_{law}_{protomodel}_{p}_{m}.inp"
    
    filepath = run_dir + sep + usermatfile
    subprocess.run(["rm", filepath])


def sim_finished(test, input_type, p=0, m=0):
    test_file = f"temp_{test}_{input_type}_{p}_{m}.lck"
    test_filepath = run_dir + sep + test_file
    if os.path.exists(test_filepath):
        return(False)
    return(True)

def has_failed(batch):
    out_file = "slurm-" + batch + ".out"
    out_filepath = run_dir + sep + out_file
    with open(out_filepath, "r") as f :
        lines = f.readlines()
        if lines[-2] == "Abaqus Error: Problem during linking - Abaqus/Standard User Subroutines\n":
            return(True)
        return(False)
    
def batch_running(batch):
    out_file = "slurm-" + batch + ".out"
    out_filepath = run_dir + sep + out_file
    return(os.path.exists(out_filepath))

def del_lck_file(test, input_type, p=0, m=0):
    test_file = f"temp_{test}_{input_type}_{p}_{m}.lck"
    test_filepath = run_dir + sep + test_file
    subprocess.run(['rm', test_filepath])

def launch_run(tests, func, material, degree, law, protomodel, input_type, p=0, m=0):
    #Run les tests disponibles experimentalement

    time_limit = 360
    n = len(tests)
    
    for test in tests:
        #Creating and changing paras file
        change_paras(test, material)

        #Changing the user material file
        change_usermat(test, func, material, degree, law, protomodel, input_type, p, m)

    batchs = []

    #Launching simulations job
    for test in tests:
        batch = simulate(test, func, input_type, law, p, m)
        batchs.append(batch)

    #Waiting for simulations to start
    for i in range(n):
        while not(batch_running(batchs[i])):
            time.sleep(1)

    t0 = time.time()
    tests_failure = []
    tests_success = []

    #Average time to be able to see if one simulation crashed
    time.sleep(30)

    for i in range(n): 
        if has_failed(batchs[i]):
            tests_failure.append(tests[i])
        else:
            tests_success.append(tests[i])

    #Relaunch the failed simulations
    if len(tests_failure) != 0:
        launch_run(tests_failure, func, material, degree, law, protomodel, input_type, p, m)
    
    #Waiting for simulations to end
    for test in tests_success:
        while not(sim_finished(test, input_type, p, m)):
            dt = time.time() - t0
            print(dt)
            if dt < time_limit :
                time.sleep(1)
            else :
                del_lck_file(test, input_type, p, m)

    #del_usermatfile(material, func, degree, law, protomodel, p, m)

def create_csv(tests, material, input_type, p=0, m=0):

    #Tests post processing
    for test in tests:
        post_process(test, material, input_type, p, m)
        del_temp_file(test, input_type, p, m)

if __name__ == "__main__":
    p = read_param()
    material = p["material"]
    func = p["func"]    
    gseed = int(p["gseed"])
    input_type = p["input_type"]
    enu = float(p["enu"])
    density = float(p["density"])
    nb_virtual_pt = int(p["nb_virtual_pt"])
    degree = int(p["degree"])
    sh = int(p["sh"])
    nt6 = int(p["nt6"])
    weight_ut = float(p["weight_ut"])
    weight_vir = float(p["weight_vir"])
    weight_e2 = float(p["weight_e2"])
    protomodel = p["protomodel"]
    law = p["law"]
    mat_exp = analyze_exp_data(material)

    tests = []
    for type_test in mat_exp.keys():
        for ori in mat_exp[type_test]:
            if type_test != "UT":
                test = type_test + "_" + ori
                tests.append(test)

    coeff_polyN_mini = np.load(polyN_cali_dir + sep + material + "_poly" + str(degree) + "_mini_coeff.npy")
    coeff_law = np.load(polyN_cali_dir + sep + material + f"_{law}_mini_coeff.npy")

    a = coeff_law[0]
    b = coeff_law[1]
    c = coeff_law[2]
    ymod = coeff_law[3]
    nmon = len(coeff_polyN_mini) - 2
    powers = get_param_polyN_mini(degree)

    write_coeff_abq_mini(coeff_polyN_mini, a, b, c, ymod, enu, protomodel, degree, material, law, density, powers, 10, 10)
    launch_run(["NT6_45"], func, material, degree, law, protomodel, input_type, 10, 10)

    #launch_run(tests, func, material, degree, law, protomodel, input_type, 10, 10)
    #create_csv(tests, material, input_type, 10, 10)
    #compare_large_strain(material, func, degree, input_type, 10, 10)