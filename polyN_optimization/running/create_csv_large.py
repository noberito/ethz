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

load_points = {"UT" :"7", "CH" :"25791", "SH":"8676", "NT6" : "7382", "NT20" : "1683"}

def post_process(test):
    """
        Run the abaqus simulation of the test_polyN.inp and generate a csv file in results_sim folder
        Input :
            - test : string, name of the test (ex : UT_00)
    """
    job = f'{test}_{input_type}'

    input_file = f"{job}.inp"
    odb = "{}.odb".format(job)

    copy_odb = f'copy temp_{odb} {odb}'
    subprocess.call(copy_odb, shell=True, cwd=exec_dir)

    delete_cmd = f'del temp_{job}*'
    subprocess.call(delete_cmd, shell=True, cwd=exec_dir)


    #REPORTING
    print("Creating report {}".format(test))

    if test.split("_")[0] == "UT":
        field = "S,LE,SDV_EPBAR"
        abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
        command = f"abq2023 odbreport job={job} odb={odb} field={field} components"
        full_command = f"{abq} {command}"
        result = subprocess.run(full_command, shell=True, cwd=exec_dir)

        print("Report {} ended".format(job))

        #POST PROCESSING
        print("Post processing report {}".format(test))
        filepath = exec_dir + sep + r"{}.rep".format(job) 

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
    
        filename = "{}.csv".format(job)
        if not(os.path.exists(results_sim_dir)):
            os.makedirs(results_sim_dir)
        filepath = results_sim_dir + sep + filename
        print("Post processing {} ended".format(job))
        df.to_csv(filepath)
   

    else :

        load_pt = load_points[test.split("_")[0]]
        histregion=f'"Node PART-1-1.{load_pt}"'
        history = "U2,RF2"
        abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
        command = f"abq2023 odbreport job={job} odb={odb} histregion={histregion} history={history} components"
        full_command = f"{abq} {command}"
        result = subprocess.run(full_command, shell=True, capture_output=True, text=True, cwd=exec_dir)
        print("Report {} ended".format(job))

        #POST PROCESSING
        print("Post processing report {}".format(test))
        filepath = exec_dir + sep + r"{}.rep".format(job) 

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
        print("Post processing {} ended".format(job))
        df.to_csv(filepath)         



if __name__ == "__main__":
    mat_tests = analyze_exp_data(material)
    tests = []

    for type_test in mat_tests.keys():
        for ori in mat_tests[type_test]:
            if type_test != "UT":
                test = type_test + "_" + ori
                tests.append(test)
    if(1):
        tests = ["CH_00", "CH_45",
         "NT6_00", "NT6_45", "NT6_90",
         "NT20_00", "NT20_45",
         "SH_000", "SH_090", "SH_p45"]
        
    pool = multiprocessing.Pool()
    pool.map(post_process, tests)
    pool.close()
    pool.join()
