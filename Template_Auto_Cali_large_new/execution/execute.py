import subprocess
import os
import multiprocessing
import pandas as pd
import time
import numpy as np

# Get the directory where the Python exec is located



current_dir = "./"  # Assuming current directory
dir = "/"
if os.name == "nt":
    current_dir = ".\\"
    dir = "\\"

exec_dir = os.path.dirname(os.path.abspath(__file__))
template_dir = os.path.dirname(exec_dir)
res_sim_dir = template_dir + dir + "results_sim"

jobs = ["UT_00_polyN", "UT_15_polyN", "UT_30_polyN", "UT_45_polyN", "UT_60_polyN", "UT_75_polyN","UT_90_polyN", "UT_EBT_polyN"]
""""CH_00_polyN", "CH_45_polyN",
"NT6_00_polyN", "NT6_45_polyN", "NT6_90_polyN",
"NT20_00_polyN", "NT20_45_polyN",
"SH_000_polyN", "SH_090_polyN", "SH_p45_polyN"]"""

def simulate(job):
    print("Simulating {}".format(job))
    subroutine = "abqUMAT_PolyN_3D.for"
    input_file = "{}.inp".format(job)

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'

    command = f"abq2023 job={job} input={input_file} double user={subroutine} && exit"
    full_command = f"{abq} {command}"
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True, cwd=exec_dir)

    print("Standard Output:\n", result.stdout)
    print("Standard Error:\n", result.stderr)

    print("Simulation {} ended".format(job))

def report(job):
    print("Creating report {}".format(job))
    odb = "{}.odb".format(job)
    field = "S,LE,SDV_EPBAR"

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    command = f"abq2023 odbreport job={job} odb={odb} field={field} components"
    full_command = f"{abq} {command}"
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True, cwd=exec_dir)

    print("Standard Output:\n", result.stdout)
    print("Standard Error:\n", result.stderr)
    print("Report {} ended".format(job))

def post_process_report(job):
    print("Post processing report {}".format(job))

    filename = exec_dir + dir + r"{}.rep".format(job) 
    t = time.time()
    with open(filename, "r") as file:
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
    filepath = res_sim_dir + dir + filename
    print("Post processing {} ended".format(job))
    df.to_csv(filepath)

if __name__ == "__main__":
    pool = multiprocessing.Pool()
    pool.map(simulate, jobs)
    pool.map(report, jobs)
    pool.map(post_process_report, jobs)

    pool.close()
    pool.join()
