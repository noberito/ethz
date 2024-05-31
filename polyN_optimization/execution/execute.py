import subprocess
import os
import multiprocessing
import pandas as pd
import time
import numpy as np

# Get the directory where the Python exec is located

dir = os.sep

material = "DP600"
law = "voce"
degree = 6

exec_dir = os.path.dirname(os.path.abspath(__file__))
polyN_cali_dir = os.path.dirname(exec_dir)
results_sim_dir = polyN_cali_dir + dir + "results_sim" + dir + material

sim_params = ["UT_00", "UT_15", "UT_30", "UT_45", "UT_60", "UT_75","UT_90", "UT_EBT"]
""""CH_00", "CH_45",
"NT6_00", "NT6_45", "NT6_90",
"NT20_00", "NT20_45",
"SH_000", "SH_090", "SH_p45"]"""

tests = ["UT_00", "UT_15", "UT_30", "UT_45", "UT_60", "UT_75","UT_90", "UT_EBT"]
""""CH_00", "CH_45",
"NT6_00", "NT6_45", "NT6_90",
"NT20_00", "NT20_45",
"SH_000", "SH_090", "SH_p45"]"""

dt = 1

def simulate(sim_param):
    test = sim_param
    print(f"Simulating {test}")
    job = f'{test}_polyN'
    subroutine = f"abqUMAT_PolyN_3D_v2_{law}.for"
    input_file = f"{job}.inp"
    cp_input_file = f'temp_{input_file}'
    odb = "{}.odb".format(job)

    copy_sim_cmd = f'copy {input_file} {cp_input_file} '
    subprocess.call(copy_sim_cmd, shell=True, cwd=exec_dir)

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    job_command = f"abq2023 job=temp_{job} double interactive user={subroutine} && exit"
    sim_command = f"{abq} {job_command}"
    
    subprocess.call(sim_command, shell=True, cwd=exec_dir)

    copy_odb = f'copy temp_{odb} {odb}'
    subprocess.call(copy_odb, shell=True, cwd=exec_dir)

    delete_cmd = f'del temp_{job}*'
    #subprocess.call(delete_cmd, shell=True, cwd=exec_dir)

    print("Simulation {} ended".format(job))

    #REPORTING
    print("Creating report {}".format(test))
    job = f'{test}_polyN'
    sim_dir = exec_dir + dir + test

    odb = "{}.odb".format(job)
    field = "S,LE,SDV_EPBAR"
    
    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    command = f"abq2023 odbreport job={job} odb={odb} field={field} components"
    full_command = f"{abq} {command}"
    result = subprocess.run(full_command, shell=True, capture_output=True, text=True, cwd=exec_dir)

    print("Standard Output:\n", result.stdout)
    print("Standard Error:\n", result.stderr)
    print("Report {} ended".format(job))

    #POST PROCESSING
    print("Post processing report {}".format(test))
    job = f'{test}_polyN'
    filepath = exec_dir + dir + r"{}.rep".format(job) 

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
    filepath = results_sim_dir + dir + filename
    print("Post processing {} ended".format(job))
    df.to_csv(filepath)

if __name__ == "__main__":
    """for test in tests :
        simulate(test)"""
    pool = multiprocessing.Pool()
    pool.map(simulate, sim_params)
    pool.close()
    pool.join()
