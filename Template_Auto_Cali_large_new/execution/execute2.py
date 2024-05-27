import subprocess
import os
import multiprocessing
import pandas as pd
import time
import numpy as np

# Get the directory where the Python exec is located
current_dir = "./"  # Assuming current directory
dir_separator = "\\" if os.name == "nt" else "/"

exec_dir = os.path.dirname(os.path.abspath(__file__))
template_dir = os.path.dirname(exec_dir)
res_dir = template_dir + dir_separator + "results"
os.makedirs(res_dir, exist_ok=True)

jobs = [
    "UT_00_polyN", "UT_15_polyN", "UT_30_polyN", "UT_45_polyN", 
    "UT_60_polyN", "UT_75_polyN", "UT_90_polyN", "UT_EBT_polyN"
]
# "CH_00_polyN", "CH_45_polyN",
# "NT6_00_polyN", "NT6_45_polyN", "NT6_90_polyN",
# "NT20_00_polyN", "NT20_45_polyN",
# "SH_000_polyN", "SH_090_polyN", "SH_p45_polyN"

def run_command(full_command, timeout=300):
    try:
        result = subprocess.Popen(full_command, shell=True, capture_output=True, text=True, cwd=exec_dir, timeout=timeout)
        print("Standard Output:\n", result.stdout)
        print("Standard Error:\n", result.stderr)
        return result.returncode
    except subprocess.TimeoutExpired:
        print(f"Command timed out: {full_command}")
        return 1  # Return non-zero to indicate timeout/failure

def simulate(job):
    print("Simulating {}".format(job))
    subroutine = "abqUMAT_PolyN_3D.for"
    input_file = "{}.inp".format(job)

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    command = f"abq2023 job={job} input={input_file} double user={subroutine}"
    full_command = f"{abq} {command}"

    process = subprocess.Popen(full_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=exec_dir)
    stdout, stderr = process.communicate()  # Wait for command to complete

    if process.returncode == 0:
        print("Standard Output:\n", stdout.decode())
        print("Standard Error:\n", stderr.decode())
        print("Simulation {} ended successfully".format(job))
    else:
        print("Simulation {} failed".format(job))

    # Add delay
    time.sleep(10)  # Adjust the delay time as needed


def report(job):
    print("Creating report {}".format(job))
    odb = "{}.odb".format(job)
    field = "S,LE,SDV_EPBAR"

    abq = r'"C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\env\vars.bat" -arch intel64 vs2019 &'
    command = f"abq2023 odbreport job={job} odb={odb} field={field} components"
    full_command = f"{abq} {command}"

    if run_command(full_command) == 0:
        print("Report {} ended successfully".format(job))
    else:
        print("Report {} failed".format(job))

    # Add delay
    time.sleep(10)  # Adjust the delay time as needed

def post_process_report(job):
    print("Post processing report {}".format(job))

    filename = exec_dir + dir_separator + r"{}.rep".format(job)
    try:
        with open(filename, "r") as file:
            file_no_space = file.read().replace("\n", " ")
            for _ in range(10):
                file_no_space = file_no_space.replace("  ", " ")
            content = file_no_space.split(" ")
    except FileNotFoundError:
        print(f"Report file not found: {filename}")
        return

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

    # DELETE HEADER
    while content[i] != "Frame":
        i += 1
    
    pos_start_file = i - 1
    content = content[pos_start_file:]

    # SIZE OF A FRAME
    i = 2  # ---- AND FIRST FRAME IGNORED

    while content[i] != "Frame":
        i += 1

    len_frame = i - 1

    i = 0
    # TIME POSITION IN THE STEP
    while content[i] != "Time":
        i += 1
    
    pos_time = i + 2

    # NUMBER OF FIELD
    while i < len_frame:
        if content[i] == "Field":
            n_field += 1
            pos_start_field.append(i)
        
        if content[i] == "type":
            i += 2
            if content[i] == "SCALAR":
                n_comp_field.append(1)
                n_comp_cumul += 1
            else:
                n_comp_field.append(6)
                n_comp_cumul += 6
            n_comp_field_cumul.append(n_comp_cumul)
        
        if content[i] == "IP":
            i += 1
            n_comp = n_comp_field[n_field - 1]
            if n_comp == 6:
                for k in range(n_comp):
                    labels.append(content[i])
                    i += 1
            else:
                pos_label = pos_start_field[n_field - 1] + 2
                labels.append(content[pos_label][1:-1])
    
        i += 1
    
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
    filepath = res_dir + dir_separator + filename
    print("Post processing {} ended".format(job))
    df.to_csv(filepath)
    
    # Add delay
    time.sleep(10)  # Adjust the delay time as needed

if __name__ == "__main__":
    # Create a pool of workers with a smaller pool size to limit concurrent processes
    pool_size = 4 # Adjust pool size based on your system's capability
    pool = multiprocessing.Pool(pool_size)
    
    # Run simulations
    pool.map(simulate, jobs)
    
    # Run reports
    pool.map(report, jobs)
    
    # Post-process reports
    pool.map(post_process_report, jobs)

    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()
