import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

sigma = "\u03c3"
epsilon = "\u03B5"

###In the same directory, put a folder of an experiment which contains a directory "DATA" containing
###CSV files named "test_angle_exp.csv". Create also a directory "PLOTS" to put plots

##folder of the data
experience="2018_11_09_G_DP-DP780_1d50mm"

##filename
filename="UT_EBT_exp"

##columns name
variables=["Eng. axial strain[-]", "Eng. axial stress[MPA]"]

def ut_plot(filename, experience, variables):
    output_directory = ".//" + experience + "//PLOTS//"
    filepath = ".//" + experience + "//DATA//" + filename + ".csv"
    db=pd.read_csv(filepath, names=variables, usecols=range(len(variables)))
    plt.plot(db[variables[0]], db[variables[1]])
    plt.xlabel(variables[0])
    plt.ylabel(variables[1])
    plt.title(filename)
    plot_filename = os.path.join(output_directory, filename)
    plt.savefig(plot_filename)

ut_plot(filename, experience, variables)