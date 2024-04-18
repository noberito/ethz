import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

###In the same directory, put a folder of an experiment which contains a directory "DATA" containing
###CSV files named "test_angle_exp.csv". Create also a directory "PLOTS" to put plots

sigma = "\u03c3"
epsilon = "\u03B5"
##folder of the data
experience="2018_11_09_G_DP-DP780_1d50mm"

##filename
filenames=["UT_00_exp", "UT_15_exp", "UT_30_exp", "UT_45_exp", "UT_60_exp", "UT_75_exp", "UT_90_exp"]

##columns name
variables=["Displacement[mm]", "Axial Force[kN]"]

def uts_plot(filenames, experience, variables):
    for filename in filenames :
        filepath = ".\\" + experience + "\\DATA\\" + filename + ".csv"
        db=pd.read_csv(filepath, names=variables, usecols=range(len(variables)))
        plt.plot(db[variables[0]], db[variables[1]], label = filename)
        plt.xlabel(variables[0], fontsize="large")
        plt.ylabel(variables[1], fontsize="large")
    plt.legend()
    plt.show()

uts_plot(filenames, experience, variables)