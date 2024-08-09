import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

polyN_dir = os.path.dirname(os.path.abspath(__file__))
filename = "param.txt"
sep = os.sep

def read_param():
    """
        Read the parameters in the "param.txt" file which has to be located in the same directory as this script
        Output :
            - dic : dictionnary {"variable" : "value"}
    """
    dic = {}
    filepath = polyN_dir + sep + filename
    with open(filepath, "r") as f:
        params = f.readlines()[1:]
        for l in params:
            if not(l == "\n"):
                try:
                    v, _, value = l.strip("\n").split(" ")
                    dic[v] = value
                except ValueError :
                    pass
                
    return(dic)

def readdata_exp(material):
    """
        Returns a dataframe containing the experimental points from UT from the given material 
        and virtual points using the given protomodel.
        df.columns = [#TODO]
    
        Input :
            - material : string
            - protomodel : string

        Output :
            - df : dataFrame
    """

    folderpath = f"{polyN_dir}{sep}calibration_data{sep}{material}"
    filename_e = "data_exp_" + material + ".csv"
    filepath_e = folderpath + sep + filename_e
    df_e = pd.read_csv(filepath_e)

    df_e["LoadAngle"] = 2 * np.pi / 360 * df_e["LoadAngle"]
    sigma0 = df_e["YieldStress"].iloc[0]
    
    df_e["s11"] = df_e["YieldStress"] / sigma0 * (np.square(np.cos(df_e["LoadAngle"])) + df_e["q"] * np.square(np.sin(df_e["LoadAngle"])))
    df_e["s22"] = df_e["YieldStress"] / sigma0 * (np.square(np.sin(df_e["LoadAngle"])) + df_e["q"] * np.square(np.cos(df_e["LoadAngle"])))
    df_e["s33"] = df_e["YieldStress"] * 0
    df_e["s12"] = df_e["YieldStress"] / sigma0 * (1 - df_e["q"]) * np.sin(df_e["LoadAngle"]) * np.cos(df_e["LoadAngle"])
    df_e["s13"] = df_e["YieldStress"] * 0
    df_e["s23"] = df_e["YieldStress"] * 0
    
    return(df_e)

def readData_2d(material, protomodel):
    """
        Returns a dataframe containing the experimental points from UT from the given material 
        and virtual points using the given protomodel.
        df.columns = [#TODO]
    
        Input :
            - material : string
            - protomodel : string

        Output :
            - df : dataFrame
    """

    df_e = readdata_exp(material)

    folderpath = f"{polyN_dir}{sep}calibration_data{sep}{material}"

    try:
        filename_v = "data_virtual_" + material + "_" + protomodel + ".csv"
        filepath_v = folderpath + sep + filename_v
        df_v = pd.read_csv(filepath_v)
    except:
        print(f"protomodel {protomodel} not available for the moment")
        df_v = pd.DataFrame()
    
    if 0:
        if protomodel == "bezier":
            data = df_v[["s11", "s22", "s12"]].values

            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            X = data[:,0]
            Y = data[:,1]
            Z = data[:,2]
            ax.scatter(X, Y, Z, marker="o", s=1)

            plt.show()

    df = pd.concat([df_e, df_v])
    df["Norm"] = np.linalg.norm(df[["s11", "s22", "s33", "s12", "s13", "s23"]].values, axis=1)

    return(df)

def get_coeff_mini(material, degree):
    filename = f"{material}_poly{degree}_mini.npy"
    filedir = polyN_dir + sep + "coeff"
    filepath = filedir + sep + filename
    coeff_mini = np.load(filepath)
    return(coeff_mini)

def get_coeff_mini_opti(material, degree, opti):
    filename = f"{material}_poly{degree}_mini_{opti}.npy"
    filedir = polyN_dir + sep + "coeff"
    filepath = filedir + sep + filename
    coeff_mini = np.load(filepath)
    return(coeff_mini)

def get_coeff_law(material, law):
    filename = f"{material}_{law}.npy"
    filedir = polyN_dir + sep + "coeff"
    filepath = filedir + sep + filename
    coeff = np.load(filepath)
    coeff_law = coeff[:-1]
    ymod = coeff[-1]
    return(coeff_law, ymod)

def get_coeff_mises(degree):
    if degree==2:
        C = np.array([1, -1, 1, 3, 3, 3])
    elif degree==4:
        C = np.array([1, -2, 3, -2, 1, 6, -6, 6, 9, 3, 3])
    elif degree==6:
        C = np.array([1, -3, 6, -7, 6, -3, 1, 12, -24, 30, -24, 12, 14, -2, 14, 30, 3, 3])
    elif degree==8:
        C = np.array([1, -4, 14, -30, 40, -30, 14, -4, 1, 11, -21, 17, -1, 17, -21, 11, 40, -4, 18, -4, 40, 42, -6, 42, 120, 3, 3])
    else:
        C = []
    return(C)