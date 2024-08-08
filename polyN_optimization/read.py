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