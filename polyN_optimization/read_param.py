import os

file_dir = os.path.dirname(os.path.abspath(__file__))
filename = "param.txt"
sep = os.sep

def read_param():
    """
        Read the parameters in the "param.txt" file which has to be located in the same directory as this script
        Output :
            - dic : dictionnary {"variable" : "value"}
    """
    dic = {}
    filepath = file_dir + sep + filename
    with open(filepath, "r") as f:
        params = f.readlines()[1:]
        for l in params:
            if not(l == "\n"):
                v, _, value = l.strip("\n").split(" ")
                dic[v] = value
    return(dic)