import os

file_dir = os.path.dirname(os.path.abspath(__file__))
sep = os.sep

files = os.listdir(file_dir)
for file in files:
    if file[-4::] == ".inp":
        with open(file_dir + sep +file, "r") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                if lines[i] == "*Static\n":
                    lines[i+1] = "<MINDT>, <DTIME>, ,<MAXDT>\n"
        
        with open(file_dir + sep + file, "w") as f:
            f.writelines(lines)

