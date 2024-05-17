import pandas as pd
import numpy as np

df = pd.read_csv('output.txt', delimiter="               ", names=["TIME", "EPBAR", "S22", "YF", "HF", "AA", "BB", "CC"], dtype="float")
df["VOCE"] = df["AA"] - df["BB"] * (1 - np.exp(- df["CC"] * df["EPBAR"]))
filename = "output.xlsx"
df.to_excel(filename, index=False, float_format="%.8f")