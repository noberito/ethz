import pandas as pd

df = pd.read_csv('output.txt', delimiter="                    ", names=["TIME", "HF", "YF", "EPBAR", "SIGMA(1)"])
filename = "output.xlsx"

df.to_excel(filename, index=False)