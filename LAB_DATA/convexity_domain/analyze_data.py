import pandas as pd
import sklearn
import sklearn.preprocessing
import matplotlib.pyplot as plt
import numpy as np
import time
import multiprocessing

degree = 4

X_coeff = np.load("X_convex_{}.npy".format(degree))
Y_coeff = np.load("Y_convex_{}.npy".format(degree))

print(X_coeff[Y_coeff > 0.5])