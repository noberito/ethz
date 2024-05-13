import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.preprocessing

A = np.eye(6,6)
B = np.tile(np.expand_dims(A, 0), (22, 1, 1))

print(B[1])