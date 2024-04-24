import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

a = np.array([0, 1, 2, 3, 4])
p = np.array([1, 1, 1, 1, 1])
ndata = 10000
weigth = np.ones(ndata)
X = np.ones((ndata,3))
M = np.zeros((3, 3))

print(a ** (p * 2))



for i in range(ndata):
    Xi = np.expand_dims(X[i], axis=0)
    M = M + weigth[i] * np.dot(Xi.T, Xi)
M = M / 2

V = np.sum(np.expand_dims(weigth, axis=1) * X, axis = 0)

D = np.sum(weigth, axis = 0)
print(M, V, D)
