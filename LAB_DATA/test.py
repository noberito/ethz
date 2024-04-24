import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.preprocessing

degree = 4
ndata = 10000
weight_exp = 1
weigth = np.ones(ndata)

print(np.prod([2, 2]))

a = np.array([0, 1, 2, 3, 4])
b = np.zeros(5)
b[0] = 1
p = np.array([1, 1, 1, 1, 1])

X = np.zeros((ndata,5))
X[0] = np.ones(5) * 1
M = np.zeros((5, 5))

polyN = sklearn.preprocessing.PolynomialFeatures((degree, degree), include_bias=False)
X = polyN.fit_transform(X[:, :5],)
powers = polyN.powers_

selected_indices = []
for i in range(len(powers)):
    k,l,m,n,p = powers[i]
    if ((m==n) and (n==p)) or ((m%2==0) and (n%2==0) and (p%2==0)):
        selected_indices.append(i)

powers = powers[selected_indices].T
sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
powers = powers.T[sorted_indices]

X = X[:,sorted_indices]
ndata, nmon = X.shape
M = np.zeros((nmon, nmon))
"""
for i in range(ndata):
    Xi = np.expand_dims(X[i], axis=0)
    M = M + weigth[i] * np.dot(Xi.T, Xi)
M = M / 2

V = np.sum(np.expand_dims(weigth, axis=1) * X, axis = 0)

D = np.sum(weigth, axis = 0)
print(M.shape, V.shape, D.shape)

def J(a):
    return np.dot(np.dot(a, M), a) - np.dot(V, a) + D

def Grad_J(a):
    return 2 * np.dot(M, a) - V


print(J(b))"""