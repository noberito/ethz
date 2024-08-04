import pandas as pd
import numpy as np
import os

polyN_dir = os.path.dirname(os.path.abspath(__file__))
sep = os.sep

def hill48(S, coeff):
    if S.ndim == 1:
        S = np.expand_dims(S, axis=0)

    p12 = coeff[0]
    p22 = coeff[1]
    p33 = coeff[2]

    P = np.zeros((3,3))

    P[0,0] = 1
    P[0,1] = p12
    P[1,0] = p12
    P[1,1] = p22
    P[2,2] = p33

    PS = np.einsum("ik,jk->ji", P, S)
    PSS = np.einsum("ik, ik -> i", PS, S)

    PSS = np.expand_dims(PSS, axis=1)

    return(np.sqrt(PSS))

def grad_hill48(S, coeff):
    if S.ndim == 1:
        S = np.expand_dims(S, axis=0)
    
    g12 = coeff[3]
    g22 = coeff[4]
    g33 = coeff[5]

    G = np.zeros((3,3))

    G[0,0] = 1
    G[0,1] = g12
    G[1,0] = g12
    G[1,1] = g22
    G[2,2] = g33

    print(G)
    GS = np.einsum("ik,jk->ji", G, S)
    GSS = np.einsum("ik, ik -> i", GS, S)

    h = np.sqrt(GSS)
    h = np.expand_dims(h, axis=1)

    grad_f = 2 * GS
    grad_h = (1/2) * grad_f / h

    return(grad_h)
    
def load_coeff_hill48(material):
    filename = "_Hill48naFr_pre.csv"
    filepath = polyN_dir + sep + "results_exp" + sep + material + sep + filename
    try :
        df = pd.read_csv(filepath)
        coeff = df.loc[3].values
        return coeff
    except:
        print("YLD2000 not calibrated")

def test():
    material = "DP780"
    coeff = load_coeff_hill48(material)
    S = np.zeros(3)
    S[2] = 2
    print(grad_hill48(S, coeff))
