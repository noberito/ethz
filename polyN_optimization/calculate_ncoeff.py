import sklearn.preprocessing
import numpy as np

def get_param_polyN_mini(degree):
    """
        Returns polyN minimalistic powers for the given degree

        Input :
            - degree : integer, degree of the polyN function

        Output :
            - powers : ndarray of shape (nmon, 5), powers[i, j] is the power of the variable j in the monomial i
    """
    x0 = np.zeros((2,3))
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X = polyN.fit_transform(x0)
    powers = polyN.powers_
    nmon_degree = len(powers)

    selected_indices = []
    for i in range(nmon_degree):
        k,l,m = powers[i]
        if (k + l + m == degree) and (m%2==0) :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2]))
    powers = powers.T[sorted_indices]
    X = X[:,sorted_indices]

    return(powers)

def get_param_polyN(degree):
    """
        Returns the parameters of polyN according to the degree
        Input :
            - degree : integer, degree of the polyN function
        Output :
            - powers : ndarray of shape (nmon, 5), powers[i, j] is the power of the variable j in the monomial i
    """
    x0 = np.zeros((2,6))
    polyN = sklearn.preprocessing.PolynomialFeatures(degree, include_bias=False)
    X = polyN.fit_transform(x0)
    powers = polyN.powers_
    nmon_degree = len(powers)

    selected_indices = []
    for i in range(nmon_degree):
        k,l,m,n,p,q = powers[i]
        if (k + l + m + n + p + q == degree)  :
            selected_indices.append(i)

    powers = powers[selected_indices].T
    X = X[:,selected_indices]

    sorted_indices = np.lexsort((powers[1], powers[2], powers[3], powers[4]))
    powers = powers.T[sorted_indices]
    
    X = X[:,sorted_indices]


    return(powers)


for d in [2,4,6,8]:
    print(len(get_param_polyN_mini(d)))

for d in [2,4,6,8]:
    print(len(get_param_polyN(d)))
