## The optimization algorithm

import numpy as np

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

def readData(data):
    pass

n = 10
data=np.array([[4,5,1,1,1], [4,3,1,1,1], [4,1,1,1,1], [4,2,1,1,1]])

# Sample data (shape (ndata, nvariables))

def format_data(data, n):
    """
        Format data in order to be optimized
        Output :
            - data_formatted : n_data * n_monomials
    """
    
    #y = np.ones(len(X))
    # Define the degree of the polynomial and the number of polynomials of degree equal to n
    poly = PolynomialFeatures(degree=(n,n), include_bias=False)
    X_poly = poly.fit_transform(data)
    nMonN = len(poly.powers_)
    selected_indices = []

    # Polynomials that satisfy the constraints on shear components
    for i in range(nMonN):
        powers = poly.powers_[i]
        [k,l,m,n,p] = powers
        if ((m%2==0) and (n%2==0) and (p%2==0)) or ((m==n) and (n==p)):
            selected_indices.append(i)

    X_powers = poly.powers_[selected_indices]
    X_poly_selected = X_poly[:, selected_indices]

    #Sort matrixes so that it corresponds to the article
    X_powers_t = X_powers.transpose()
    sorted_indices = np.lexsort((X_powers_t[1],X_powers_t[2], X_powers_t[3], X_powers_t[4]))

    X_powers = X_powers[sorted_indices]
    X_poly_selected = X_poly_selected[:,sorted_indices]

    return(X_poly_selected)

format_data(data, n)

## Data is well formatted with X_poly_selected_data[i,j] = phi_j(sigma[i])



"""# Fit the model to the data
model.fit(X, y)

# Predict using the model
predictions = model.predict(X)

# Print the coefficients and intercept
print("Coefficients:", model.named_steps['linearregression'].coef_)
print("Intercept:", model.named_steps['linearregression'].intercept_)
print("Predictions:", predictions)"""
