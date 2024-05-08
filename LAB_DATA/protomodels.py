import numpy as np

def mises(sigma):
    """
    Vectorized Mises function
    """
    sigma = np.atleast_2d(sigma)
    if sigma.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    s11 = sigma[:, 0]
    s22 = sigma[:, 1]
    s33 = sigma[:, 2]
    s12 = sigma[:, 3]
    s13 = sigma[:, 4]
    s23 = sigma[:, 5]

    res = np.sqrt(0.5 * (s22 - s33)**2 + 0.5 * (s33 - s11)**2 + 0.5 * (s11 - s22)**2 + 3 * (s23**2 + s13**2 + s12**2))
    return res

def tresca(sigma):
    sigma = np.atleast_2d(sigma)
    N = len(sigma)
    if sigma.shape[-1] != 6:
        raise ValueError("Input array must have shape (N, 6)")

    sigma_matrix = np.zeros((N, 3, 3))
    
    index_x = np.array([0, 1, 2, 0, 0, 1])
    index_y = np.array([0, 1, 2, 1, 2, 2])
    sigma_matrix[:, index_x, index_y] = sigma
    sigma_matrix = sigma_matrix + np.transpose(sigma_matrix, (0, 2, 1))
    sigma_matrix[:,np.arange(3),np.arange(3)] = sigma_matrix[:,np.arange(3),np.arange(3)] // 2

    ps = np.linalg.eigvals(sigma_matrix)
    dps = np.zeros((N, 3))

    dps[:,0] = np.abs(ps[:, 0] - ps[:, 1])
    dps[:,1] = np.abs(ps[:, 0] - ps[:, 2])
    dps[:,2] = np.abs(ps[:, 1] - ps[:, 2])

    return(np.max(dps, axis=1))