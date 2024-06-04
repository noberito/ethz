
**MODELLING AND CALIBRATING A POLYN YIELD SURFACE FUNCTION**

---

In the context of my research assistant internship at ETH Zurich, I have to implement and find an operationnal strategy to calibrate a yield surface function based on the polyN model which is a homogenous polynomial of degree N.

To make the calibration operationnal please copy paste a folder "{material}_results" containing the csv files for a given material in the lab_data directory and then change the parameters you want in the polyN_optimization/param.txt. Then run polyN_optimization/optimize_polyN.py. The coefficients of the polyN function will be accessible in the polyN_coeff.npy file.
