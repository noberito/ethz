# Framework polyN : Parameters guide

### Materials parameters

**material** (str) : Name of the material given

**enu** (float) : Poisson ratio of the material given

**density** (float) : Density of the material given

**degree** (int) : Degree of the polyN function (not every degree available lack of starting points)

### First optimization based on yield stress points and r-values

**nt6** (0 or 1) : Add the nt6 points at each iteration of the first optimization

**sh** (0 or 1) : Add the shear points to the database for first optimization

**weight_ut** (float in [0, 1]) : Weight of the stress points against r-value points in 1st optimization

**weight_e2** (float in [0, 1]) : Weight of NT6 and SH points against UT point

**weight_vir** (float in [0, 1]) : Weight of virtual points against UT point

**nb_virtual_pt** (int) : Number of points extracted from a protomodel to complete the dataset

**protomodel** (str) : Protomodel used (only mises available)

**law** (string) : Hardening law (framework only tested with Swift (not Voce nor Swift-Voce))

**gen_v_data** (0 or 1) : Generating virtual data (must be set to 1 after changing the number of virtual points otherwise set it to 0 not to lose time)

**gen_e_data** (0 or 1) : Generating calibration data (must be set to 1 after adding a material otherwise set it to 0 not to lose time)

### Second optimization using FEA

**n_opti** (int) : Number of iterations maximum of the gradient descent in the framework

**input_type** (str in ["UMAT", "VUMAT"]) : Type of abaqus analysis (only works for UMAT)

**var_optim** (int separated by a comma !, ex : "1,4,5" not "1, 4,5")) : Variables of the yield function to optimize (starts at 1, 0 is the last coefficient, etc)

opti = 1
loadcoeff = 0
