#####################################################################
# This is a readme for the locproj.m, locproj_conf.m, and locproj_cv.m matlab programs
#####################################################################

A demo implementing these programs can be found in demo.m


######################################################################
file: locproj.m
######################################################################

locproj creates a struct object. Once you create the object, you can see it in the workspace. 
You can see all the properties the object contains by clicking on it. 

obj = locproj(y,x,w,H_min,H_max,type)
obj = locproj(y,x,w,H_min,H_max,type,r,lambda)

Description
obj = locproj(y,x,w,H_min,H_max,type) returns a local projection struct with the 
instrumental impulse response coefficients for a normal local projection.

obj = locproj(y,x,w,H_min,H_max,type,r,lambda) returns a local projection struct with
instrumental impulse response coefficients for a smooth local projection.

Input Arguments
y - Response (LHS) vector
x - Endogenous impulse vector
w - Optional matrix of all other RHS vectors, including the constant, instrumental variables, and lagged variables (if applicable)
H_min - First period the shock hits (0 or 1)
H_max - Number of periods the local projection spans
type - tring value that dictates whether the local projection ran is regular or smooth ('reg' or 'smooth')
r - Used in smooth local projection: order of the limit polynomial  (for example, r=2 implies the impulse response is shrunk towards a line)
lambda - Used in smooth local projection: numeric value for the shrinking parameter lambda

Output Arguments
obj.T - Number of observations
obj.H_min - First period the shock hits (0 or 1)
obj.H_max - Number of periods the local projection spans
obj.HR - Number of periods of impulse responses (H_max + 1 - H_min)
obj.K - Number of columns of the basis spline
obj.B - Basis spline used in smooth local projection
obj.P - Penalty matrix
obj.lambda - lambda used above as an input
obj.type - Type of local project, smooth or reg
obj.delta - numeric value containing the residual standard error from the first stage OLS regression of the endogenous impulse variable on other endogenous and exogenous RHS variables
obj.idx - Observation horizon matrix
obj.Y - Shifted response (LHS) vector ran in local projection simulatenously over all horizons
obj.X - Shifted RHS matrix ran in local projection simultaenously over all horizons (includes all RHS vectors including the constant)
obj.theta - Theta matrix containing coefficient estimates for all variables, for each value of labmda
obj.IR - matrix containing estimations of the instrumented impulse coefficient for each value of lambda


######################################################################
file: locproj_cv.m
######################################################################

locproj_cv

obj = locproj_cv(obj, K, lambda)

Description
locproj_conf takes in a struct object (a previously created smooth local projection),
and runs a cross validation for the optimal choice of the shrinking parameter lambda.

Input Arguments
obj - A previously created smooth local projection struct
K - The size of the cross-validation sample, (we choose K=5).
lambda - Options for the shrinking parameter lambda, can be a single value or a vector of values.

Output Arguments
obj.rss - Residual sum of squares for each choice of lambda
obj.lambda - lambda input

######################################################################
file: locproj_conf.m
######################################################################

locproj_conf

obj = locproj_conf(obj, H)
obj = locproj_conf(obj, H, lambda)

Description
locproj_conf takes in a struct object (a previously created local projection),
and creates confidence intervals around the impulse response coefficients.

Input Arguments
obj - A previously created local projection struct
H - Number of horizons in the local projection
lambda - Used in smooth local projection: numeric value for the shrinking parameter lambda

Output Arguments
obj.ster - Vector of standard errors for the impulse response in each period
obj.conf - Matrix of 95 percent confidence interval values around the impulse response in each period
