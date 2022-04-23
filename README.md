# Parsimonious-MVTIN-mixtures
This repository contains the code for fitting 196 parsimonious mixtures of MVTIN and MVN distributions.
In the following you can find a description of the functions (and their arguments) contained in the Functions.R file.

##### Pars_init.paral #####

Runs the initialization of the EM-based algorithms used for fitting parsimonious mixtures of "MVN" or "MVTIN" distributions. Parallel computing is implemented and highly recommended for a faster calculation.

X : An array of dimensions p x r x N, where p and r are the variables in the rows and columns, respectively, while N refers to the statistical units.

k : A positive integer or a vector containing the numbers of groups to be tried.

nstartR : An positive integer indicating the number of random starts.

nThreads: A positive integer indicating the number of cores used for running in parallel.

density : The matrix-variate distribution to be used for the mixture model. Possible values are: "MVN" for the normal distribution, "MVTIN" for the tail-inflated normal distribution.

##### Pars_mixt.paral ##### 

Fits, by using EM-based algorithms, parsimonious mixtures of "MVN" or "MVTIN" distributions to the given data. Parallel computing is implemented and highly recommended for a faster model fitting. 

X : An array of dimensions p x r x N, where p and r are the variables in the rows and columns, respectively, while N refers to the statistical units.

k : Must be the same used for the Pars_init.paral() function.

init.par : The output of the Pars_init.paral() function for initializing the fitting algorithm.

mod.row : A character vector indicating the parsimonious structure of the row covariance (or scale) matrices. Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV" or "all". When "all" is used, all of the 14 parsimonious structures are considered.

mod.col : A character vector indicating the parsimonious structure of the column covariance (or scale) matrices. Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV" or "all". When "all" is used, all of the 7 parsimonious structures are considered.

mod.theta : A character vector indicating the parsimonious structure of the tailedness parameters. Possible values are: "E", "V" or "all". When "all" is used, both parsimonious structures are considered.

nThreads: A positive integer indicating the number of cores used for running in parallel.

density : The same density used in the Pars_init.paral() function.

##### r_data.gen ##### 

Random number generation for mixtures of "MVN" or "MVTIN" distributions.

num : A positive integer indicating the sample size of the data to be generated.

pi : A numeric vector of length k representing the probability of belonging to the k groups for each data point.

M : An array for the mean matrices of dimensions p x r x k, where p and r are the variables in the rows and columns, respectively, while k refers to the number of groups.

U : An array for the row covariance (or scale) matrices of dimensions p x p x k, where p identifies the variables in the rows and k refers to the number of groups.

V : An array for the column covariance (or scale) matrices of dimensions r x r x k, where r identifies the variables in the columns and k refers to the number of groups.

theta : A vector of length k representing the tailedness parameters.

density : The matrix-variate distribution to be used for the mixture model. Possible values are: "MVN" for the normal distribution, "MVTIN" for the tail-inflated normal distribution.

##### extract.bestM ##### 

This functions extracts the best fitting model according to the Bayesian information criterion (BIC).

results : The output of the Pars_mixt.paral() function.

mod.row : The same character vector used in the Pars_mixt.paral() function.

mod.col : The same character vector used in the Pars_mixt.paral() function.

mod.theta : The same character vector used in the Pars_mixt.paral() function.

density : The same density used in the Pars_mixt.paral() function.


