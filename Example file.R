library(MatManlyMix)

data("IMDb")
Y <- IMDb$Y

# We run the analysis by using 14 cores for parallel computing.
# Modify this value if you have/want to use a lower number of cores on your PC.

initMVN <- Pars_init.paral(X = Y, k = 1:5, density = "MVN", nThreads = 14)
resMVN <- Pars_mixt.paral(X = Y, k = 1:5, init.par = initMVN, density = "MVN", nThreads = 14)
winMVN <- extract.bestM(resMVN, density = "MVN")

initMVTIN <- Pars_init.paral(X = Y, k = 1:5, density = "MVTIN", nThreads = 14)
resMVTIN <- Pars_mixt.paral(X = Y, k = 1:5, init.par = initMVTIN, density = "MVTIN", nThreads = 14)
winMVTIN <- extract.bestM(resMVTIN, density = "MVTIN")
