# Efficient constrained Gaussian processes for shape restricted function estimation using elliptical slice sampling (ESS)
This repository contains R codes and functions for implementing the efficient constrained Gaussian processes for shape restricted function estimate using Elliptical Slice Sampling (ESS) based on Maatouk et al. (2024b), and Maatouk and Bay (2017). Furthermore, it includes two functions for generating very large Gaussian vectors extracted from stationary Gaussian processes. The approach developed in this repository can incorporate multiple shape constraints such as monotonicity, boundedness, and convexity, and estimate the hyperparameters.

# Description of the associated R files:
1. 'all_base_functions_1D.R' file contains all the required base functions to generate the efficient constrained GP models, including functions like
Fast.LS (draws a very large Gaussian vector prior based on Maatouk et al. [2023b]), LS.KLE (draws a very large Gaussian vector prior based on Maatouk et al. [2023a]),
samp.WC (draws a sample based on sampling scheme developed in Wood and Chan [1994]), ESS (draws a sample using elliptical slice sampler by Murray et al. [2010]).
One can find detailed description of each functions in the file.
2. 'all_models_1D.R' file contains all models that implemented the efficient constrained Gaussian processes for shape-restricted function estimation using elliptical slice sampling like 'linCGP.ESS' (the proposed highly efficient approach based on the elliptical slice sampling), 'linCGP.HMC' (the approach based on the Hamiltonian Monte Carlo sampler) and 'linCGP.WC.ESS' (the approach based on the Fast Fourier Transform). The proposed function 'linCGP.ESS' can incorporate multiple shape constraints, estimate the hyperparameters and and utilize various efficient samplers to generate the prior. This function has been tested in cases of extreme constraints, a task that becomes impossible for the HMC sampler.
3. 'illustr_exam_1D.R' file contains code for several numerical examples, demonstrating both boundedness and multiple shape constraints (monotonicity, boundedness and convexity).


   For more details on the codes or the functions, refer to the associated R files.

# Note:
Part of this work was conducted with the support of the consortium in Applied Mathematics CIROQUO, gathering partners in technological research (BRGM, CEA, IFPEN, IRSN, Safran, Storengy) and academia (CNRS, Ecole Centrale de Lyon, Mines Saint-Etienne, University of Grenoble, University of Nice, University of Toulouse) around advanced methods for Computer Experiments [Link]( https://doi.org/10.5281/zenodo.65812)

# Authors:
Hassan Maatouk (EBI, Cergy-Pontoise).

# Maintainer: 
Hassan Maatouk, h.maatouk@hubebi.com

# References
Maatouk, H. and Bay, X. (2017). "Gaussian process emulators for computer experiments with inequality constraints". Mathematical Geosciences, 49(5): 557-582 [doi](https://link.springer.com/article/10.1007/s11004-017-9673-2)

Maatouk, H. and Rullière, D. and Bay, X. (2023a). "Sampling large hyperplane‐truncated multivariate normal distributions", Computational Statistics. [doi](https://link.springer.com/article/10.1007/s00180-023-01416-7)

Maatouk, H. and Rullière, D. and Bay, X. (2023b). "Large-scale constrained Gaussian processes for shape-restricted function estimation". [preprint](https://hal.science/hal-04348962/document)


Maatouk, H. and Rullière, D. and Bay, X. (2024a). "Bayesian analysis of constrained Gaussian processes". Bayesian Analysis, doi: 10.1214/24-BA1429
[doi](https://projecteuclid.org/journals/bayesian-analysis/advance-publication/Bayesian-Analysis-of-Constrained-Gaussian-Processes/10.1214/24-BA1429.full)

Maatouk, H. and Rullière, D. and Bay, X. (2024b). "Efficient constrained Gaussian processes using elliptical slice sampling". [preprint](https://hal.science/hal-04496474)

Murray I., Adams R., MacKay D. (2010). "Elliptical slice sampling". In: Proceedings of the thirteenth international conference on artificial intelligence and statistics, JMLR Workshop and Conference Proceedings, pp 541–548

Wood, A. and Chan, G. (1994). "Simulation of Stationary Gaussian Processes in [0,1]^d". Journal of Computational and Graphical Statistics. [doi](https://www.jstor.org/stable/1390903)

Ray, P. and Pati, D. and Bhattacharya, A. (2020). "Efficient Bayesian shape-restricted function estimation with constrained Gaussian process priors". Statistics Computing 30, 839–853. [doi](https://link.springer.com/article/10.1007/s11222-020-09922-0) 
