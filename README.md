# pdphmc
Piecewise deterministic processes using Hamiltonian dynamics

**This package has been discontinued - the current development of numerical generalized randomized HMC is at [https://github.com/torekleppe/pdmphmc]( https://github.com/torekleppe/pdmphmc), which constitute an almost complete rewrite of the methology**

This is the development version of the R-package pdphmc, which is used in the paper "Connecting the Dots: Towards Continuous Time Hamiltonian Monte Carlo" by Tore Selland Kleppe [https://arxiv.org/abs/2005.01285](https://arxiv.org/abs/2005.01285)


See also [https://github.com/torekleppe/PDPHMCpaperCode]( https://github.com/torekleppe/PDPHMCpaperCode) for the model- and data files, and the specific version of the package used in the above paper.

The R-package has currently been tested on mac (using the clang compiler and R 4.x), linux (using native gcc and R 4.x) and windows 10 (R 4.x and using the mingw64 compiler that ships with Rtools4.0). In general, it requires that you have a working installation of R>=4.0.0 and **a working installation of the R-package rstan**.

Installation of the package is most easily done using the `devtools` function `install_github()`, e.g.: ```devtools::install_github("https://github.com/torekleppe/pdphmc")```



