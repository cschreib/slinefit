# slinefit

```slinefit``` is a simple software that can be used to derive spectroscopic redshifts from 1D spectra, and measure line properties (fluxes, velocity width, velocity offset, ...). Here is the list of the main features:

* Brute force redshift scan on pre-determined interval (```z0=...``` and ```dz=...``` will search redshifts ```z0-dz``` and ```z0+dz```).
* Can fit one or more continuum templates simultaneously with emission and/or absorption lines.
* Any number of emission/absorption lines can be used in the fit (can be zero, or just one, or all the known lines, or any combination of lines).
* The best solution and probability distribution are determined from the chi2.
* Lines velocity widths and offsets can be varied, either on a fixed grid, or using a Levenberg-Markwardt solver.
* Line widths can be constrained to be the same for all lines, or to vary independently.
* Accurate uncertainties on line properties (flux, width, offset) can be derived using Monte Carlo simulations.
* Can use multithreading to speed up the Monte Carlo simulations.
* Can output results in FITS and/or ASCII format.

On my desktop computer, the code takes 7 seconds to find the redshift of a galaxy with one continnuum template and five emission lines. This is using a 1D spectrum with about 1000 spectral elements (R~3000), with 8000 possible redshifts and 2 possible line widths.


# Install instructions

To install ```slinefit```, you first need to build it. For this, you will need a recent C++ compiler, [git](https://git-scm.com/) and [CMake](https://cmake.org/). All of these tools are available by default on most modern linux distributions, and for MacOS you may need to install them with MacPorts. From there installing ```slinefit``` is very easy. Navigate to a directory of your choosing where the code will be downloaded (it will be placed in a subdirectory called ```slinefit```), then execute the following commands:
```
# Download the code and dependencies
git clone https://github.com/cschreib/slinefit.git
cd slinefit

# Compile
mkdir build
cd build
cmake ../
make install
```

This will create an executable called ```slinefit``` in the ```slinefit/bin``` directory, which you can use immediately.
