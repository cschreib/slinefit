# slinefit

[![Build Status](https://travis-ci.com/cschreib/slinefit.svg?branch=master)](https://travis-ci.com/cschreib/slinefit)

<!-- MarkdownTOC depth=0 -->

- [Description](#description)
- [Install instructions](#install-instructions)
- [Format of input spectrum](#format-of-the-input-spectrum)
- [Typical usage](#typical-usage)
- [In depth tutorial](#in-depth-tutorial)

<!-- /MarkdownTOC -->

# Description

```slinefit``` is a simple software that can be used to derive spectroscopic redshifts from 1D spectra, and measure line properties (fluxes, velocity width, velocity offset, ...). Here is the list of the main features and facts:

* Reads the 1D flux and error spectra from a FITS file, wavelengths from the WCS data.
* Frequencies are supported too!
* Performs a brute-force redshift scan over a fixed interval.
* Lines can be chosen from a pre-defined list; more can be added from the command line.
* Lines can be emission or absorption lines.
* Lines have Gaussian velocity profiles, formally integrated over each spectral element.
* Line velocity widths can be varied, optionally constrained to be the same for all lines.
* Lines can optionally have variable velocity offsets.
* Lines can optionally have multiple velocity components, with different widths and offsets.
* One or more continuum templates can be fit simultaneously with the lines.
* The best redshift and its probability distribution are determined from the chi2.
* Accurate uncertainties on line properties can be derived using Monte Carlo simulations.
* Can use multi-threading to speed up the Monte Carlo simulations.
* Outputs redshift, line fluxes, widths, offsets and rest-frame equivalent widths.
* Can output results in both FITS and ASCII format.
* Has been used to fit spectra from VIMOS, KMOS, X-SHOOTER, MOSFIRE, ALMA and VLA.


# Install instructions

To install slinefit, you first need to build it. For this, you will need a recent C++ compiler, [git](https://git-scm.com/) and [CMake](https://cmake.org/). All of these tools are available by default on most modern linux distributions, and for MacOS you may need to install them with MacPorts. From there installing slinefit is very easy. Navigate to a directory of your choosing where the code will be downloaded (it will be placed in a subdirectory called ```slinefit```), then execute the following commands:
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

This will create an executable called ```slinefit``` in the ```slinefit/bin``` directory, which you can use immediately. In the ```slinefit/bin/templates``` directory, you will find a set of stellar continuum templates that can be used to fit UV-NIR spectra of galaxies. These templates correspond to the template set of EAzY (Brammer et al. 2008), for which I increased the spectral resolution by fitting them with FAST (Kriek et al. 2009).


# Format of input spectrum

Because spectroscopy is a precision measurement, it is crucial that the input data are laid out in a clear, unambiguous way. If you care about your results, please read the following carefully.

Because precision and clarify are paramount, slinefit only supports one type of input file: multi-extension FITS files, with wavelength/frequency information precisely defined in the WCS coordinate system. A number of instrument pipelines already use this format as their default output format for 1D spectra, so with luck you will not have to do any conversion. But in case this is needed, detailed instructions are provided below.

The multi-extension FITS file is expected to contain at least two "HDU" (header data units) containing, respectively, the flux and the error spectra. This can be, for example, the flux data in the primary HDU, and the error data in the first extension. Or the flux in the first extension, and the error in the second extension, etc. The code will not attempt to automatically determine this; you have to specify where the data is located in the command line options. This is done with the ```flux_hdu``` and ```error_hdu``` options. For example, setting ```flux_hdu=0``` tells the code to look for the flux data in the primary HDU. The default is ```flux_hdu=1``` and ```error_hdu=2```, so both flux and error are placed in FITS extensions with an empty primary HDU (this is the standard ESO format).

For the code to return accurate results, both the flux and error spectra must have the same unit. The error spectrum is also assumed to give the "1 sigma" uncertainty of each spectral element (not the variance). Missing data must be indicated with, either, zero or negative uncertainty, infinite or NaN uncertainty, or infinite or NaN flux. It is not necessary to perform a continuum subtraction before running the fit, as slinefit can model the continuum emission when needed. Lastly, if you are only interested in measuring a redshift and you do not care about the line fluxes, then it is possible (although not recommended) to fit spectra that are not flux-calibrated .

The wavelength or frequency axis is determined from the WCS keywords in the FITS header of the flux extension. There are two mandatory keywords.

First, ```CUNIT1``` specifies the unit of the axis, which can be any of the following (case insensitive):

* wavelengths: angstrom, nm, um, micron, cm, m
* frequencies: hz, khz, mhz, ghz

Any non-recognized unit will trigger an error. If the ```CUNIT1``` keyword is missing, a warning will be issued and the code will assume the unit is ```um```.

Second, ```CTYPE1``` indicates the format of the axis. There are three formats: linear, logarithmic, and tabulated. If ```CTYPE1``` starts with ```LOG...```, the code will assume a logarithmic scaling. If ```CTYPE1``` starts with ```TAB...```, the code will assume a tabulated axis. In all other cases, the code will assume linear scaling (in which case, it is recommended but not mandatory to set ```CTYPE1='WAVE'``` for a wavelength axis and ```CTYPE1='FREQ'``` for frequency data).

For both the linear and logarithmic formats, the code will then use the standard ```CRPIX1```, ```CRVAL1```, and ```CDELT1``` keywords. If ```i``` is the index of the spectral element in the data (starting at 1 for the first element, as per the FITS convention), the following equations are used to compute the spectral coordinate at the _center_ of the spectral element:

```
# For a linear axis:
axis = (i - CRPIX1)*CDELT1 + CRVAL1
# For a logarithmic axis:
axis = 10^((i - CRPIX1)*CDELT1 + CRVAL1)
```

For example, for a spectrum of 100 elements with linear coverage from 450nm to 650nm, you would set ```CTYPE1='WAVE'```, ```CUNIT1='nm'```, ```CRPIX1=1```, ```CRVAL1=450```, ```CDELT1=2```.

The tabulated format should be used for axes that are neither linear not logarithmic. It is a very explicit format: the code attributes a range of wavelength/frequency to each spectral element, and you must provide the beginning and end of this range for each element. This allows for example to have disjoint spectral elements (e.g, gaps), or even overlapping elements. To use this format, you must set the ```CTYPE1``` keyword with the format ```'TABxx'```, where ```xx``` can be any one or two letters (e.g., ```TABW```). You must then specify two keywords: ```xxLOWEXT``` and ```xxUPEXT``` (e.g., ```WLOWEXT``` and ```WUPEXT``` if ```xx``` is ```W```), which specify the FITS extensions in which the code will read the beginning and end of the spectral range of each spectral elements, respectively. These extensions must contain 1D data with as many elements as the flux spectrum, and the unit is set by the ```CUNIT1``` keyword in the flux header.

Lastly, it is possible to fit multiple spectra at once. These spectra can cover either the same wavelength/frequency range, or a completely different range. The axes can be any combination of wavelength or frequency axes, with any combination of units. However the fluxes and uncertainties must have the same units. To do this, simply create a text file (extension ```.txt``` or ```.cat``` or ```.dat```, anything but ```.fits```) containing the filename (or full path) to each spectrum on a separate line, and give this text file as input spectrum to slinefit.


# Typical usage

Usually you will want to run a "first pass" to measure the redshift, which is pretty fast. For example, on my desktop computer, the code takes 7 seconds to find the redshift of a galaxy with one continuum template (typically the best-fit template obtained from a fit to the photometry with your favorite photometric redshift code) and five emission lines. This is using a 1D spectrum with about 1000 spectral elements (R~3000), scanning 8000 possible redshifts between z=2 and z=5 with 2 possible line widths, and without Monte Carlo simulations. This is the command that was used:
```
slinefit spectrum.fits z0=3.5 dz=1.5 delta_z=1 width_min=60 width_max=300 delta_width=2 \
    same_width use_global_chi2 full_range fit_continuum_template verbose \
    lines=[em_halpha,em_n2_6583,em_o3_5007,em_hbeta,em_o2_3727]

# Notes:
# - 'z0=3.5 dz=1.5' means we will scan redshifts from 2 to 5.
#
# - 'delta_z' and 'delta_width' are specified in terms of resolution elements. So 'delta_z=1'
#   means that we will build the redshift grid with a step such that an emission line would
#   move by roughly one spectral element from one redshift to the next.
#
# - 'use_global_chi2' and 'full_range' force the fit to use the entire spectrum to compute
#   the chi2 and fit the continuum, not only the regions surrounding potential lines. This
#   is good for a redshift search.
#
# - 'fit_continuum_template' allows the program to fit one ore more templates for the
#   continuum level, the templates must be in the "templates" sub-folder with *.dat extension.
#
# - 'lines' specifies which lines to use, here Halpha, [NII], [OIII], Hbeta and [OII]. These
#   are the brightest lines, and will be enough to find the redshift if some lines are indeed
#   present in the spectrum.
```

Once the redshift is known, you will want to do a "second pass" where line properties are varied on a finer grid, where fainter lines are added, or where the continuum is modeled with more templates. For example, on the same computer, it takes about 2 minutes to measure the properties of 14 lines, with 6 continuum templates, 25 different redshifts (```dz=0.002```, centered on the known zspec), 9 possible line widths (with a step of 50 km/s), and 41 velocity offsets (step of 50 km/s). The most time-consuming step is the Monte Carlo simulation (performing the fit itself takes about one second).
```
slinefit spectrum.fits z0=2.369 dz=0.002 delta_z=0.5 width_min=50 width_max=500 delta_width=0.5 \
    same_width allow_offsets offset_max=1000 delta_offset=0.5 fit_continuum_template verbose \
    num_mc=200 threads=8

# Notes:
# - Here we have decreased 'delta_z' and 'delta_width', so the grid became finer.
#
# - We have reduced 'dz' because we know the redshift, we just want to refine it.
#
# - We have increased the range of possible line widths using 'width_min' and 'width_max'.
#
# - We now 'allow_offsets', so bright lines can be shifted in velocity up to 'offset_max'.
#
# - The list of emission lines has been omitted, so the program will try to fit all the lines
#   that it knows about (run "slinefit list_lines" to see the list).
#
# - We have enabled Monte Carlo simulations with 'num_mc=200', so the program will perturb
#   the input spectrum 200 times and refit the lines each time.
#
# - We allowed the program to use up to 8 concurrent threads with the 'threads' option.
```

More detailed information can be obtained by calling ```slinefit``` without arguments. There are a lots of small tweaks you can apply using command line options, do take a look!

Below are two examples of fits obtained with ```slinefit```. The first is a fit of the [OII] line on top of strong continuum with Balmer absorption (the absorption lines are part of the continuum template, which is was created with the Bruzual & Charlot 2003 models). The second is a rest-FUV spectrum with the CII] emission line and several FeII absorption lines.

![Demo OII](demo_oii.png) ![Demo UV](demo_uv.png)


# In depth tutorial

A tutorial is available in the ```docs``` directory, where the various features of the code are explained and illustrated on a real spectrum.
