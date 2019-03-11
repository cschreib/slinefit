# slinefit tutorial

This tutorial will guide you through the analysis of an example spectrum observed with VIMOS on the VLT (taken from the VANDELS DR2). We will get a first fit running, and then see how it can be improved by using a number of command line options.

<!-- MarkdownTOC autolink="true" -->

- [Setting up slinefit](#setting-up-slinefit)
- [Setting up the initial run](#setting-up-the-initial-run)
    - [The input file, and FITS extensions](#the-input-file-and-fits-extensions)
    - [The redshift and line grids](#the-redshift-and-line-grids)
    - [Fixing the FITS header](#fixing-the-fits-header)
    - [A first \(not very good\) run](#a-first-not-very-good-run)
- [Improving the fit](#improving-the-fit)
    - [Visualizing the fit](#visualizing-the-fit)
    - [Cutting corrupted edges](#cutting-corrupted-edges)
    - [Fitting the continuum](#fitting-the-continuum)

<!-- /MarkdownTOC -->


# Setting up slinefit

To install slinefit, please follow the instructions given in the README. In what follows, I will assume you have installed slinefit in the ```/home/user/programs/slinefit/``` directory, but you can of course install the program anywhere else, provided you update the path to the program in the instructions below. Consequently, the program slinefit is located in ```/home/user/programs/slinefit/bin/slinefit```.


# Setting up the initial run

## The input file, and FITS extensions

In this directory, you will find a FITS file called ```sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits```, which was taken straight from the VANDELS DR2 data base. To run slinefit on this spectrum, we must supply the name of the FITS file as first argument to the program:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits
```

Do not attempt to run this yet, however. We have to set a number of important command line options for the fit to proceed correctly.

Firt, this file contains a number of FITS extensions each containing different data. For more information on the format of the FITS files that slinefit can deal with, please look at the README file. In particular, this file contains the flux in the primary FITS HDU, and the uncertainty (or noise, or error) spectrum in the third extension. But slinefit assumes a different format by default, so we have to let it know where to look for the data. This is done with the ```flux_hdu``` and ```error_hdu``` options:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits \
    flux_hdu=0 error_hdu=3
```

## The redshift and line grids

If you try to run the command above, you will see that the program fails with several error messages:

```bash
error: please provide the fiducial redshift z0=...
error: please provide the name of an emission line(s) to fit with lines=[...]
note: available lines:
...
```

Indeed, the way slinefit works is that you must provide it with a first guess of the redshift (```z0```), a range of redshifts to explore around this value (```dz```), and a set of emission lines to look for (```lines```). Let's do this step by step.

This VANDELS spectrum already has a redshift estimate in the header (```HIERARCH PND Z``` keyword), which is ```z=2.8895```. In this tutorial we will ignore this redshift and try to determine our own estimate. Suppose we had a photometric redshift for this object that gave us ```z=2.9 +/- 0.2```. We will therefore use ```z0=2.9``` as a starting point, and to make sure we cover all the possibilities allowed by the photometric redshift we will scan a range of ```dz=0.6``` around this value. This will scan redshifts from ```2.3``` to ```3.5```. Therefore we now have the command:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6
```

By default, the program will scan a fine grid of redshifts: the default redshift step is such that a line would move by just one fifth of a pixel (or spectral element). This is controlled by the ```delta_z``` keyword. Such a small step is good for getting the most precise measurement of the redshift, but at this stage we just want to have a first estimate. We will therefore scan with a step of two pixels (given the spectral resolution of the spectrum, ```R=3000```, this is roughly a step of ```0.0003``` in redshift). We therefore set ```delta_z=1```:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1
```

Now we need to define which emission line to look for. This is done by listing the lines in the ```lines=[...]``` option, where the ```...``` in brackets must be a comma-separated list (without spaces) of emission line names. The list of all emission line name available by default in the program can be seen by calling the slinefit program without any argument:

```bash
/home/user/programs/slinefit/bin/slinefit
```

This will display the command line help for the program. Among other things, it will also print the list of all the supported lines:

```bash
Available lines:
  - abs_al2_1671, lambda=0.16708 (Aluminum II  -- AlII)
  - abs_al3_1855, lambda={0.18547, 0.18628}, ratios={1, 1} (Aluminum III -- AlIII (doublet))
  - abs_c2_1335, lambda=0.13345 (Carbon II -- CII)
  - abs_c3_1176, lambda=0.11755 (Carbon III -- CIII)
  - abs_fe2_1608, lambda={0.16085, 0.16112}, ratios={1, 1} (Iron II -- FeII (doublet))
  - abs_fe2_2344, lambda=0.23442 (Iron II -- FeII)
  - abs_fe2_2380, lambda={0.23745, 0.23828}, ratios={1, 1} (Iron II -- FeII (doublet))
  - abs_fe2_2600, lambda={0.25867, 0.26002}, ratios={1, 1} (Iron II -- FeII (doublet))
  - abs_o1_1302, lambda=0.13022 (Oxygen I -- OI)
  - abs_o4_1342, lambda=0.13416 (Oxygen IV -- OIV)
  - abs_s5_1502, lambda=0.15018 (Sulfur V -- SV)
  - abs_si2_1260, lambda=0.12604 (Silicium II -- SiII)
  - abs_si2_1304, lambda=0.13043 (Silicium II -- SiII)
  - abs_si2_1526, lambda=0.15267 (Silicium II -- SiII)
  - em_c1_370, lambda=370.42 (Carbon I -- CI)
  - em_c1_609, lambda=609.14 (Carbon I -- CI)
  - em_c2_157, lambda=157.741 (Carbon II --- [CII])
  - em_c2_2326, lambda={0.23235, 0.23247}, ratios={1, 1} (Carbon II -- CII] (doublet))
  - em_c3_1909, lambda=0.19087 (Carbon III -- CIII])
  - em_c4_1550, lambda={0.15482, 0.15508}, ratios={1, 1} (Carbon IV -- CIV (doublet))
  - em_co10, lambda=2600.75 (Carbon Monoxyde -- CO(1-0))
  - em_co21, lambda=1300.4 (Carbon Monoxyde -- CO(2-1))
  - em_co32, lambda=866.96 (Carbon Monoxyde -- CO(3-2))
  - em_co43, lambda=650.25 (Carbon Monoxyde -- CO(4-3))
  - em_co54, lambda=520.23 (Carbon Monoxyde -- CO(5-4))
  - em_co65, lambda=433.57 (Carbon Monoxyde -- CO(6-5))
  - em_co76, lambda=371.65 (Carbon Monoxyde -- CO(7-6))
  - em_co87, lambda=325.23 (Carbon Monoxyde -- CO(8-7))
  - em_co98, lambda=298.12 (Carbon Monoxyde -- CO(9-8))
  - em_halpha, lambda=0.65628 (Hydrogen alpha)
  - em_hbeta, lambda=0.48613 (Hydrogen beta)
  - em_hdelta, lambda=0.41017 (Hydrogen delta)
  - em_he2_1640, lambda=0.16404 (Helium II -- HeII)
  - em_hgamma, lambda=0.43405 (Hydrogen gamma)
  - em_lyalpha, lambda=0.12157 (Lyman alpha)
  - em_mg2_2799, lambda={0.27964, 0.28035}, ratios={1, 1} (Magnesium II  -- MgII (doublet))
  - em_n2_6583, lambda={0.65835, 0.6548}, ratios={1, 0.3} (Nitrogen II -- [NII] (doublet))
  - em_n5_1240, lambda={0.12388, 0.12428}, ratios={1, 1} (Nitrogen V -- NV (doublet))
  - em_ne3_3869, lambda=0.38688 (Neon III -- NeIII)
  - em_ne4_2422, lambda=0.24218 (Neon IV -- NeIV)
  - em_ne5_3426, lambda=0.34259 (Neon V -- NeV)
  - em_o1_145, lambda=145.525 (Oxygen I --- [OI])
  - em_o2_3727, lambda={0.3726, 0.37288}, ratios={1, 1} (Oxygen II -- [OII] (doublet))
  - em_o3_5007, lambda={0.50068, 0.49589}, ratios={1, 0.3} (Oxygen III -- [OIII] (doublet))
  - em_palpha, lambda=1.87513 (Pashen alpha)
  - em_s2_6717, lambda={0.67164, 0.67308}, ratios={1, 0.75} (Sulfur II -- [SII] (doublet))
  - em_si4_1400, lambda={0.13938, 0.14028}, ratios={1, 1} (Silicium IV -- SiIV (doublet))
```

The name of the line is the first element of each item on this list. For example ```em_o2_3727``` is the code name for the [OII] emission line (```em``` is for lines that are typically found in emission, and ```abs``` for absorption, but this is only indicative). In the line list used for the fit, we can specify as many lines as we wish. The program is smart enough to detect when a line is not covered, and it will simply ignore it. It would therefore be tempting to include all the lines, and let the code figure out which ones are indeed detected in the spectrum. This is fine if we already know the redshift. However, for a redshift search this can be dangerous: the program has no knowledge of which line is usually brighter or fainter. Therefore, if you specify too many possible lines and you are out of luck, the program may find a wrong redshift solution by fitting a strong emission line (e.g., [OII]) with an a priori weak line, which we know should not be the brightest line in the spectrum (e.g., [HeII] at 1640A).

In our case, this is an optical spectrum of a ```z=3``` galaxy; we know that there is little chance to see emission lines like H-alpha, because at this redshift they are moved into the near infrared domain. Instead, our spectrum is probing the rest-frame UV domain, from roughly 1200 to 2600 Angstroms. In this wavelength regime, there are only a few bright emission lines: Ly-alpha, [MgII]2799, and [OII]3727. There are also a couple of strong absorption lines (e.g., Steidel et al. 2010): SiII 1260, OI 1302, CII 1335, SiIV 1400, SiII 1526, CIV 1550. We will therefore go with this list:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1 \
    lines=[em_lyalpha,em_mg2_2799,em_o2_3727,abs_si2_1260,abs_o1_1302,abs_c2_1335,em_si4_1400,abs_si2_1526,em_c4_1550]
```

The program will then attempt to fit each of these lines with a variable flux and line width. As for the redshift, the line widths are varied on a grid with a fixed step, which is controlled by the ```delta_width``` option. By default the step in line widths is relatively small to allow for precise measurements (also a fifth of a pixel, as for the redshift). But again, for our first run this step is too small, so we will increase it a bit to one pixel (```delta_width=1```, which corresponds to a step of ```110``` km/s in the velocity dispersion for this spectrum with ```R=3000```):

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1 delta_width=1 \
    lines=[em_lyalpha,em_mg2_2799,em_o2_3727,abs_si2_1260,abs_o1_1302,abs_c2_1335,em_si4_1400,abs_si2_1526,em_c4_1550]
```

## Fixing the FITS header

Finally, we will also add the ```verbose``` option, so the code will print a summary of its progress in the terminal as the calculation goes on:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1 delta_width=1 verbose \
    lines=[em_lyalpha,em_mg2_2799,em_o2_3727,abs_si2_1260,abs_o1_1302,abs_c2_1335,em_si4_1400,abs_si2_1526,em_c4_1550]
```

If we run this command, however, we get another error message:
```bash
note: read input spectrum...
error: reading sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits
error: could not find unit of wavelength axis (missing CUNIT1 keyword)
```

It looks like the input spectrum is broken, and the program cannot read the unit of the spectral axis. This if of course crucial, because the unit could be anything from Angstroms to meters. It turns out, in this case this happens because the spectrum is missing the ```CUNIT1``` FITS keyword. This is a mistake of the person who produced this spectrum, and ideally we would ask them to fix their pipeline to always add the ```CUNIT``` keyword. But for now we need to add it ourselves. There are a multitude of ways of doing this. Here we will use a brute force method, by replacing a useless keyword:

```bash
# Warning: make sure you copy the entire command at once to run it in Bash
# and that you get to the end of the "sed ..." line (there are a lot of blank spaces
# but they are important).
cat sc_CDFS006664_P1M2Q4_P2M1Q4_003_1.fits | \
    sed "s/HIERARCH PND FILE ID = 'ce7489726d481715e83779460408ce60'/CUNIT1  = 'Angstrom'                                     /g" > \
    sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed.fits
```

This creates a new spectrum ```sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed.fits``` with the proper keyword set. We can now run the command on the fixed spectrum:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1 delta_width=1 verbose \
    lines=[em_lyalpha,em_mg2_2799,em_o2_3727,abs_si2_1260,abs_o1_1302,abs_c2_1335,em_si4_1400,abs_si2_1526,em_c4_1550]
```

## A first (not very good) run

You can now run the command, and get your first fit running! Congratulations. On my computer, this takes only one second to run. It returns a best fit redshift of ```z=2.67568```, which is surprisingly different from the redshift that was given to us in the header, but also the reduced chi squared value is horribly wrong: ```56027.7```. For a good fit, this value should be close to one. Obviously, something has gone wrong... Let's see how we can fix this.


# Improving the fit

## Visualizing the fit

In all cases, it is always advisable to visually inspect the spectrum and the best fit model together, to see if all is as we expect. By default slinefit does not write the best fit model, but we can ask it to do so by adding the ```save_model``` option:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1 delta_width=1 verbose save_model \
    lines=[em_lyalpha,em_mg2_2799,em_o2_3727,abs_si2_1260,abs_o1_1302,abs_c2_1335,em_si4_1400,abs_si2_1526,em_c4_1550]
```

This will create a file called ```sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed_slfit_model.fits```. It is a FITS file containing the best model 1D spectrum in the first extension. To visualize this, you can use any program of your choice (Python, IDL, etc.) provided you can open FITS files and make a simple plot. I will use IDL in the following, but there's nothing specific to IDL in any of this.

```idl
filebase = 'sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed'
; Read the flux spectrum
f = mrdfits(filebase+'.fits', 0, hdr, /silent)
; Read the error spectrum
e = mrdfits(filebase+'.fits', 3, /silent)
; Read the error spectrum
m = mrdfits(filebase+'_slfit_model.fits', 1, /silent)
; Create the wavelength axis
l = (dindgen(n_elements(f))+1 - sxpar(hdr, 'CRPIX1'))*sxpar(hdr, 'CDELT1') + sxpar(hdr, 'CRVAL1')

; Plot the spectrum in white and the model in red
plot, l, f
oplot, l, m, color='ff'x
```

This produces the following:

![Tutorial image 01](tutorial_01.png)

Something is obviously wrong with the data, there seems to be a strong negative spike at 10300A and the model attempts to reproduce it as best it can. This could be just noise though. Let's look at it in more detail by zooming in on the area of interest:

```idl
# Plot the flux in white, the uncertainty in yellow, the model in red
plot, l, f, xrange=[10250,10300], psym=-5
oplot, l, e, color='ffff'x
oplot, l, m, color='ff'x
```

![Tutorial image 02](tutorial_02.png)

The spike is relatively narrow, just two spectral elements, but the error spectrum is very small at that location; we would not expect such a strong negative value just from the noise. The data was probably corrupted here.

## Cutting corrupted edges

This is a relatively frequent case, where the edges of a spectrum are corrupted or extremely noisy. For this reason, slinefit allows you to discard a chosen number of spectral elements at the edge of the spectrum. This is done with the ```lambda_pad``` option. By default it is set at 5, which excludes the first and last 5 elements of the spectrum. Given the plot above, it seems we should probably increase this value, for example to 11:

```bash
/home/user/programs/slinefit/bin/slinefit sc_CDFS006664_P1M2Q4_P2M1Q4_003_1_fixed.fits \
    flux_hdu=0 error_hdu=3 z0=2.9 dz=0.6 delta_z=1 delta_width=1 verbose save_model lambda_pad=11 \
    lines=[em_lyalpha,em_mg2_2799,em_o2_3727,abs_si2_1260,abs_o1_1302,abs_c2_1335,em_si4_1400,abs_si2_1526,em_c4_1550]
```

This definitely changed things. The best fit redshift is now ```z=3.37703```, and the reduced chi squared has reduced to ```187.015```. This is much better, but still terrible!

## Fitting the continuum

Let's look again at the model and the spectrum:

```idl
# Make sure you read again the files... (not shown here for clarity)
# Plot the flux in white, the model in red
# (we truncate the Y range to the range covered by the model, to avoid the big negative spike)
plot, l, f, yrange=[min(m), max(m)]
oplot, l, m, color='ff'x
```

![Tutorial image 03](tutorial_03.png)

Ah. It seems there is some clear continuum flux in this spectrum, but currently we are only modeling lines. Let's improve this.
