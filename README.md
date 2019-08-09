# visLC
A widget which allows you to play with the analysis choices made when generating emission line light curves from spectral data to make informed decisions about which parameters to use.  After loading in your spectral data it will call the same functions used in [OzDES_makeLC](https://github.com/jhoormann/OzDES_makeLC) to make the emission light curve.  You can then interactivly play with the line integration windows, local continuum subtraction windows, and continuum uncertainty regions to determine how these choices impact the the resulting light curve and uncertainties.  

# Starting visLC
To run execute >> python visCalc.py.  You will then see the following screen. 

![](outputData/exampleFigs/OpenWindow.png)

This program assumes the input spectra are saved in the fits format outputted after calibration and coadding obtained using [OzDES_calibSpec](https://github.com/jhoormann/OzDES_calibSpec).  These fits files are read in using the class defined in SpectralClass.py.  If you data is in a different form all you need to do is modify this class to handle your files.  As long as you don't change the overall class structure the rest of the code should work without issue.

# Run Requirements
This code was tested using the following (as stated in requirements.txt)

python == 3.7.3

matplotlib == 3.1.0

numpy == 1.16.4

astropy == 3.2.1

scipy == 1.3.0

# Reference
If you are using this program please link to this github repository and cite the paper where this light curve extraction methodology was first presented. 

[Hoormann et al 2019, MNRAS 487, 3:3650](https://ui.adsabs.harvard.edu/abs/2019MNRAS.487.3650H/abstract)
