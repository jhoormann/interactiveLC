# visLC
A widget which allows you to play with the analysis choices made when generating emission line light curves from AGN spectral data to make informed decisions about which parameters to use.  After loading in your spectral data it will call the same functions used in [OzDES_makeLC](https://github.com/jhoormann/OzDES_makeLC) to make the emission light curve.  You can then interactivly play with the line integration windows, local continuum subtraction windows, and continuum uncertainty regions to determine how these choices impact the the resulting light curve and uncertainties.  

# Starting visLC
To run execute >> python visCalc.py.  You will then see the following screen where you can specify some initial parameters.

![](outputData/exampleFigs/OpenWindow.png)

## Input List File
This is the path and name of the file that lists the names of all the spectral data files you want to analyze.  You don't need to include the path to the data just the name of each file (ie Test_AGN.fits).

## Spectral Data Location
This is the path to where you are keeping your data which includes all the data files listed in the above field.

This program assumes the input spectra are saved in the fits format outputted after calibration and coadding obtained using [OzDES_calibSpec](https://github.com/jhoormann/OzDES_calibSpec).  These fits files are read in using the class defined in SpectralClass.py.  If you data is in a different form all you need to do is modify this class to handle your files.  As long as you don't change the overall class structure the rest of the code should work without issue.

## Output Data File
The chosen windows for each AGN/emission line will be saved to the location/file specified in this field.

## Max Window Shift
This is how much in each direction you can shift the windows in each direction using sliders.  You may choose a large number if you want to try out a more distant continuum window but you may want a smaller number if you want finer control over where the window goes. For the line integration windows one bound will aways be the center of the emission line.  

## Number of Bootstrap Iterations
In order to determine the uncertainty in the line flux measurements due to continnum subtraction the analysis method performs bootstrap resampling where the continuum subtraction region is randomly chosen within a larger potentially clean continuum region.  This is the number of times you want to move the continuum subtraction window to find these errors.  This is the most time consuming part of the analysis so the default is a fairly low number (15) but increasing it does not drastically effect the errors given.  

## Flux Scale Factor
One of my life mottos is that there is no point in making a plot if you can't read it.  This means I do not like axis labels that contain exponents because they tend to be unpleasent to read.  This is a factor used to scale the flux values to eliminate exponents.  If you don't care as much as I do about this you can just set it to be 1!

## Considered Lines
This allows you to select which emission lines you want to focus on (if present in your spectra).  If no boxes are check it will assume you want to look at everything.  If you want to look at another feature not specified on this list you can use 'Other' combined with a large Max Window Shift.  This will place a set of windows in the middle of the plot which you can move around to whatever feature you want to look out.  In this case the line integration window is not bounded by the location of the emission line being considered.

## Analyze/Exit
Once you have modified the fields to your specifications click 'Analyze' to start looking at the data.  If instead you have decided you don't really want to analyze data right now 'Exit' will close the program down.

## Default Parameters
The defaults for each field will allow you to analyze the example data provided here.  If your data is somewhere else and you don't want to change these everytime you run the code you can modify the defaults.  They are defined at the top of visLC.py.  Here you can also add in any other emission lines that you may be interested in.  

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
