# ----------------------------------------------------------- #
# ----------------------- visCalc.py ------------------------ #
# ----------- https://github.com/jhoormann/visLC ------------ #
# ----------------------------------------------------------- #
# These are the functions used in visLC.py.  Much of the      #
# light curve generation code is from                         #
# https://github.com/jhoormann/OzDES_makeLC plus a few extra  #
# function to make the calculations for this specific         #
# application easier.  Code for the interactive widgets are   #
# defined in the main code.  Unless otherwise noted this code #
# was written by Janie Hoormann.                              #
# ----------------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from matplotlib.widgets import RadioButtons


# -------------------------------------------------- #
# ------------------- magToFlux -------------------- #
# -------------------------------------------------- #
# Reads in magnitude, error, and pivot wavelength    #
# and converts to f_lambda in units of ergs/s/cm^2/A #
# -------------------------------------------------- #
def magToFlux(mag, err, pivot):
    flux = (3*pow(10,18)/pow(pivot,2))*pow(10, -(2./5.)*(mag + 48.6))
    flux_err = abs(flux*(-2./5.)*2.30259*err)
    return flux, flux_err


# -------------------------------------------------- #
# -------------------- findBin ----------------------#
# -------------------------------------------------- #
# Finds the bin of the given vector (wavelength)     #
# where the specified quantity (line) is located.    #
# -------------------------------------------------- #
def findBin(line, wavelength):
    bin = 0
    for i in range(len(wavelength)-1):
        if line >= wavelength[i] and line <= wavelength[i+1]:
            bin = i
            i = len(wavelength)
        if line > wavelength[-1]:
            bin = len(wavelength)-1
            i = len(wavelength)
    return bin


# -------------------------------------------------- #
# ---------------- interpolateVals ------------------#
# -------------------------------------------------- #
# Interpolates a linear line between two points and  #
# propagates the uncertainty.                        #
# -------------------------------------------------- #
def interpolateVals(x, y, s, val):
    # uncertainty is variance

    interp = y[0] + (val - x[0]) * (y[1] - y[0]) / (x[1] - x[0])

    interp_var = s[0] + (s[0] + s[1]) * ((val - x[0]) / (x[1] - x[0])) ** 2.

    return interp, interp_var


# -------------------------------------------------- #
# ------------------ meanUncert ---------------------#
# -------------------------------------------------- #
# Finds the uncertainty corresponding to the mean    #
# of a set of numbers.                               #
# -------------------------------------------------- #
def meanUncert(variance):
    length = len(variance)
    var = 0
    num = 0
    for i in range(length):
        if np.isnan(variance[i]) == False:
            var = var + variance[i]
            num += 1

    sigma2 = (var / (num ** 2))

    return sigma2


# -------------------------------------------------- #
# ---------------- cont_fit_reject ------------------#
# -------------------------------------------------- #
# Interpolate a linear line through the mean of the  #
# continuum subtraction windows to represent the     #
# continuum and subtract this line.  Modifies the    #
# given flux/variance vectors.                       #
# -------------------------------------------------- #
def cont_fit_reject(wavelength, fluxes, variances, minWin, maxWin):

    # Define the wavelength range for the continuum model, between the mean of both windows
    wave = np.array([np.nanmean(minWin), np.nanmean(maxWin)])
    nBins = len(wavelength)

    # Determine how many epochs there are to continuum subtract
    number = int(fluxes.size / nBins)

    contArray = np.zeros([nBins, number])

    for epoch in range(number):
        if number == 1:
            flux = fluxes
            variance = variances
        else:
            flux = fluxes[:, epoch]
            variance = variances[:, epoch]

        # Calculate the average flux at each extreme of the wave vector (ie mean of the continuum subtraction window)
        fvals = np.array([np.nanmean(flux[findBin(minWin[0], wavelength):findBin(minWin[1], wavelength)]),
                          np.nanmean(flux[findBin(maxWin[0], wavelength):findBin(maxWin[1], wavelength)])])

        # Calculate the average uncertainty at each extreme of the wave vector
        svals = np.array([meanUncert(variance[findBin(minWin[0], wavelength):findBin(minWin[1], wavelength)]),
                          meanUncert(variance[findBin(maxWin[0], wavelength):findBin(maxWin[1], wavelength)])])

        cont = np.zeros(nBins)
        contVar = np.zeros(nBins)

        # Find the interpolated linear continuum model
        for i in range(nBins):
            cont[i], contVar[i] = interpolateVals(wave, fvals, svals, wavelength[i])

        # Subtract the continuum from the flux and add the error of the model in quadrature with the spectral error
        flux -= cont
        variance += contVar

        contArray[:, epoch] = cont

    return contArray


# -------------------------------------------------- #
# ----------------- computeABmag ------------------- #
# -------------------------------------------------- #
# computes the AB magnitude for given transmission   #
# functions and spectrum (f_lambda).  Returns the    #
# magnitude and variance.                            #
# -------------------------------------------------- #
def computeABmag(trans_flux, trans_wave, tmp_wave, tmp_flux, tmp_var):
    # Takes and returns variance
    # trans_ : transmission function data
    # tmp_ : spectral data

    # trans/tmp not necessarily defined over the same wavelength range
    # first determine the wavelength range over which both are defined
    minV = min(trans_wave)
    if minV < min(tmp_wave):
        minV = min(tmp_wave)
    maxV = max(trans_wave)
    if maxV > max(trans_wave):
        maxV = max(trans_wave)

    interp_wave = []
    tmp_flux2 = []
    tmp_var2 = []

    # Make new vectors for the flux just using that range (assuming spectral binning)

    for i in range(len(tmp_wave)):
        if minV < tmp_wave[i] < maxV:
            interp_wave.append(tmp_wave[i])
            tmp_flux2.append(tmp_flux[i])
            tmp_var2.append(tmp_var[i])

    # interpolate the transmission function onto this range
    # the transmission function is interpolated as it is generally much smoother than the spectral data
    trans_flux2 = interp1d(trans_wave, trans_flux)(interp_wave)

    # And now calculate the magnitude and uncertainty

    c = 2.992792e18  # Angstrom/s
    Num = np.nansum(tmp_flux2 * trans_flux2 * interp_wave)
    Num_var = np.nansum(tmp_var2 * (trans_flux2 * interp_wave) ** 2)
    Den = np.nansum(trans_flux2 / interp_wave)

    with np.errstate(divide='raise'):
        try:
            magAB = -2.5 * np.log10(Num / Den / c) - 48.60
            magABvar = 1.17882 * Num_var / (Num ** 2)
        except FloatingPointError:
            magAB = 99.
            magABvar = 99.

    return magAB, magABvar


# --------------------------------------------------- #
# --------------- uncertainty_cont ------------------ #
# --------------------------------------------------- #
# This function finds the uncertainty in line flux    #
# and width measurements.  For line flux you can      #
# input a range of potential continuum windows and    #
# it will randomly pick regions to use for continuum  #
# subtraction. You can also input a region over which #
#  to randomly choose the integration window.  These  #
# all also include flux randomization in order to     #
# consider the effect of the variance spectrum.       #
# You can also look at the effect flux randomization  #
# has on the line width measurements FWHM and         #
# velocity dispersion.  You can also specify to look  #
# at the RMS spectrum (flag='rms') for the line width #
# measurements, the default is to look at the provided#
# spectrum as is.  The error is calculated through    #
# bootstrap resampling using strapNum iterations.     #
# The standard deviation of the calculated quantity   #
# is then the associated error.                       #
# --------------------------------------------------- #
def uncertainty_cont(wavelength, flux, variance, strapNum, z, line, pivotLC, winLimMin, winLimMax, winsizeMin,
                     winsizeMax, scale, calc='cont', flag='mean', res=0):

    # calc = cont -> continuum subtraction
    # calc = win -> integration window
    # calc = fwhm -> FWHM line width: can specify flag=rms
    # calc = sigma -> line velocity dispersion: can specify flag=rms

    # Performs bootstrap resampling in the range of potentially clean continuum to determine
    # uncertainties on the flux measurement

    # Continuum window in Angstroms - will be scaled according to redshift

    # Winsize means the continuum subtraction windows are all the same size, just the locations shift
    winsizeMin = winsizeMin/(1+z)
    winsizeMax = winsizeMax/(1+z)

    lineMin = line[0]
    lineMax = line[1]

    # Option for potentially clean continuum region pass in bootstrap

    # Calculate the width of the bootstrapping region on each side of the line
    lowW = (winLimMin[1]-winLimMin[0])/(1+z)
    highW = (winLimMax[1]-winLimMax[0])/(1+z)

    # Check edge conditions: if the bootstraping region goes outside the region of the spectrograph use the spectrograph
    # bounds as the edges
    if winLimMin[0] < wavelength[0]:
        winLimMin[0] = wavelength[0]
        winLimMin[1] = (winLimMin[0] / (1 + z) + lowW) * (1 + z)
    if winLimMin[1] > wavelength[line[0]]:
        winLimMin[1] = wavelength[line[0]]
    if winLimMax[1] > wavelength[-1]:
        winLimMax[1] = wavelength[-1]
        winLimMax[0] = (winLimMax[1] / (1 + z) - highW) * (1 + z)
    if winLimMax[0] < wavelength[line[1]]:
        winLimMax[0] = wavelength[line[1]]

    # Wavelengths to choose in each window in steps of 0.5A
    winMinVect = np.arange(winLimMin[0], winLimMin[1] - (winsizeMin - 0.5) * (1 + z), 0.5 * (1 + z))
    winMaxVect = np.arange(winLimMax[0], winLimMax[1] - (winsizeMax - 0.5) * (1 + z), 0.5 * (1 + z))

    # Array of random continuum window starting points
    randVectMin = len(winMinVect) * np.random.rand(strapNum)
    randVectMin = randVectMin.astype(int)

    randVectMax = len(winMaxVect) * np.random.rand(strapNum)
    randVectMax = randVectMax.astype(int)

    # An array of values obtained through bootstrapping to determine uncertainties
    vals = np.zeros(strapNum)

    for i in range(strapNum):

        if calc == 'win':
            # subtracts from standard continuum but changes integration window, in this case feed in potential
            # integration windows instead of bootstrapping regions

            lineMinNew = findBin(winMinVect[randVectMin[i]], wavelength)
            lineMaxNew = findBin(winMaxVect[randVectMax[i]], wavelength)

            # Performs flux resampling to account for variance spectrum.  Flux values shifted by Gaussian scaled by
            # variance
            varC = np.copy(variance)
            fluxC = flux + np.random.normal(size=flux.shape) * (variance ** 0.5)

            # Continuum Subtract this new vector
            cont_fit_reject(wavelength, fluxC, varC, winLimMin, winLimMax)

            # Calculate the flux
            lc_mag, lc_mag_err = computeABmag(np.ones(len(wavelength[lineMinNew:lineMaxNew])),
                                              wavelength[lineMinNew:lineMaxNew], wavelength[lineMinNew:lineMaxNew],
                                              fluxC[lineMinNew:lineMaxNew]*scale, varC[lineMinNew:lineMaxNew]*
                                              pow(scale,2))

            vals[i], lc_mag_err = magToFlux(lc_mag, lc_mag_err**0.5, pivotLC)

        if calc == "cont":
            # changes cont region
            minWin = [winMinVect[randVectMin[i]], winMinVect[randVectMin[i]] + winsizeMin * (1 + z)]
            maxWin = [winMaxVect[randVectMax[i]], winMaxVect[randVectMax[i]] + winsizeMax * (1 + z)]

            # Performs flux resampling to account for variance spectrum.  Flux values shifted by Gaussian scaled by
            # variance
            varC = np.copy(variance)
            fluxC = flux + np.random.normal(size=flux.shape) * (variance ** 0.5)

            # Continuum Subtract this new vector
            cont_fit_reject(wavelength, fluxC, varC, minWin, maxWin)

            # Calculate the flux
            lc_mag, lc_mag_err = computeABmag(np.ones(len(wavelength[lineMin:lineMax])),wavelength[lineMin:lineMax],
                                              wavelength[lineMin:lineMax], fluxC[lineMin:lineMax]*scale,
                                              varC[lineMin:lineMax]*pow(scale, 2))

            vals[i], lc_mag_err = magToFlux(lc_mag, lc_mag_err**0.5, pivotLC)

        if calc == "fwhm":
            # Determine uncertainty in FWHM line measurement
            # do flux randomization and continuum subtraction
            varC = np.copy(variance)
            fluxC = flux + np.random.normal(size=flux.shape) * (variance ** 0.5)
            cont_fit_reject(wavelength, fluxC, varC, winLimMin, winLimMax)

            if flag == 'rms':
                # first calculate the RMS spectrum if requested
                fluxC, varC = rmsSpec(fluxC, varC)

            vals[i] = fwhm(wavelength[lineMin:lineMax], fluxC[lineMin:lineMax], res)

        if calc == "sigma":
            # Determine uncertainty in velocity dispersion measurement
            # do flux randomization and continuum subtraction
            varC = np.copy(variance)
            fluxC = flux + np.random.normal(size=flux.shape) * (variance ** 0.5)
            cont_fit_reject(wavelength, fluxC, varC, winLimMin, winLimMax)

            if flag == 'rms':
                # first calculate the RMS spectrum if requested
                fluxC, varC = rmsSpec(fluxC, varC)
            vals[i] = lineVar(wavelength[lineMin:lineMax], fluxC[lineMin:lineMax], res)

    stddev_bs = np.nanstd(vals)
    return stddev_bs

# --------------------------------------------------- #
# These are the rest of the functions written         #
# specifically for visLC.py   There are a bunch of    #
# others dealing with widget functionality that are   #
# defined in the main code.                           #
# --------------------------------------------------- #

# --------------------------------------------------- #
# ------------------- redoLCAxis -------------------- #
# --------------------------------------------------- #
# I have to clear and replot the light curve axis a   #
# lot.  This will reformat everything correctly.      #
# --------------------------------------------------- #

title_font = {'size': '15', 'color': 'black', 'weight': 'normal', 'verticalalignment': 'bottom'}
axis_font = {'size': '15'}


def redoLCAxis(ax, dates, scale):
    ax.clear()
    for label in (ax.get_xticklabels() + ax.get_yticklabels()):
        label.set_fontsize(15)
    ax.set_ylabel(r"Line Strength (" + '{:0.0e}'.format(scale/10) + " erg/s/cm$^2$/$\AA$)", **axis_font)
    ax.set_xlabel(r"MJD", **axis_font)
    ax.set_xlim([dates[0] - 100, dates[-1] + 100])
    return


# --------------------------------------------------- #
# ----------------- resize_buttons ------------------ #
# --------------------------------------------------- #
# Make radio buttons a good size.                     #
# --------------------------------------------------- #
def resize_buttons(r, f):
    "Resize all radio buttons in `r` collection by fractions `f`"
    [c.set_radius(c.get_radius()*f) for c in r.circles]


# --------------------------------------------------- #
# -------------------- findLines -------------------- #
# --------------------------------------------------- #
# Creates a dictionary that lets you know if you will #
# be able to find a given emission line in your       #
# spectrum.                                           #
# --------------------------------------------------- #

def findLines(wavelength, lineName, contWinBSMin, contWinBSMax):
    # decide which emission lines are available in the spectrum
    availLines = np.zeros(len(lineName)).astype(bool)
    availLines = dict(zip(lineName, availLines))

    for l in range(len(lineName)):
        # for a line to be in the spectrum you need to include the continuum subtraction windows as well.  This can
        # be limiting but as we need continuum subtracted spectra it is necessary.
        minWave = min(contWinBSMin[lineName[l]])
        maxWave = max(contWinBSMax[lineName[l]])

        if minWave > wavelength[0] and maxWave < wavelength[-1]:
            availLines[lineName[l]] = True
        elif lineName[l] == 'Other':
            availLines[lineName[l]] = True

    return availLines

