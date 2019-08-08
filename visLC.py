# ---------------------------------------------------------- #
# ------------------------- visLC.py ----------------------- #
# ----------- https://github.com/jhoormann/visLC ----------- #
# ---------------------------------------------------------- #
# This is a program I wrote which will allow you to          #
# interactively play with changing line integration windows  #
# and continuum subtraction windows and see the effect they  #
# have on the light curve.                                   #
# ---------------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox, CheckButtons
import visCalc as calcfcn
import SpectralClass as readSpec


# Here are the list of all the emission lines that will be considered.  If you have any others add the relevant info
# to the dictionaries here.  If you just want to analyse a subset of the lines just update the lineName array to
# just consider those of interest to you.

lineName = np.array(['CIV', 'MgII', 'Hbeta', 'Halpha', 'Other'])

lineLoc = {'CIV': 1549, 'MgII': 2798, 'Hbeta': 4861, 'Halpha': 6563, 'Other': 0}
lineInt = {'CIV': [1470, 1595], 'MgII': [2700, 2920], 'Hbeta': [4810, 4940], 'Halpha': [6420, 6680], 'Other': [0, 0]}

contWinMin = {'CIV': [1450, 1460], 'MgII': [2190, 2210], 'Hbeta': [4760, 4790], 'Halpha': [6190, 6210], 'Other': [0, 0]}
contWinMax = {'CIV': [1780, 1790], 'MgII': [3007, 3027], 'Hbeta': [5100, 5130], 'Halpha': [6810, 6830], 'Other': [0, 0]}

contWinBSMin = {'CIV': [1435, 1480], 'MgII': [2180, 2240], 'Hbeta': [4700, 4800], 'Halpha': [6120, 6220],
                'Other': [0, 0]}
contWinBSMax = {'CIV': [1695, 1820], 'MgII': [2987, 3057], 'Hbeta': [5080, 5180], 'Halpha': [6800, 6900],
                'Other': [0, 0]}


# If you want to change the default values for the files/locations/parameters do so here
inputVals = {}
inputVals['inputFile'] = 'inputData/input.txt'
inputVals['specLoc'] = 'inputData/'
inputVals['outputFile'] = 'outputData/output.txt'
inputVals['outputFig'] = 'outputData/'
inputVals['lineWidth'] = '1000'
inputVals['strapNum'] = '15'
inputVals['scale'] = '1e-16'
inputVals['lineCon'] = np.zeros(len(lineName)).astype(bool)
inputVals['lines'] = lineName


# You shouldn't need to change anything else in the rest of the code

# ---------------------------------------------------------- #
# ---------------------------------------------------------- #
# I have just kept all of the functions used to dictate the  #
# behaviour of widget functions here.                        #
# ---------------------------------------------------------- #
# ---------------------------------------------------------- #
# ------------------------ nightfcn ------------------------ #
# ---------------------------------------------------------- #
# This will change the epoch of spectrum plotted depending   #
# on which radio button is chosen and highlight the          #
# corresponding light curve data point.                      #
# ---------------------------------------------------------- #
def nightfcn(label):
    nightdict = dict(zip(radioLabels, radioData))
    nightlab = np.linspace(-1, numEpochs-1, numEpochs+1).astype(int)
    epochFind = dict(zip(radioLabels, nightlab))
    ydata = nightdict[label]
    s.set_ydata(ydata)
    calcfcn.redoLCAxis(ax2, spectra.dates, scale)

    if sum(mutable_object['lc_err']) == 0:
        ax2.scatter(dates, mutable_object['lc'], color='black')
        if epochFind[label] > -1:
            ax2.scatter(dates[epochFind[label]], mutable_object['lc'][epochFind[label]], color='firebrick')
            mutable_object['night'] = epochFind[label]
            cc.set_ydata(mutable_object['contArray'][:, epochFind[label]])
        else:
            cc.set_ydata(mutable_object['cont'])

    else:
        ax2.errorbar(dates, mutable_object['lc'], yerr=mutable_object['lc_err'], fmt='o', ms=6, elinewidth=1.5,
                     color='black')
        if epochFind[label] > -1:
            ax2.errorbar(dates[epochFind[label]], mutable_object['lc'][epochFind[label]],
                         yerr=mutable_object['lc_err'][epochFind[label]], fmt='o', ms=6, elinewidth=1.5,
                         color='firebrick')
            mutable_object['night'] = epochFind[label]
            cc.set_ydata(mutable_object['contArray'][:, epochFind[label]])

        else:
            cc.set_ydata(mutable_object['cont'])

    plt.draw()


# ---------------------------------------------------------- #
# ------------------------ calcFlux ------------------------ #
# ---------------------------------------------------------- #
# The will calculate the flux and uncertainty (if the button #
# is clicked!) for the new integration windows chosen via    #
# the sliders.                                               #
# ---------------------------------------------------------- #
def calcFlux(wavelength, fluxesVect, variancesVect, scale, flag, v1=[0,0]):

    if v1 == [0,0]:
        x1 = calcfcn.findBin(sLineMinWin.val, wavelength)
        x2 = calcfcn.findBin(sLineMaxWin.val, wavelength)
    else:
        x1 = v1[0]
        x2 = v1[1]

    lc = np.zeros(numEpochs)
    lc_err = np.zeros(numEpochs)
    pivotLC = pow(np.nansum(wavelength[x1:x2]) / np.nansum(1 / wavelength[x1:x2]), 0.5)

    for i in range(numEpochs):
        lc[i], lc_err[i] = calcfcn.computeABmag(np.ones(len(wavelength[x1:x2])), wavelength[x1:x2], wavelength[x1:x2],
                                               fluxesVect[x1:x2, i]*scale, variancesVect[x1:x2, i]*pow(scale,2))

        lc[i], lc_err[i] = calcfcn.magToFlux(lc[i], lc_err[i]**0.5, pivotLC)

        if flag == True:
            contwinBS1 = [sMinContBSLoc.val - sMinContBSWidth.val / 2, sMinContBSLoc.val + sMinContBSWidth.val / 2]

            contwinBS2 = [sMaxContBSLoc.val - sMaxContBSWidth.val / 2, sMaxContBSLoc.val + sMaxContBSWidth.val / 2]
            lc_err[i] = calcfcn.uncertainty_cont(wavelength, fluxesVect[:, i], variancesVect[:, i], strapNum, z,
                                                [x1, x2], pivotLC, contwinBS1, contwinBS2, sMinContWidth.val,
                                                sMaxContWidth.val, scale)
        else:
            lc_err = np.zeros(numEpochs)

    lc = lc / (scale/10)
    lc_err = lc_err / (scale/10)

    return lc, lc_err


# ---------------------------------------------------------- #
# ------------------------ redoCalc ------------------------ #
# ---------------------------------------------------------- #
# The will calculate the flux and uncertainty (if the button #
# is clicked!) for the new windows chosen via the sliders    #
# after redoing continuum subtraction for the new windows.   #
# ---------------------------------------------------------- #
def redoCalc(wavelength, fluxesVect, variancesVect, fluxCoadd, varCoadd, scale, flag):

    contwin1 = [sMinContLoc.val - sMinContWidth.val / 2, sMinContLoc.val + sMinContWidth.val / 2]

    contwin2 = [sMaxContLoc.val - sMaxContWidth.val / 2, sMaxContLoc.val + sMaxContWidth.val / 2]

    fluxes_cont = np.copy(fluxesVect)
    variances_cont = np.copy(variancesVect)

    contArray = calcfcn.cont_fit_reject(wavelength, fluxes_cont, variances_cont, contwin1, contwin2)

    fluxCoadd_cont = np.copy(fluxCoadd)
    varCoadd_cont = np.copy(varCoadd)

    cont = calcfcn.cont_fit_reject(wavelength, fluxCoadd_cont, varCoadd_cont, contwin1, contwin2)

    lc, lc_err = calcFlux(wavelength, fluxes_cont, variances_cont, scale, flag)

    return lc, lc_err, cont, contArray


# ---------------------------------------------------------- #
# ----------------------- update/C/BS ---------------------- #
# ---------------------------------------------------------- #
# From Slider data updates windows on spectrum plot and      #
# updates the light curve given these new values.            #
# ---------------------------------------------------------- #
def update(val):
    minline.set_xdata([sLineMinWin.val, sLineMinWin.val])
    maxline.set_xdata([sLineMaxWin.val, sLineMaxWin.val])

    mutable_object['lc'], mutable_object['lc_err'], mutable_object['cont'], mutable_object['contArray'] = \
        redoCalc(wavelength, fluxes, variances, fluxCoadd, varCoadd, scale, False)

    # I can't use update y_data because it doesn't support scatter plots, now I just need to redefine the axis
    calcfcn.redoLCAxis(ax2, spectra.dates, scale)

    ax2.scatter(dates, mutable_object['lc'], color='black')
    if mutable_object['night'] > -1:
        ax2.scatter(dates[mutable_object['night']], mutable_object['lc'][mutable_object['night']], color='firebrick')
        cc.set_ydata(mutable_object['contArray'][:, mutable_object['night']])

    else:
        cc.set_ydata(mutable_object['cont'])

    buttoni.color = axcolor
    fig.canvas.draw_idle()

    return


def updateC(val):
    cMinWidth = sMinContWidth.val
    cMinLoc = sMinContLoc.val

    minwincont1 = cMinLoc - cMinWidth/2
    maxwincont1 = cMinLoc + cMinWidth/2

    minCont1.set_xdata([minwincont1, minwincont1])
    maxCont1.set_xdata([maxwincont1, maxwincont1])

    cMaxWidth = sMaxContWidth.val
    cMaxLoc = sMaxContLoc.val

    minwincont2 = cMaxLoc - cMaxWidth/2
    maxwincont2 = cMaxLoc + cMaxWidth/2

    minCont2.set_xdata([minwincont2, minwincont2])
    maxCont2.set_xdata([maxwincont2, maxwincont2])

    mutable_object['lc'], mutable_object['lc_err'], mutable_object['cont'], mutable_object['contArray'] = \
        redoCalc(wavelength, fluxes, variances, fluxCoadd, varCoadd, scale, False)

    # I can't use update y_data because it doesn't support scatter plots, now I just need to redefine the axis
    calcfcn.redoLCAxis(ax2, spectra.dates, scale)

    ax2.scatter(dates, mutable_object['lc'], color='black')
    if mutable_object['night'] > -1:
        ax2.scatter(dates[mutable_object['night']], mutable_object['lc'][mutable_object['night']], color='firebrick')
        cc.set_ydata(mutable_object['contArray'][:, mutable_object['night']])

    else:
        cc.set_ydata(mutable_object['cont'])

    buttoni.color = axcolor
    fig.canvas.draw_idle()

    return


def updateBS(val):
    cMinBSWidth = sMinContBSWidth.val
    cMinBSLoc = sMinContBSLoc.val

    minwincontBS1 = cMinBSLoc - cMinBSWidth/2
    maxwincontBS1 = cMinBSLoc + cMinBSWidth/2

    minContBS1.set_xdata([minwincontBS1, minwincontBS1])
    maxContBS1.set_xdata([maxwincontBS1, maxwincontBS1])

    cMaxBSWidth = sMaxContBSWidth.val
    cMaxBSLoc = sMaxContBSLoc.val

    minwincontBS2 = cMaxBSLoc - cMaxBSWidth/2
    maxwincontBS2 = cMaxBSLoc + cMaxBSWidth/2

    minContBS2.set_xdata([minwincontBS2, minwincontBS2])
    maxContBS2.set_xdata([maxwincontBS2, maxwincontBS2])

    buttoni.color = axcolor
    fig.canvas.draw_idle()


# ---------------------------------------------------------- #
# ------------------------- errfcn ------------------------- #
# ---------------------------------------------------------- #
# If you decide you want to know the errors this will        #
# calculate and plot them.                                   #
# ---------------------------------------------------------- #
def errfcn(event):
    mutable_object['lc'], mutable_object['lc_err'], mutable_object['cont'], mutable_object['contArray'] =\
        redoCalc(wavelength, fluxes, variances, fluxCoadd, varCoadd, scale, True)

    calcfcn.redoLCAxis(ax2, spectra.dates, scale)

    ax2.errorbar(dates, mutable_object['lc'], yerr=mutable_object['lc_err'], fmt='o', ms=6, elinewidth=1.5,
                 color='black')

    if mutable_object['night'] > -1:
        ax2.errorbar(dates[mutable_object['night']], mutable_object['lc'][mutable_object['night']],
                     yerr=mutable_object['lc_err'][mutable_object['night']], fmt='o', ms=6, elinewidth=1.5,
                     color='firebrick')
        cc.set_ydata(mutable_object['contArray'][:, mutable_object['night']])

    else:
        cc.set_ydata(mutable_object['cont'])

    buttoni.color = axcolor
    fig.canvas.draw_idle()


# ---------------------------------------------------------- #
# ------------------- reset/c1/c2/BS1/BS2 ------------------ #
# ---------------------------------------------------------- #
# If clicked this button will reset the windows to their     #
# default values.                                            #
# ---------------------------------------------------------- #
def reset(event):
    sLineMinWin.reset()
    sLineMaxWin.reset()
    buttoni.color = axcolor


def resetc1(event):
    sMinContLoc.reset()
    sMinContWidth.reset()
    buttoni.color = axcolor


def resetc2(event):
    sMaxContLoc.reset()
    sMaxContWidth.reset()
    buttoni.color = axcolor


def resetcBS1(event):
    sMinContBSLoc.reset()
    sMinContBSWidth.reset()
    buttoni.color = axcolor


def resetcBS2(event):
    sMaxContBSLoc.reset()
    sMaxContBSWidth.reset()
    buttoni.color = axcolor


# ---------------------------------------------------------- #
# ------------------------- ignore ------------------------- #
# ---------------------------------------------------------- #
# Is this data and the windows not worth saving, this        #
# controls the button where you can indicate that.           #
# ---------------------------------------------------------- #
def ignore(event):
    if mutable_object['ignore'] == 0:
        mutable_object['ignore'] = 1
        buttoni.color = 'darkgrey'
    else:
        mutable_object['ignore'] = 0
        buttoni.color = axcolor


# ---------------------------------------------------------- #
# ------------------------- flagfcn ------------------------ #
# ---------------------------------------------------------- #
# Is their something peculiar about the data you want to     #
# make note of.  This flag will remind you to take another   #
# look.                                                      #
# ---------------------------------------------------------- #
def flagfcn(event):
    if mutable_object['key'] == 0:
        mutable_object['key'] = 1
        buttonflag.color = 'darkgrey'
    else:
        mutable_object['key'] = 0
        buttonflag.color = axcolor


# ---------------------------------------------------------- #
# -------------------------- exit -------------------------- #
# ---------------------------------------------------------- #
# Exits the program.                                         #
# ---------------------------------------------------------- #
def exitfcn(event):
    if mutable_object['ignore'] == 0:
        saveData()

    output.close()
    plt.close()
    raise SystemExit

def exitfcnns(event):
    plt.close()
    raise SystemExit


# ---------------------------------------------------------- #
# ------------------------- nextfcn ------------------------ #
# ---------------------------------------------------------- #
# Saves data and moves on the next line/AGN.                 #
# ---------------------------------------------------------- #
def nextfcn(event):
    if mutable_object['ignore'] == 0:
        saveData()

    plt.close()


# ---------------------------------------------------------- #
# ------------------------- savePlt ------------------------ #
# ---------------------------------------------------------- #
# Saves the figure, it will just increment a version number  #
# for each line/AGN combo so you can save many versions of   #
# analysis on the same data.                                 #
# ---------------------------------------------------------- #
def savePlt(event):
    plt.savefig(inputVals['outputFig'] + str(AGN) + "_" + line + "_version_" + str(mutable_object['fignum']) + ".png")
    mutable_object['fignum'] = mutable_object['fignum'] + 1


# ---------------------------------------------------------- #
# ------------------------- saveData ----------------------- #
# ---------------------------------------------------------- #
# Saves the final line/continuum windows.                    #
# ---------------------------------------------------------- #
def saveData():
    cMinWidth = sMinContWidth.val
    cMinLoc = sMinContLoc.val

    minwincont1 = cMinLoc - cMinWidth / 2
    maxwincont1 = cMinLoc + cMinWidth / 2

    cMaxWidth = sMaxContWidth.val
    cMaxLoc = sMaxContLoc.val

    minwincont2 = cMaxLoc - cMaxWidth / 2
    maxwincont2 = cMaxLoc + cMaxWidth / 2

    cMinBSWidth = sMinContBSWidth.val
    cMinBSLoc = sMinContBSLoc.val

    minwincontBS1 = cMinBSLoc - cMinBSWidth / 2
    maxwincontBS1 = cMinBSLoc + cMinBSWidth / 2

    cMaxBSWidth = sMaxContBSWidth.val
    cMaxBSLoc = sMaxContBSLoc.val

    minwincontBS2 = cMaxBSLoc - cMaxBSWidth / 2
    maxwincontBS2 = cMaxBSLoc + cMaxBSWidth / 2

    if mutable_object['key'] == 1:
        flag = True
    else:
        flag = False

    output.write("%s %s %.3f %.0f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f \n" % (AGN, line, z, flag,
                                                                                           sLineMinWin.val,
                                                                                           sLineMaxWin.val, minwincont1,
                                                                                           maxwincont1, minwincont2,
                                                                                           maxwincont2, minwincontBS1,
                                                                                           maxwincontBS1, minwincontBS2,
                                                                                           maxwincontBS2))

    return


# ---------------------------------------------------------- #
# All the text boxes and check box for line for intro        #
# analysis set up.                                           #
# ---------------------------------------------------------- #
def setInput(text):
   inputVals['inputFile'] = text


def setSpecLoc(text):
   inputVals['specLoc'] = text


def setOutputLoc(text):
   inputVals['outputFig'] = text


def setOutputFile(text):
   inputVals['outputFile'] = text


def setLineWidth(text):
   inputVals['lineWidth'] = text


def setStrapNum(text):
   inputVals['strapNum'] = text


def setScale(text):
   inputVals['scale'] = text


def linefcn(event):
   indexL = np.linspace(0, len(lineName) - 1, len(lineName)).astype(int)
   linelab = dict(zip(lineName, indexL))
   if inputVals['lineCon'][linelab[event]] == True:
      inputVals['lineCon'][linelab[event]] = False
   else:
      inputVals['lineCon'][linelab[event]] = True

   inputVals['lines'] = [line for line in lineName if dict(zip(lineName, inputVals['lineCon']))[line] == True]


def gofcn(event):
   plt.close()


title_font = {'size': '15', 'color': 'black', 'weight': 'normal', 'verticalalignment': 'bottom'}
axis_font = {'size': '15'}

# ---------------------------------------------------------- #
# ---------------------------------------------------------- #
# This will allow you to input/output files/directories and  #
# specify some analysis parameters.                          #
# ---------------------------------------------------------- #
# ---------------------------------------------------------- #

# Define the interactive window

fig = plt.gcf()
fig.set_size_inches(8, 8, forward=True)
ax = fig.add_subplot(111)
plt.axis('off')

axcolor = 'whitesmoke'

lineConsidered = np.zeros(len(lineName)).astype(bool)

box_input = plt.axes([0.35, 0.77, 0.45, 0.03])
text_input = TextBox(box_input, 'Input List File', initial=inputVals['inputFile'], color=axcolor, hovercolor='darkgrey')
text_input.on_submit(setInput)

box_specLoc = plt.axes([0.35, 0.72, 0.45, 0.03])
text_specLoc = TextBox(box_specLoc, 'Spectral Data Location', initial=inputVals['specLoc'], color=axcolor,
                       hovercolor='darkgrey')
text_specLoc.on_submit(setSpecLoc)


box_outputFile = plt.axes([0.35, 0.67, 0.45, 0.03])
text_outputFile = TextBox(box_outputFile, 'Output Data File', initial=inputVals['outputFile'], color=axcolor,
                         hovercolor='darkgrey')
text_outputFile.on_submit(setOutputFile)

box_outputLoc = plt.axes([0.35, 0.62, 0.45, 0.03])
text_outputLoc = TextBox(box_outputLoc, 'Output Figure Location', initial=inputVals['outputFig'], color=axcolor,
                         hovercolor='darkgrey')
text_outputLoc.on_submit(setOutputLoc)

box_lineWidth = plt.axes([0.35, 0.57, 0.45, 0.03])
text_lineWidth = TextBox(box_lineWidth, 'Max Window Shift ($\mathrm{\AA}$)', initial=inputVals['lineWidth'],
                         color=axcolor, hovercolor='darkgrey')
text_lineWidth.on_submit(setLineWidth)

box_strapNum = plt.axes([0.35, 0.52, 0.45, 0.03])
text_strapNum = TextBox(box_strapNum, 'Number of Bootstrap Iterations', initial=inputVals['strapNum'], color=axcolor,
                        hovercolor='darkgrey')
text_strapNum.on_submit(setStrapNum)

box_scale = plt.axes([0.35, 0.47, 0.45, 0.03])
text_scale = TextBox(box_scale, 'Flux Scale Factor', initial=inputVals['scale'], color=axcolor, hovercolor='darkgrey')
text_scale.on_submit(setScale)

rax = plt.axes([0.35, (0.45 - 0.035 * (len(lineName) + 1)), 0.12, 0.035 * (len(lineName) + 1)], facecolor=axcolor)
lineButtons = CheckButtons(rax, lineName, lineConsidered)
lineButtons.on_clicked(linefcn)

plt.text(0.18, 0.42, 'Considered Lines', horizontalalignment='center', verticalalignment='center',
         transform=ax.transAxes)

plt.text(0.5, 0.95, 'visLC', horizontalalignment='center', verticalalignment='center', fontsize=20,
         transform=ax.transAxes)

goax = plt.axes([0.35, (0.45 - 0.035 * (len(lineName) + 1)) - 0.05, 0.12, 0.03])

buttongo = Button(goax, 'Analyze', color=axcolor, hovercolor='darkgrey')

buttongo.on_clicked(gofcn)

exitax = plt.axes([0.48, (0.45 - 0.035 * (len(lineName) + 1)) - 0.05, 0.12, 0.03])
buttonex = Button(exitax, 'Exit', color=axcolor, hovercolor='darkgrey')

buttonex.on_clicked(exitfcnns)

plt.show()

# Now that the info as been inputted we will now use that data

# The file with the list all of the IDs that you want to analyze.
listFile = inputVals['inputFile']
idList = np.loadtxt(listFile, dtype=str)
if idList.size == 1:
    idList = np.array([idList])

# The output file is just the name of the input file with -windows at the end.
outName = inputVals['outputFile']

# Define names of the columns of info you might want to save
output = open(outName, 'a+')
output.write('ID    Line    z    Flag    LineMin    LineMax    ContMin1    ContMax1    ContMin2    ContMax2    '
             'ContBSMin1    ContBSMax1    ContBSMin2    ContBSMax2\n')

# location of the spectroscopic data to consider.
specBase = inputVals['specLoc']

lineWidth = int(inputVals['lineWidth'])  # Angstroms, value in each direction you can move each line/cont window
strapNum = int(inputVals['strapNum'])  # Number of bootstrapping iterations used for uncertainties.
scale = float(inputVals['scale'])  # value you want to scale fluxes by to make it prettier on the axis labels
linesCon = inputVals['lines']  # Maybe you only want to look at H-beta lines, perhaps CIV here is the list to focus on

# Loop through all the AGN in your list.
for AGN in idList:
    fileName = specBase + AGN
    spectra = None
    spectra = readSpec.SpectrumCoadd(fileName)

    z = spectra.redshift
    numEpochs = spectra.numEpochs
    dates = spectra.dates
    wavelength = spectra.wavelength

    nBins = len(wavelength)

    fluxCoadd = spectra.fluxCoadd/scale
    varCoadd = spectra.varianceCoadd/pow(scale, 2)

    fluxes = spectra.flux/scale
    variances = spectra.variance/pow(scale, 2)

    # Find which emission lines are present.
    availLines = calcfcn.findLines(wavelength/(1+z), linesCon, contWinBSMin, contWinBSMax)

    # Do the analysis for each present emission line.
    for line in linesCon:
        if availLines[line] == True:

            # Define the figure
            fig, (ax, ax2) = plt.subplots(2)
            for label in (ax.get_xticklabels() + ax.get_yticklabels()):
                label.set_fontsize(14)
            for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
                label.set_fontsize(14)
            fig = plt.gcf()
            fig.set_size_inches(19, 10.5, forward=True)
            plt.subplots_adjust(left=0.45)

            # First axis is for the spectra
            ax.set_title(line + " line   in   " + str(AGN) + "   at   z = " + str(round(z, 3)), **title_font)
            ax.set_ylabel(r"Flux (" + str(scale) + " erg/s/cm$^2$/$\AA$)", **axis_font)
            ax.set_xlabel(r"Wavelength ($\AA$)", **axis_font)

            xmin = wavelength[0]
            xmax = wavelength[4999]
            ymin = min(fluxCoadd)-2
            ymax = max(fluxCoadd)+2

            ax.set_xlim([xmin, xmax])
            ax.set_ylim([ymin, ymax])
            yaxis = [ymin, ymax]

            # The second axis is for the emission line
            # For the plot updating I need to continually clear and redefine axis labels so I have just created a
            # function to do it.
            calcfcn.redoLCAxis(ax2, spectra.dates, scale)

            left, width = .05, .95
            bottom, height = .25, .9
            right = left + width
            top = bottom + height

            # Now I am going to start plotting, First just plot the possible spectra
            s, = ax.plot(wavelength, fluxCoadd, lw=2, color='grey')

            axcolor = 'whitesmoke'
            rax = plt.axes([0.04, .65, 0.09, 0.016*(numEpochs+1)], facecolor=axcolor)

            radioLabels = ["" for x in range(numEpochs+1)]
            radioLabels[0] = 'Coadd'

            for i in range(numEpochs):
                radioLabels[i+1] = "Night " + str(i)
            radio = RadioButtons(rax, radioLabels, activecolor='darkgrey')
            calcfcn.resize_buttons(radio, 1.5)

            radioData = np.zeros((numEpochs+1, nBins))
            radioData[0,:] = fluxCoadd
            for i in range(numEpochs):
                radioData[i+1,:] = fluxes[:,i]

            radio.on_clicked(nightfcn)

            # Now I am going to define where the line and continuum windows are for this case

            if line == 'Other':
                lLoc = (wavelength[0] + (wavelength[-1] - wavelength[0])/2)
                lWin = np.array([lLoc - 200, lLoc + 200])
                cMin = np.array([lLoc - 350, lLoc - 300])
                cMax = np.array([lLoc + 300, lLoc + 350])
                cBSMin = np.array([lLoc - 400, lLoc - 250])
                cBSMax = np.array([lLoc + 250, lLoc + 400])
            else:
                lLoc = lineLoc[line] * (1+z)
                lWin = np.array(lineInt[line]) * (1+z)
                cMin = np.array(contWinMin[line]) * (1+z)
                cMax = np.array(contWinMax[line]) * (1+z)
                cBSMin = np.array(contWinBSMin[line]) * (1+z)
                cBSMax = np.array(contWinBSMax[line]) * (1+z)


            x1 = calcfcn.findBin(lWin[0], wavelength)
            x2 = calcfcn.findBin(lWin[1], wavelength)

            cMinWidth = cMin[1] - cMin[0]
            cMinLoc = cMin[0] + cMinWidth / 2

            cMaxWidth = cMax[1] - cMax[0]
            cMaxLoc = cMax[0] + cMaxWidth / 2

            cMinBSWidth = cBSMin[1] - cBSMin[0]
            cMinBSLoc = cBSMin[0] + cMinBSWidth / 2

            cMaxBSWidth = cBSMax[1] - cBSMax[0]
            cMaxBSLoc = cBSMax[0] + cMaxBSWidth / 2

            # Do continuum subtraction and calculate line flux
            fluxes_cont = np.copy(fluxes)
            variances_cont = np.copy(variances)
            fluxCoadd_cont = np.copy(fluxCoadd)
            varCoadd_cont = np.copy(varCoadd)

            mutable_object = {}
            mutable_object['key'] = 0
            mutable_object['ignore'] = 0
            mutable_object['night'] = -1
            mutable_object['fignum'] = 0
            mutable_object['cont'] = np.zeros(nBins)
            mutable_object['contArray'] = np.zeros([nBins, numEpochs])
            mutable_object['lc'] = np.zeros(numEpochs)
            mutable_object['lc_err'] = np.zeros(numEpochs)

            mutable_object['cont'] = calcfcn.cont_fit_reject(wavelength, fluxCoadd_cont, varCoadd_cont, cMin, cMax)

            mutable_object['contArray'] = calcfcn.cont_fit_reject(wavelength, fluxes_cont, variances_cont, cMin, cMax)

            cc, = ax.plot(wavelength, mutable_object['cont'], color='black')

            mutable_object['lc'], mutable_object['lc_err'] = calcFlux(wavelength, fluxes_cont, variances_cont, scale,
                                                                      False, [x1, x2])

            # Lets start plotting these things
            minline, = ax.plot([lWin[0], lWin[0]], yaxis, color='firebrick')
            maxline, = ax.plot([lWin[1], lWin[1]], yaxis, color='firebrick')

            ax2.scatter(dates, mutable_object['lc'], color='black')


            ignoreax = plt.axes([0.18, 0.80, 0.05, 0.03])
            buttoni = Button(ignoreax, 'Ignore', color=axcolor, hovercolor='darkgrey')

            buttoni.on_clicked(ignore)

            flagax = plt.axes([0.24, 0.80, 0.05, 0.03])
            buttonflag = Button(flagax, 'Flag', color=axcolor, hovercolor='darkgrey')

            buttonflag.on_clicked(flagfcn)

            flagsave = plt.axes([0.24, 0.86, 0.05, 0.03])
            buttonsave = Button(flagsave, 'Save Fig', color=axcolor, hovercolor='darkgrey')

            buttonsave.on_clicked(savePlt)

            calcErr = plt.axes([0.18, 0.86, 0.05, 0.03])
            buttonErr = Button(calcErr, 'Calc Errors', color=axcolor, hovercolor='darkgrey')

            buttonErr.on_clicked(errfcn)

            nextB = plt.axes([0.18, 0.92, 0.05, 0.03])
            buttonNext = Button(nextB, 'Next', color=axcolor, hovercolor='darkgrey')

            buttonNext.on_clicked(nextfcn)

            exitB = plt.axes([0.24, 0.92, 0.05, 0.03])
            buttonExit = Button(exitB, 'Exit', color=axcolor, hovercolor='darkgrey')

            buttonExit.on_clicked(exitfcn)

            # Add in slider for emission line
            axLineMaxWin = plt.axes([0.18, 0.69, 0.18, 0.02], facecolor=axcolor)
            axLineMinWin = plt.axes([0.18, 0.72, 0.18, 0.02], facecolor=axcolor)

            minWinLoc = lWin[0]-lineWidth
            if minWinLoc < wavelength[0]:
                minWinLoc = wavelength[0]
            maxWinLoc = lWin[1]+lineWidth
            if minWinLoc > wavelength[-1]:
                minWinLoc = wavelength[-1]

            if line == 'Other':
                sLineMinWin = Slider(axLineMinWin, 'Line Min', minWinLoc, maxWinLoc, valinit=lWin[0], color='firebrick')
                sLineMaxWin = Slider(axLineMaxWin, 'Line Max', minWinLoc, maxWinLoc, valinit=lWin[1], color='firebrick')
            else:
                sLineMinWin = Slider(axLineMinWin, 'Line Min', minWinLoc, lLoc, valinit=lWin[0], color='firebrick')
                sLineMaxWin = Slider(axLineMaxWin, 'Line Max', lLoc, maxWinLoc, valinit=lWin[1], color='firebrick')
            sLineMinWin.on_changed(update)
            sLineMaxWin.on_changed(update)

            resetax = plt.axes([0.18, 0.65, 0.05, 0.03])
            button = Button(resetax, 'Reset', color=axcolor, hovercolor='darksalmon')

            button.on_clicked(reset)

            # Now I am going to add in the continuum subtraction window
            # First the left continuum subtraction window
            minCont1, = ax.plot([cMin[0], cMin[0]], yaxis, color='mediumblue', linestyle='--')
            maxCont1, = ax.plot([cMin[1], cMin[1]], yaxis, color='mediumblue', linestyle='--')

            axMinContWidth = plt.axes([0.18, 0.53, 0.18, 0.02], facecolor=axcolor)
            axMinContLoc = plt.axes([0.18, 0.56, 0.18, 0.02], facecolor=axcolor)

            minwincont1 = cMin[0]
            maxwincont1 = cMin[1]

            cMinLocMin = cMin[0] - lineWidth
            if cMinLocMin < wavelength[0]:
                cMinLocMin = wavelength[0]
            cMinLocMax = cMin[1] + lineWidth
            if cMinLocMax > wavelength[-1]:
                cMinLocMax = wavelength[-1]

            sMinContWidth = Slider(axMinContWidth, 'Left Continuum Width', 0, 4*cMinWidth, valinit=cMinWidth,
                                   color='mediumblue')
            sMinContLoc = Slider(axMinContLoc, 'Left Continuum Location', cMinLocMin, cMinLocMax, valinit=cMinLoc,
                                 color='mediumblue')
            sMinContLoc.on_changed(updateC)
            sMinContWidth.on_changed(updateC)

            resetaxc1 = plt.axes([0.18, 0.49, 0.05, 0.03])
            buttonc1 = Button(resetaxc1, 'Reset', color=axcolor, hovercolor='lightsteelblue')

            buttonc1.on_clicked(resetc1)

            # Now let's do the right continuum subtraction window

            minCont2, = ax.plot([cMax[0], cMax[0]], yaxis, color='rebeccapurple', linestyle='--')
            maxCont2, = ax.plot([cMax[1], cMax[1]], yaxis, color='rebeccapurple', linestyle='--')

            axMaxContWidth = plt.axes([0.18, 0.40, 0.18, 0.02], facecolor=axcolor)
            axMaxContLoc = plt.axes([0.18, 0.43, 0.18, 0.02], facecolor=axcolor)

            minwincont2 = cMax[0]
            maxwincont2 = cMax[1]

            cMaxLocMin = cMax[0] - lineWidth
            if cMaxLocMin < wavelength[0]:
                cMaxLocMin = wavelength[0]
            cMaxLocMax = cMax[1] + lineWidth
            if cMaxLocMax > wavelength[-1]:
                cMaxLocMax = wavelength[-1]

            sMaxContWidth = Slider(axMaxContWidth, 'Right Continuum Width', 0, 4*cMaxWidth, valinit=cMaxWidth,
                                   color='rebeccapurple')
            sMaxContLoc = Slider(axMaxContLoc, 'Right Continuum Location', cMaxLocMin, cMaxLocMax, valinit=cMaxLoc,
                                 color='rebeccapurple')

            sMaxContLoc.on_changed(updateC)
            sMaxContWidth.on_changed(updateC)

            resetaxc2 = plt.axes([0.18, 0.36, 0.05, 0.03])
            buttonc2 = Button(resetaxc2, 'Reset', color=axcolor, hovercolor='thistle')

            buttonc2.on_clicked(resetc2)

            # Now it is time to start thinking about uncertainties. I going to allow the user to move the bootstrapping
            # windows with the sliders but I am not going to update error bars on the fly - there will be a button for
            # that!

            minContBS1, = ax.plot([cBSMin[0], cBSMin[0]], yaxis, color='forestgreen', linestyle=':')
            maxContBS1, = ax.plot([cBSMin[1], cBSMin[1]], yaxis, color='forestgreen', linestyle=':')

            axMinContBSWidth = plt.axes([0.18, 0.27, 0.18, 0.02], facecolor=axcolor)
            axMinContBSLoc = plt.axes([0.18, 0.30, 0.18, 0.02], facecolor=axcolor)

            minwincontBS1 = cBSMin[0]
            maxwincontBS1 = cBSMin[1]

            cMinBSLocMin = cBSMin[0] - lineWidth
            if cMinBSLocMin < wavelength[0]:
                cMinBSLocMin = wavelength[0]
            cMinBSLocMax = cBSMin[1] + lineWidth
            if cMinBSLocMax > wavelength[-1]:
                cMinBSLocMax = wavelength[-1]

            cMinWidth = sMinContWidth.val
            cMinLoc = sMinContLoc.val

            lowMinWidth = 1.5*sMinContWidth.val

            sMinContBSWidth = Slider(axMinContBSWidth, 'Left Error Width', 1.5*sMinContWidth.val, 4*cMinBSWidth,
                                     valinit=cMinBSWidth, color='forestgreen')
            sMinContBSLoc = Slider(axMinContBSLoc, 'Left Error Location', cMinBSLocMin, cMinBSLocMax, valinit=cMinBSLoc,
                                   color='forestgreen')
            sMinContBSLoc.on_changed(updateBS)
            sMinContBSWidth.on_changed(updateBS)

            resetaxcBS1 = plt.axes([0.18, 0.23, 0.05, 0.03])
            buttoncBS1 = Button(resetaxcBS1, 'Reset', color=axcolor, hovercolor='darkseagreen')

            buttoncBS1.on_clicked(resetcBS1)

            minContBS2, = ax.plot([cBSMax[0], cBSMax[0]], yaxis, color='darkorange', linestyle=':')
            maxContBS2, = ax.plot([cBSMax[1], cBSMax[1]], yaxis, color='darkorange', linestyle=':')

            axMaxContBSWidth = plt.axes([0.18, 0.14, 0.18, 0.02], facecolor=axcolor)
            axMaxContBSLoc = plt.axes([0.18, 0.17, 0.18, 0.02], facecolor=axcolor)

            minwincontBS2 = cBSMax[0]
            maxwincontBS2 = cBSMax[1]

            cMaxBSLocMin = cBSMax[0] - lineWidth
            if cMaxBSLocMin < wavelength[0]:
                cMaxBSLocMin = wavelength[0]
            cMaxBSLocMax = cBSMax[1] + lineWidth
            if cMaxBSLocMax > wavelength[-1]:
                cMaxBSLocMax = wavelength[-1]

            sMaxContBSWidth = Slider(axMaxContBSWidth, 'Right Error Width', 1.5*cMaxWidth, 4*cMaxBSWidth,
                                     valinit=cMaxBSWidth, color='darkorange')
            sMaxContBSLoc = Slider(axMaxContBSLoc, 'Right Error Location', cMaxBSLocMin, cMaxBSLocMax,
                                   valinit=cMaxBSLoc, color='darkorange')
            sMaxContBSLoc.on_changed(updateBS)
            sMaxContBSWidth.on_changed(updateBS)

            resetaxcBS2 = plt.axes([0.18, 0.10, 0.05, 0.03])
            buttoncBS2 = Button(resetaxcBS2, 'Reset', color=axcolor, hovercolor='navajowhite')

            buttoncBS2.on_clicked(resetcBS2)

            plt.show()
