#!/usr/bin/env python

from __future__ import division, print_function

import ROOT

import numpy, sys, os, math

noiseFracDirectories = ["noiseContribToRechitsPassingThresholdInfo_variableNoise", "noiseContribToRechitsPassingThresholdInfo_constantNoise"]
noiseFracFilesCommonPrefix = ["noiseFrac", "noiseFrac"]
overrideDefaultRanges = [[0.00001, 10.], []]
overrideDefaultErrorsRanges = [[0.000001, 2.], []]
outputNames = ["variableNoise", "constantNoise"]
nFracsToCompare = len(noiseFracDirectories)
outputDirectory = "noiseFracComparison"
listOfThresholds = [0.5, 5.0]
startLayer = 36
stopLayer = 51
rMin = 900.
rMax = 2700.
nRadialBins = 16
totNLayers = 52
drawLogScale = True
paintTextFormat = ".0e"

if (len(outputDirectory) > 0):
    if (os.path.isdir("{outputDirectory}".format(**locals()))):
        if (os.system("rm -f {outputDirectory}/*".format(**locals())) != 0): sys.exit("Unable to clean up outputDirectory: {outputDirectory}".format(**locals()))
    else:
        if (os.system("mkdir {outputDirectory}".format(**locals())) != 0): sys.exit("Unable to create outputDirectory: {outputDirectory}".format(**locals()))
else:
    sys.exit("Please enter nonempty string for outputDirectory")
    
def readArraysFromFile(inputFileName, dataType):
    return numpy.loadtxt(inputFileName, dtype={'names': ('radialBinNumber', 'data'), 'formats': (numpy.intc, dataType)}, usecols=(0,1), unpack=True)

def save2DMapToFile(inputMap, outputDirectory, outputName, useLogScale, plotText, overrideDefaultRange):
    # ROOT.gStyle.Reset()
    outputCanvas = ROOT.TCanvas(inputMap.GetName(), inputMap.GetTitle(), 1024, 768)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogy(0)
    ROOT.gPad.SetLogz(0)
    if useLogScale: ROOT.gPad.SetLogz()
    if (plotText):
        ROOT.gStyle.SetPaintTextFormat(paintTextFormat)
        inputMap.Draw("TEXTCOLZ0")
    else: inputMap.Draw("COLZ0")
    if (len(overrideDefaultRange) > 0):
        print("Overriding default range to %s"%(overrideDefaultRange))
        inputMap.GetZaxis().SetRangeUser(overrideDefaultRange[0], overrideDefaultRange[1])
    ROOT.gPad.Update()
    outputCanvas.SaveAs(outputDirectory + "/" + outputName + ".png")

def initialize2DMapTo0(map2D):
    nbins_x = map2D.GetNbinsX()
    nbins_y = map2D.GetNbinsY()
    for binxcounter in range(0, 2+nbins_x):
        for binycounter in range(0, 2+nbins_y):
            map2D.SetBinContent(binxcounter, binycounter, 0.)

for thresholdCounter in range(len(listOfThresholds)):
    threshold = listOfThresholds[thresholdCounter]
    print ("Analyzing results for threshold = %.1f"%(threshold))
    noiseFracMaps = {toCompareIndex:ROOT.TH2F("noiseContribToRechits_2DMap_threshold%.1f_%s"%(threshold,outputNames[toCompareIndex]), "Noise Contribution to Rechit Energy: %s;Layer;r(mm)"%(outputNames[toCompareIndex]), stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nFracsToCompare)}
    noiseFracErrorsMaps = {toCompareIndex:ROOT.TH2F("noiseContribToRechitsErrors_2DMap_threshold%.1f_%s"%(threshold,outputNames[toCompareIndex]), "Error: Noise Contribution to Rechit Energy: %s;Layer;r(mm)"%(outputNames[toCompareIndex]), stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nFracsToCompare)}
    for toCompareIndex in range(nFracsToCompare):
        initialize2DMapTo0(noiseFracMaps[toCompareIndex])
        initialize2DMapTo0(noiseFracErrorsMaps[toCompareIndex])

    for layer in range(0, totNLayers):
        layerBinNumber = (noiseFracMaps[0]).GetXaxis().FindFixBin(1.0*layer)
        print ("Layer: {layer}".format(**locals()))
        totalCells = []
        noiseFracs = {toCompareIndex:[] for toCompareIndex in range(nFracsToCompare)}
        noiseFracErrors = {toCompareIndex:[] for toCompareIndex in range(nFracsToCompare)}

        for toCompareIndex in range(nFracsToCompare):
            radialBinNumbers, noiseFracs[toCompareIndex] = readArraysFromFile("%s/threshold_%.1f/%sInRadialBins_layer%d.dat"%(noiseFracDirectories[toCompareIndex], threshold, noiseFracFilesCommonPrefix[toCompareIndex], layer), numpy.double)
            radialBinNumbers, noiseFracErrors[toCompareIndex] = readArraysFromFile("%s/threshold_%.1f/%sErrorsInRadialBins_layer%d.dat"%(noiseFracDirectories[toCompareIndex], threshold, noiseFracFilesCommonPrefix[toCompareIndex], layer), numpy.double)
            
        for radialBinNumber in range(1, 1+nRadialBins):
            toCompareNoise = {toCompareIndex:0. for toCompareIndex in range(nFracsToCompare)}
            toCompareNoiseError = {toCompareIndex:0. for toCompareIndex in range(nFracsToCompare)}
            for toCompareIndex in range(nFracsToCompare):
                toCompareNoise[toCompareIndex] = (noiseFracs[toCompareIndex])[radialBinNumber]
                toCompareNoiseError[toCompareIndex] = (noiseFracErrors[toCompareIndex])[radialBinNumber]

            if (layer >= startLayer and layer <= stopLayer):
                for toCompareIndex in range(nFracsToCompare):
                    noiseFracMaps[toCompareIndex].SetBinContent(layerBinNumber, radialBinNumber, toCompareNoise[toCompareIndex])
                    noiseFracErrorsMaps[toCompareIndex].SetBinContent(layerBinNumber, radialBinNumber, toCompareNoiseError[toCompareIndex])
                        
    for toCompareIndex in range(nFracsToCompare):
        save2DMapToFile(noiseFracMaps[toCompareIndex], outputDirectory, "noiseFrac2DMap_threshold%.1f_%s"%(threshold, outputNames[toCompareIndex]), True, True, overrideDefaultRanges[toCompareIndex])
        save2DMapToFile(noiseFracErrorsMaps[toCompareIndex], outputDirectory, "noiseFracErrors2DMap_threshold%.1f_%s"%(threshold, outputNames[toCompareIndex]), True, True, overrideDefaultErrorsRanges[toCompareIndex])
