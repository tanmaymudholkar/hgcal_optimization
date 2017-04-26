#!/usr/bin/env python

from __future__ import division, print_function

import ROOT

import numpy, sys, os, math

totalCellsDirectory = "totalNumberOfCellsInRadialBinsInfo"
baseOccupiedCellsDirectory = "occupiedNumberOfCellsInRadialBinsInfo_base"
overrideDefaultBaseRange = [0.0001, 0.05]

toCompareOccupiedCellsDirectories = ["occupiedNumberOfCellsInRadialBinsInfo_constantNoise", "occupiedNumberOfCellsInRadialBinsInfo_variableNoise", "occupiedNumberOfCellsInRadialBinsInfo_noiseEstimate_constantNoise", "occupiedNumberOfCellsInRadialBinsInfo_noiseEstimate_variableNoise"]
namesOfProfilesToCompare = ["constantNoise", "variableNoise", "noiseEstimate_constantNoise", "noiseEstimate_variableNoise"]
overrideDefaultRanges = [[], [], [0.00002, 0.03], [0.00006, 0.05]]
overrideDefaultErrorsRanges = [[], [], [], []]
fileNamePrefixesToCompare = ["occupiedCells", "occupiedCells", "occupiedCells", "occupiedCells"]
findRatioToCompare = [True, True, True, True]
# toCompareOccupiedCellsDirectories = ["occupiedNumberOfCellsInRadialBinsInfo_variableNoise"]
# namesOfProfilesToCompare = ["variableNoise"]
# overrideDefaultRanges = [[0.0001, 0.05]]
# overrideDefaultErrorsRanges = [[]]
# fileNamePrefixesToCompare = ["occupiedCells"]
# findRatioToCompare = [True]

nProfilesToCompare = len(toCompareOccupiedCellsDirectories)
outputDirectory = "rootpyOutput"
listOfThresholds = [0.5]
startLayer = 36
stopLayer = 51
rMin = 900.
rMax = 2700.
nRadialBins = 16
totNLayers = 52
listOfRadialBinsToPlot = [1, 4, 7, 10, 13, 16]
colorsForRadialBins = {1: ROOT.kRed, 4: ROOT.kGreen, 7: ROOT.kBlue, 10: ROOT.kYellow, 13: ROOT.kCyan, 16: ROOT.kMagenta}
namesForRadialBins = {1: r"900.0 < #cbar#kern[0.3]{R}#cbar < 1012.5", 4: r"1237.0 < #cbar#kern[0.3]{R}#cbar < 1350.0", 7: r"1575.0 < #cbar#kern[0.3]{R}#cbar < 1687.0", 10: r"1912.0 < #cbar#kern[0.3]{R}#cbar < 2025.0", 13: r"2250.0 < #cbar#kern[0.3]{R}#cbar < 2362.5", 16: r"2587.0 < #cbar#kern[0.3]{R}#cbar < 2700.0"}
# occupancyMax = -1.
drawLogScale = True
paintTextFormat = ".2f"

if (len(outputDirectory) > 0):
    if (os.path.isdir("{outputDirectory}".format(**locals()))):
        if (os.system("rm -f {outputDirectory}/*".format(**locals())) != 0): sys.exit("Unable to clean up outputDirectory: {outputDirectory}".format(**locals()))
    else:
        if (os.system("mkdir {outputDirectory}".format(**locals())) != 0): sys.exit("Unable to create outputDirectory: {outputDirectory}".format(**locals()))
else:
    sys.exit("Please enter nonempty string for outputDirectory")
    
def readArraysFromFile(inputFileName, dataType):
    return numpy.loadtxt(inputFileName, dtype={'names': ('radialBinNumber', 'nTotalCells'), 'formats': (numpy.intc, dataType)}, usecols=(0,1), unpack=True)

def save2DMapToFile(inputMap, outputDirectory, outputName, useLogScale, plotText, overrideDefaultRange):
    # ROOT.gStyle.Reset()
    outputCanvas = ROOT.TCanvas(inputMap.GetName(), inputMap.GetTitle(), 1024, 768)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetLogx(0)
    ROOT.gPad.SetLogy(0)
    ROOT.gPad.SetLogz(0)
    if useLogScale: ROOT.gPad.SetLogz()
    inputMap.GetYaxis().SetTitleOffset(1.3)
    inputMap.GetXaxis().SetNdivisions(16)
    for layerNumber in range(startLayer, 1+stopLayer):
        labelText = ""
        offset = 0
        if layerNumber <= 39:
            labelText = "FH"
            offset = 27
        elif layerNumber <= 51:
            labelText = "BH"
            offset = 39
        labelText += "%s"%(layerNumber-offset)
        inputMap.GetXaxis().SetBinLabel(layerNumber-startLayer+1, labelText)
        # inputMap.GetXaxis().SetBinLabel(layerNumber-startLayer+1, labelText)
        
    if (plotText):
        ROOT.gStyle.SetPaintTextFormat(paintTextFormat)
        inputMap.Draw("TEXTCOLZ0")
    else: inputMap.Draw("COLZ0")
    if (len(overrideDefaultRange) > 0):
        print("Overriding default range to %s"%(overrideDefaultRange))
        inputMap.GetZaxis().SetRangeUser(overrideDefaultRange[0], overrideDefaultRange[1])
    ROOT.gPad.Update()
    outputCanvas.SaveAs(outputDirectory + "/" + outputName + ".png")

def customizedLegend(xstart, ystart, deltax, deltay, title, textsize):
    legend = ROOT.TLegend(xstart, ystart, xstart+deltax, ystart+deltay, title)
    legend.SetTextSize(textsize)
    return legend

def save1DProfilesToFiles(occupanciesInRadialBin, occupancyErrorsInRadialBin, outputDirectory, threshold, fileNamePrefix, reducedXRange):
    outputCanvas = ROOT.TCanvas("canvas_occupancyProfiles_threshold%.1f"%(threshold), "Occupancy profiles", 1024, 768)
    multiGraph = ROOT.TMultiGraph("occupancyProfiles_threshold%.1f"%(threshold), "Mean Fraction of Cells with Deposited Rechit Energy above %.1f mips;Layer"%(threshold))
    legend = customizedLegend(0.6, 0.6, 0.3, 0.3, "Colors for radial ranges:", 0.03)
    ROOT.gPad.SetLogy(0)
    if (drawLogScale):
        legend = customizedLegend(0.1, 0.1, 0.3, 0.3, "Colors for radial ranges:", 0.03)
        ROOT.gPad.SetLogy()
    for radialBinNumber in range(1, 1+nRadialBins):
        if radialBinNumber in listOfRadialBinsToPlot:
            occupancyValuePairs = occupanciesInRadialBin[radialBinNumber]
            occupancyErrorsValuePairs = occupancyErrorsInRadialBin[radialBinNumber]
            occupancyGraph = ROOT.TGraphErrors(len(occupancyValuePairs))
            occupancyGraph.SetName("occupancy_radialBin%d"%(radialBinNumber))
            for occupancyValuePairCounter in range(len(occupancyValuePairs)):
                occupancyValuePair = occupancyValuePairs[occupancyValuePairCounter]
                occupancyErrorsValuePair = occupancyErrorsValuePairs[occupancyValuePairCounter]
                occupancyGraph.SetPoint(occupancyValuePairCounter, occupancyValuePair[0], occupancyValuePair[1])
                occupancyGraph.SetPointError(occupancyValuePairCounter, occupancyErrorsValuePair[0], occupancyErrorsValuePair[1])
            legendEntry = legend.AddEntry(occupancyGraph, namesForRadialBins[radialBinNumber], "")
            legendEntry.SetTextColor(colorsForRadialBins[radialBinNumber])
            occupancyGraph.SetLineColor(colorsForRadialBins[radialBinNumber])
            if (radialBinNumber == listOfRadialBinsToPlot[0]):
                occupancyGraph.SetTitle("Mean Fraction of Cells with Deposited Rechit Energy above %.1f mips"%(threshold))
                occupancyGraph.GetXaxis().SetTitle("Layer")
            multiGraph.Add(occupancyGraph)
    multiGraph.Draw("AC*")
    if (reducedXRange):
        multiGraph.GetXaxis().SetRangeUser(startLayer - 0.5, stopLayer + 0.5)
    else:
        multiGraph.GetXaxis().SetRangeUser(-0.5, totNLayers - 0.5)
    multiGraph.GetYaxis().SetRangeUser(0.00001, 1.1)
    outputCanvas.Update()
    legend.Draw("SAME")
    outputCanvas.SaveAs(outputDirectory + "/" + fileNamePrefix + "_threshold%.1f"%(threshold) + ".png")

def initialize2DMapTo0(map2D):
    nbins_x = map2D.GetNbinsX()
    nbins_y = map2D.GetNbinsY()
    for binxcounter in range(0, 2+nbins_x):
        for binycounter in range(0, 2+nbins_y):
            map2D.SetBinContent(binxcounter, binycounter, 0.)

for thresholdCounter in range(len(listOfThresholds)):
    threshold = listOfThresholds[thresholdCounter]
    print ("Analyzing results for threshold = %.1f"%(threshold))
    # baseOccupancyMap = ROOT.TH2F("occupancy2DMap_threshold%.1f"%(threshold), "Occupancy Map;Layer;r(mm)", stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax)
    baseOccupancyMap = ROOT.TH2F("occupancy2DMap_threshold%.1f"%(threshold), ";Layer;r(mm)", stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax)
    baseOccupancyErrorsMap = ROOT.TH2F("occupancy2DErrorsMap_threshold%.1f"%(threshold), "Occupancy Errors Map;Layer;r(mm)", stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax)
    # toCompareOccupancyMaps = {toCompareIndex:ROOT.TH2F("toCompareOccupancy2DMap_threshold%.1f_%s"%(threshold,namesOfProfilesToCompare[toCompareIndex]), "Occupancy Map: %s;Layer;r(mm)"%(namesOfProfilesToCompare[toCompareIndex]), stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nProfilesToCompare)}
    toCompareOccupancyMaps = {toCompareIndex:ROOT.TH2F("toCompareOccupancy2DMap_threshold%.1f_%s"%(threshold,namesOfProfilesToCompare[toCompareIndex]), ";Layer;r(mm)", stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nProfilesToCompare)}
    toCompareOccupancyErrorsMaps = {toCompareIndex:ROOT.TH2F("toCompareOccupancyErrors2DMap_threshold%.1f_%s"%(threshold,namesOfProfilesToCompare[toCompareIndex]), "Occupancy Errors Map: %s;Layer;r(mm)"%(namesOfProfilesToCompare[toCompareIndex]), stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nProfilesToCompare)}
    toCompareOccupancyRatioMaps = {toCompareIndex:ROOT.TH2F("toCompareOccupancyRatio2DMap_threshold%.1f_%s"%(threshold,namesOfProfilesToCompare[toCompareIndex]), ";Layer;r(mm)", stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nProfilesToCompare)}
    toCompareOccupancyRatioErrorsMaps = {toCompareIndex:ROOT.TH2F("toCompareOccupancyRatioErrors2DMap_threshold%.1f_%s"%(threshold,namesOfProfilesToCompare[toCompareIndex]), "Occupancy Ratio Errors Map: %s;Layer;r(mm)"%(namesOfProfilesToCompare[toCompareIndex]), stopLayer-startLayer+1, -0.5+startLayer, 0.5+stopLayer, nRadialBins, rMin, rMax) for toCompareIndex in range(nProfilesToCompare)}
    occupanciesInRadialBin = {radialBinNumber:[] for radialBinNumber in range(1, 1+nRadialBins)}
    occupancyErrorsInRadialBin = {radialBinNumber:[] for radialBinNumber in range(1, 1+nRadialBins)}
    scintOnlyOccupanciesInRadialBin = {radialBinNumber:[] for radialBinNumber in range(1, 1+nRadialBins)}
    scintOnlyOccupancyErrorsInRadialBin = {radialBinNumber:[] for radialBinNumber in range(1, 1+nRadialBins)}
    initialize2DMapTo0(baseOccupancyMap)
    initialize2DMapTo0(baseOccupancyErrorsMap)
    for toCompareIndex in range(nProfilesToCompare):
        initialize2DMapTo0(toCompareOccupancyMaps[toCompareIndex])
        initialize2DMapTo0(toCompareOccupancyErrorsMaps[toCompareIndex])
        initialize2DMapTo0(toCompareOccupancyRatioMaps[toCompareIndex])
        initialize2DMapTo0(toCompareOccupancyRatioErrorsMaps[toCompareIndex])
    # nbins_x = baseOccupancyMap.GetNbinsX()
    # nbins_y = baseOccupancyMap.GetNbinsY()
    # for binxcounter in range(0, 2+nbins_x):
    #     for binycounter in range(0, 2+nbins_y):
    #         for thresholdCounter in range(len(listOfThresholds)):
    #             baseOccupancyMap.SetBinContent(binxcounter, binycounter, 0.)
    #             baseOccupancyErrorsMap.SetBinContent(binxcounter, binycounter, 0.)
    for layer in range(0, totNLayers):
        layerBinNumber = baseOccupancyMap.GetXaxis().FindFixBin(1.0*layer)
        print ("Layer: {layer}".format(**locals()))
        radialBinNumbers = []
        totalCells = []
        radialBinNumbers, totalCells = readArraysFromFile("{totalCellsDirectory}/totalNumberOfCellsInRadialBins_layer{layer}.dat".format(**locals()), numpy.intc)
        radialBinNumbersTest = []
        baseOccupiedCells = []
        baseOccupiedCellsErrors = []
        toCompareOccupiedCells = {toCompareIndex:[] for toCompareIndex in range(nProfilesToCompare)}
        toCompareOccupiedCellsErrors = {toCompareIndex:[] for toCompareIndex in range(nProfilesToCompare)}
        radialBinNumbersTest, baseOccupiedCells = readArraysFromFile("%s/threshold_%.1f/occupiedCellsInRadialBins_layer%d.dat"%(baseOccupiedCellsDirectory, threshold, layer), numpy.double)
        if not(numpy.array_equal(radialBinNumbers, radialBinNumbersTest)): sys.exit("Error: radial binning not the same in total number of cells directory as in base occupied cells directory")
        radialBinNumbersTest, baseOccupiedCellsErrors = readArraysFromFile("%s/threshold_%.1f/occupiedCellsErrorsInRadialBins_layer%d.dat"%(baseOccupiedCellsDirectory, threshold, layer), numpy.double)
        if not(numpy.array_equal(radialBinNumbers, radialBinNumbersTest)): sys.exit("Error: radial binning not the same in total number of cells directory as in base occupied cells directory")

        for toCompareIndex in range(nProfilesToCompare):
            radialBinNumbersTest, toCompareOccupiedCells[toCompareIndex] = readArraysFromFile("%s/threshold_%.1f/%sInRadialBins_layer%d.dat"%(toCompareOccupiedCellsDirectories[toCompareIndex], threshold, fileNamePrefixesToCompare[toCompareIndex], layer), numpy.double)
            if not(numpy.array_equal(radialBinNumbers, radialBinNumbersTest)): sys.exit("Error: radial binning not the same in total number of cells directory as in occupied cells directory to compare")
            radialBinNumbersTest, toCompareOccupiedCellsErrors[toCompareIndex] = readArraysFromFile("%s/threshold_%.1f/%sErrorsInRadialBins_layer%d.dat"%(toCompareOccupiedCellsDirectories[toCompareIndex], threshold, fileNamePrefixesToCompare[toCompareIndex], layer), numpy.double)
            if not(numpy.array_equal(radialBinNumbers, radialBinNumbersTest)): sys.exit("Error: radial binning not the same in total number of cells directory as in occupied cells directory to compare")

        for radialBinNumber in range(1, 1+nRadialBins):
            if (totalCells[radialBinNumber] > 0):
                baseOccupancy = baseOccupiedCells[radialBinNumber]/totalCells[radialBinNumber]
                baseOccupancyError = baseOccupiedCellsErrors[radialBinNumber]/totalCells[radialBinNumber]
                occupanciesInRadialBin[radialBinNumber].append([layer*1.0, baseOccupancy])
                occupancyErrorsInRadialBin[radialBinNumber].append([0., baseOccupancyError])
                toCompareOccupancy = {toCompareIndex:0. for toCompareIndex in range(nProfilesToCompare)}
                toCompareOccupancyError = {toCompareIndex:0. for toCompareIndex in range(nProfilesToCompare)}
                occupancyRatio = {toCompareIndex:0. for toCompareIndex in range(nProfilesToCompare)}
                occupancyRatioError = {toCompareIndex:0. for toCompareIndex in range(nProfilesToCompare)}
                for toCompareIndex in range(nProfilesToCompare):
                    toCompareOccupancy[toCompareIndex] = (toCompareOccupiedCells[toCompareIndex])[radialBinNumber]/totalCells[radialBinNumber]
                    toCompareOccupancyError[toCompareIndex] = (toCompareOccupiedCellsErrors[toCompareIndex])[radialBinNumber]/totalCells[radialBinNumber]
                    if (baseOccupancy > 0 and toCompareOccupancy[toCompareIndex] > 0):
                        occupancyRatio[toCompareIndex] = toCompareOccupancy[toCompareIndex]/baseOccupancy
                        occupancyRatioError[toCompareIndex] = occupancyRatio[toCompareIndex]*math.sqrt(math.pow(baseOccupancyError/baseOccupancy,2) + math.pow(toCompareOccupancyError[toCompareIndex]/toCompareOccupancy[toCompareIndex],2))

                if (layer >= startLayer and layer <= stopLayer):
                    baseOccupancyMap.SetBinContent(layerBinNumber, radialBinNumber, baseOccupancy)
                    baseOccupancyErrorsMap.SetBinContent(layerBinNumber, radialBinNumber, baseOccupancyError)
                    scintOnlyOccupanciesInRadialBin[radialBinNumber].append([layer*1.0, baseOccupancy])
                    scintOnlyOccupancyErrorsInRadialBin[radialBinNumber].append([0., baseOccupancyError])
                    for toCompareIndex in range(nProfilesToCompare):
                        toCompareOccupancyMaps[toCompareIndex].SetBinContent(layerBinNumber, radialBinNumber, toCompareOccupancy[toCompareIndex])
                        toCompareOccupancyErrorsMaps[toCompareIndex].SetBinContent(layerBinNumber, radialBinNumber, toCompareOccupancyError[toCompareIndex])
                        toCompareOccupancyRatioMaps[toCompareIndex].SetBinContent(layerBinNumber, radialBinNumber, occupancyRatio[toCompareIndex])
                        toCompareOccupancyRatioErrorsMaps[toCompareIndex].SetBinContent(layerBinNumber, radialBinNumber, occupancyRatioError[toCompareIndex])
    
    save2DMapToFile(baseOccupancyMap, outputDirectory, "occupancy2DMap_threshold%.1f"%(threshold), True, False, overrideDefaultBaseRange)
    save2DMapToFile(baseOccupancyErrorsMap, outputDirectory, "occupancyErrors2DMap_threshold%.1f"%(threshold), True, False, [])
    for toCompareIndex in range(nProfilesToCompare):
        save2DMapToFile(toCompareOccupancyMaps[toCompareIndex], outputDirectory, "occupancy2DMap_threshold%.1f_%s"%(threshold, namesOfProfilesToCompare[toCompareIndex]), True, False, overrideDefaultRanges[toCompareIndex])
        save2DMapToFile(toCompareOccupancyErrorsMaps[toCompareIndex], outputDirectory, "occupancyErrors2DMap_threshold%.1f_%s"%(threshold, namesOfProfilesToCompare[toCompareIndex]), True, False, [])
        if (findRatioToCompare[toCompareIndex]): save2DMapToFile(toCompareOccupancyRatioMaps[toCompareIndex], outputDirectory, "occupancyRatio2DMap_threshold%.1f_%s"%(threshold, namesOfProfilesToCompare[toCompareIndex]), False, True, [])
        if (findRatioToCompare[toCompareIndex]): save2DMapToFile(toCompareOccupancyRatioErrorsMaps[toCompareIndex], outputDirectory, "occupancyRatioErrors2DMap_threshold%.1f_%s"%(threshold, namesOfProfilesToCompare[toCompareIndex]), False, True, [])
    # save1DProfilesToFiles(occupanciesInRadialBin, occupancyErrorsInRadialBin, outputDirectory, threshold, "occupancyProfiles", False)
    # save1DProfilesToFiles(scintOnlyOccupanciesInRadialBin, scintOnlyOccupancyErrorsInRadialBin, outputDirectory, threshold, "occupancyProfiles_scint", True)
