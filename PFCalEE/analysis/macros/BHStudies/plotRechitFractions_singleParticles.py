#!/usr/bin/env python

from __future__ import division, print_function

import ROOT

import sys, os, math, argparse

def addInputFilesToTree(inputTree, inputDigiDirectory, minEta, maxEta, deltaEta, inputDigiFileNamePrefix):
    print ("Adding input files to tree...")
    currentEtaValue = minEta
    while (currentEtaValue <= maxEta):
        inputFileName = inputDigiDirectory + inputDigiFileNamePrefix + "_eta%.3f"%(currentEtaValue) + ".root"
        if (os.system("ls -alh {inputFileName} &> /dev/null".format(**locals())) != 0):
            print ("Unable to find file: {inputFileName}".format(**locals()))
        else:
            print("Adding: {inputFileName}".format(**locals()))
            inputTree.AddFile(inputFileName)
        currentEtaValue += deltaEta
    print ("Finished adding input files to tree.")

def getMaxValueFromListOfHistograms(listOfInputHistograms):
    runningMaxValue = -1.
    for inputHistogram in listOfInputHistograms:
        currentMaxValue = inputHistogram.GetMaximum()
        if (currentMaxValue > runningMaxValue): runningMaxValue = currentMaxValue
    return runningMaxValue

def setYRangesToListMax(listOfInputHistograms):
    maximumValue = getMaxValueFromListOfHistograms(listOfInputHistograms)
    for inputHistogram in listOfInputHistograms:
        inputHistogram.GetYaxis().SetRangeUser(0., 1.1*maximumValue)

def getRatioGraph(numeratorHistogram, denominatorHistogram):
    graphXAxis = numeratorHistogram.GetXaxis()
    nXBins = graphXAxis.GetNbins()
    ratioGraph = ROOT.TGraphErrors(nXBins)
    for binNumber in range(1, 1+nXBins):
        xValue = graphXAxis.GetBinCenter(binNumber)
        xBinWidth = graphXAxis.GetBinWidth(binNumber)
        numeratorValue = numeratorHistogram.GetBinContent(binNumber)
        denominatorValue = denominatorHistogram.GetBinContent(binNumber)
        ratioValue = 0.
        ratioError = 0.
        if (denominatorValue > 0.):
            ratioValue = numeratorValue/denominatorValue
            numeratorError = numeratorHistogram.GetBinError(binNumber)
            denominatorError = denominatorHistogram.GetBinError(binNumber)
            ratioError = math.sqrt(math.pow(numeratorError,2) + math.pow(denominatorError*ratioValue,2))/denominatorValue
        ratioGraph.SetPoint(binNumber-1, xValue, ratioValue)
        ratioGraph.SetPointError(binNumber-1, xBinWidth, ratioError)
    return ratioGraph

def plotFracDistributions(outputDirectory, outputFileNamePrefix, listOfFracNames, extraInfoHistograms, etaRangeTitle, fracBelowColor, fracAboveColor, plotRatios, plotRatioErrors, fracMin, fracMax, fracStep, boundaryRechitEnergyInMips, forceNoTitle):
    for frac in listOfFracNames:
        legend = ROOT.TLegend(0.6, 0.8, 0.9, 0.9)
        legend.SetTextSize(0.03)
        legend.SetHeader("Colors for fraction distributions:")
        fracDistributionBase = extraInfoHistograms[frac+"_below"]
        fracDistributionBase.SetLineColor(fracBelowColor)
        baseLegendEntry = legend.AddEntry(fracDistributionBase, "rechit energy #leq %.1f mips"%(boundaryRechitEnergyInMips))
        baseLegendEntry.SetTextColor(fracBelowColor)
        fracDistributionToCompare = extraInfoHistograms[frac+"_above"]
        fracDistributionToCompare.SetLineColor(fracAboveColor)
        toCompareLegendEntry = legend.AddEntry(fracDistributionToCompare, "rechit energy > %.1f mips"%(boundaryRechitEnergyInMips))
        toCompareLegendEntry.SetTextColor(fracAboveColor)
        setYRangesToListMax([fracDistributionBase, fracDistributionToCompare])
        fracDistributionBase.SetTitle("Rechit fractions comparison: {frac} fraction, eta range: {etaRangeTitle}".format(**locals()))
        if (forceNoTitle): fracDistributionBase.SetTitle("")
        fracDistributionBase.GetXaxis().SetTitle("Number of {frac}s #divide total number of particles".format(**locals()))

        ratioGraph = getRatioGraph(fracDistributionBase, fracDistributionToCompare)
        ratioGraph.SetName("graph_fractionsComparison_{frac}_{etaRangeTitle}".format(**locals()))
        ratioGraph.SetTitle("rechit energy #leq %.1f mips #divide rechit energy > %.1f mips"%(boundaryRechitEnergyInMips, boundaryRechitEnergyInMips))
        LineAt1 = ROOT.TLine(-0.05, 1., 1.05, 1.)
        
        outputCanvas = ROOT.TCanvas("c_rechitFractionsComparison_{frac}_{etaRangeTitle}".format(**locals()), "Comparison of fractions", 1024, 768)
        outputCanvas.cd()
        if (plotRatios):
            outputCanvas.Divide(1,2)
            outputCanvas.cd(1)
        fracDistributionBase.Draw()
        fracDistributionToCompare.Draw("SAME")
        legend.Draw("SAME")
        if (plotRatios):
            outputCanvas.cd(2)
            if (plotRatioErrors): ratioGraph.Draw()
            else: ratioGraph.Draw("XACP*")
            ratioGraph.GetXaxis().SetTitle("Fraction of rechit energy")
            ratioGraph.GetXaxis().SetRangeUser(fracMin-0.5*fracStep, fracMax+0.5*fracStep)
            ratioGraph.GetYaxis().SetRangeUser(-0.1, 2.5)
            LineAt1.Draw()
        outputCanvas.SaveAs("{outputDirectory}/{outputFileNamePrefix}_fracComparison_{frac}_etaRange_{etaRangeTitle}.png".format(**locals()))

def plotListOfProfiles(outputDirectory, outputFileNamePrefix, listOfRechitDistributions, legend, setLogX, forceNoTitle):
    setYRangesToListMax(listOfRechitDistributions)
    outputCanvas = ROOT.TCanvas("c_rechitDistributionsComparison", "Comparison of rechit distributions", 1024, 768)
    outputCanvas.cd()
    baseProfile = listOfRechitDistributions[0]
    baseProfile.SetTitle("Rechit Energy profiles comparison")
    if (forceNoTitle): baseProfile.SetTitle("")
    baseProfile.GetXaxis().SetTitle("energy (mips)")
    listOfRechitDistributions[0].Draw()
    for profileCounter in range(1, len(listOfRechitDistributions)):
        listOfRechitDistributions[profileCounter].Draw("SAME")
    legend.Draw("SAME")
    if (setLogX): ROOT.gPad.SetLogx()
    outputCanvas.SaveAs("{outputDirectory}/{outputFileNamePrefix}_rechitDistributions.png".format(**locals()))

def addEventToHistograms(rechitvec, rechitDistributionToAdd, extraInfoHistograms, startLayer, stopLayer, printDebug, listOfFracNames, boundaryRechitEnergyInMips):
    for rechit in rechitvec:
        energy = rechit.energy()
        xpos = rechit.get_x()
        ypos = rechit.get_y()
        cellid = rechit.cellid()
        layer = rechit.layer()
        gammaFrac = rechit.get_gammaFrac()
        electronFrac = rechit.get_electronFrac()
        muonFrac = rechit.get_muonFrac()
        neutronFrac = rechit.get_neutronFrac()
        protonFrac = rechit.get_protonFrac()
        hadronFrac = rechit.get_hadronFrac()
        
        if (layer >= startLayer and layer <= stopLayer and energy > 0.):
            if (printDebug): print ("Pushing back histToFill by energy={energy} at layer={layer}".format(**locals()))
            rechitDistributionToAdd.Fill(energy)
            if (energy <= boundaryRechitEnergyInMips):
                for fracToPlot in listOfFracNames:
                    extraInfoHistograms[fracToPlot+"_below"].Fill(locals()[fracToPlot + "Frac"])
            else:
                for fracToPlot in listOfFracNames:
                    extraInfoHistograms[fracToPlot+"_above"].Fill(locals()[fracToPlot + "Frac"])

def fillHistogramsWithData(inputTree, rechitDistributionToAdd, extraInfoHistograms, startLayer, stopLayer, printDebug, listOfFracNames, boundaryRechitEnergyInMips):
    totNEvents = inputTree.GetEntries()
    if (totNEvents == 0):
        sys.exit("No data available!")
    print ("Number of available events = " + str(totNEvents))
    rechitvec = ROOT.vector("HGCSSRecoHit")()
    inputTree.SetBranchAddress("HGCSSRecoHitVec", rechitvec)
    for event in range(0, totNEvents):
        inputTree.GetEntry(event)
        addEventToHistograms(rechitvec, rechitDistributionToAdd, extraInfoHistograms, startLayer, stopLayer, printDebug, listOfFracNames, boundaryRechitEnergyInMips)

def normalizeHistogam(inputHist):
    normalizationFactor = inputHist.Integral("width")
    inputHist.Scale(1./normalizationFactor)
    
def analyzeTree(inputTree, etaRangeTitle, minHistEnergy, maxHistEnergy, nEnergyBins, fracMin, fracMax, fracStep, startLayer, stopLayer, etaRangeLegendItem, lineColor, listOfRechitDistributions, extraInfoHistograms, listOfFracNames, legend, printDebug, boundaryRechitEnergyInMips):
    rechitDistributionToAdd = ROOT.TH1F("energyProfile_{etaRangeTitle}".format(**locals()), "Energy profile", nEnergyBins, minHistEnergy, maxHistEnergy)
    fracFirstBinLowEdge = fracMin - 0.5*fracStep
    fracLastBinUpEdge = fracMax + 0.5*fracStep
    nFracBins = 0
    if (((fracLastBinUpEdge-fracFirstBinLowEdge)/fracStep) - (math.floor(0.5+(fracLastBinUpEdge-fracFirstBinLowEdge)/fracStep)) > 0.00001):
        sys.exit("Unable to form integral number of bins!")
    else:
        nFracBins = int(round(math.floor(0.5+(fracLastBinUpEdge-fracFirstBinLowEdge)/fracStep)))
    extraInfoHistograms.update({"{frac}_below".format(**locals()):ROOT.TH1F("fracHistBelow_{frac}".format(**locals()), "Fraction histogram", nFracBins, fracFirstBinLowEdge, fracLastBinUpEdge) for frac in listOfFracNames})
    extraInfoHistograms.update({"{frac}_above".format(**locals()):ROOT.TH1F("fracHistAbove_{frac}".format(**locals()), "Fraction histogram", nFracBins, fracFirstBinLowEdge, fracLastBinUpEdge) for frac in listOfFracNames})
    fillHistogramsWithData(inputTree, rechitDistributionToAdd, extraInfoHistograms, startLayer, stopLayer, printDebug, listOfFracNames, boundaryRechitEnergyInMips)
    normalizeHistogam(rechitDistributionToAdd)
    for extraInfoHistogram in extraInfoHistograms.values():
        normalizeHistogam(extraInfoHistogram)
    rechitDistributionToAdd.SetLineColor(lineColor)
    legendEntry = legend.AddEntry(rechitDistributionToAdd, etaRangeLegendItem)
    legendEntry.SetTextColor(lineColor)
    listOfRechitDistributions += [rechitDistributionToAdd]

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description = "Plot rechit distributions.")
    parser.add_argument('--inputDigiDirectory', action='store', help='Input Digi directory', default="store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/", type=str)
    parser.add_argument('--inputDigiFileNamePrefix', action='store', help='Digi file name prefix', type=str, required=True)
    parser.add_argument('--outputFileNamePrefix', action='store', help='Output file name', type=str, required=True)
    parser.add_argument('--outputDirectory', action='store', help='Path to output directory', default='rootpy_extraInfoHistograms', type=str)
    parser.add_argument('--minHistEnergy', action='store', help='minimum value of energy on x-axis', type=float, required=True)
    parser.add_argument('--maxHistEnergy', action='store', help='maximum value of energy on x-axis', type=float, required=True)
    parser.add_argument('--nEnergyBins', action='store', help='number of energy bins', default=150, type=int)
    parser.add_argument('--fracMin', action='store', help='Minimum value of fraction', type=float, default=0.)
    parser.add_argument('--fracMax', action='store', help='Maximum value of fraction', type=float, default=1.)
    parser.add_argument('--fracStep', action='store', help='Stepping size of fraction value', type=float, default=0.05)
    parser.add_argument('--startLayer', action='store', help='first start to register rechits', default=36, type=int)
    parser.add_argument('--stopLayer', action='store', help='final layer to register rechits', default=51, type=int)
    parser.add_argument('--resetOutputDir', action='store_true', help='Whether or not to clean up output directory')
    parser.add_argument('--printDebug', action='store_true', help='Whether or not to print debug info')
    parser.add_argument('--setLogX', action='store_true', help='Whether or not to set log scale along x')
    parser.add_argument('--plotRatios', action='store_true', help='Whether or not to plot ratio')
    parser.add_argument('--plotRatioErrors', action='store_true', help='Whether or not to plot error bars on ratios x')
    parser.add_argument('--boundaryRechitEnergyInMips', action='store', help='Rechit energy below and above which to plot distributions', type=float, default=0.8)
    parser.add_argument('--forceNoTitle', action='store_true', help='Whether or not to display histogram title')
    inputArguments = parser.parse_args()

    inputArgumentsDict = vars(inputArguments)
    currentModule = sys.modules[__name__]
    for (argname, argvalue) in inputArgumentsDict.items():
        print ("Option {argname} set to {argvalue}".format(**locals()))
        setattr(currentModule, argname, argvalue)
        
    deltaEta = 0.01
    ROOT.gStyle.SetOptStat(0)

    # listOfEtaRanges = [[1.5, 1.7], [1.7, 1.9], [1.9, 2.1], [1.5, 2.1]]
    # listOfEtaRangeTitles = ["outer", "intermediate", "inner", "all"]
    # listOfEtaRangeLegendItems = ["1.5 #leq #cbar#kern[0.3]{#eta}#cbar < 1.7", "1.7 #leq #cbar#kern[0.3]{#eta}#cbar < 1.9", "1.9 #leq #cbar#kern[0.3]{#eta}#cbar < 2.1", "all"]
    # listOfEtaRangeColors = [ROOT.kGreen, ROOT.kRed, ROOT.kBlue, ROOT.kBlack]
    
    listOfEtaRanges = [[1.5, 2.1]]
    listOfEtaRangeTitles = ["all"]
    listOfEtaRangeLegendItems = ["1.5 #leq #cbar#kern[0.3]{#eta}#cbar < 2.1"]
    listOfEtaRangeColors = [ROOT.kBlack]
    
    listOfFracNames = ["gamma", "electron", "muon", "proton", "neutron", "hadron"]
    fracBelowColor = ROOT.kRed
    fracAboveColor = ROOT.kBlack

    if (len(outputDirectory) > 0 and resetOutputDir):
        if (os.path.isdir("{outputDirectory}".format(**locals()))):
            if (os.system("rm -f {outputDirectory}/*".format(**locals())) != 0): sys.exit("Unable to clean up outputDirectory: {outputDirectory}".format(**locals()))
        else:
            if (os.system("mkdir {outputDirectory}".format(**locals())) != 0): sys.exit("Unable to create outputDirectory: {outputDirectory}".format(**locals()))
    elif (len(outputDirectory) == 0):
        sys.exit("Please enter nonempty string for outputDirectory")

    ROOT.gSystem.Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so")
    listOfRechitDistributions = []
    legend = ROOT.TLegend(0.55, 0.6, 0.9, 0.9)
    legend.SetHeader("Colors for eta ranges:")
    for inputEtaRangesCounter in range(len(listOfEtaRanges)):
        extraInfoHistograms = {}
        etaRange = listOfEtaRanges[inputEtaRangesCounter]
        minEta = etaRange[0]
        maxEta = etaRange[1]
        inputTree = ROOT.TChain("RecoTree")
        addInputFilesToTree(inputTree, inputDigiDirectory, minEta, maxEta, deltaEta, inputDigiFileNamePrefix)
        etaRangeTitle = listOfEtaRangeTitles[inputEtaRangesCounter]
        etaRangeLegendItem = listOfEtaRangeLegendItems[inputEtaRangesCounter]
        lineColor = listOfEtaRangeColors[inputEtaRangesCounter]
        analyzeTree(inputTree, etaRangeTitle, minHistEnergy, maxHistEnergy, nEnergyBins, fracMin, fracMax, fracStep, startLayer, stopLayer, etaRangeLegendItem, lineColor, listOfRechitDistributions, extraInfoHistograms, listOfFracNames, legend, printDebug, boundaryRechitEnergyInMips)
        plotFracDistributions(outputDirectory, outputFileNamePrefix, listOfFracNames, extraInfoHistograms, etaRangeTitle, fracBelowColor, fracAboveColor, plotRatios, plotRatioErrors, fracMin, fracMax, fracStep, boundaryRechitEnergyInMips, forceNoTitle)
        inputTree.Reset()
    plotListOfProfiles(outputDirectory, outputFileNamePrefix, listOfRechitDistributions, legend, setLogX, forceNoTitle)
