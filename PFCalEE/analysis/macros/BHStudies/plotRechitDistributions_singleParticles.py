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
            inputTree.AddFile(inputFileName)
        currentEtaValue += deltaEta
    print ("Finished adding input files to tree.")

def plotListOfProfiles(outputDirectory, outputFileName, listOfProfiles, legend, setLogX):
    outputCanvas = ROOT.TCanvas("c_rechitDistributionsComparison", "Comparison of rechit distributions", 1024, 768)
    outputCanvas.cd()
    baseProfile = listOfProfiles[0]
    baseProfile.SetTitle("Rechit Energy profiles comparison")
    baseProfile.GetXaxis().SetTitle("energy (mips)")
    listOfProfiles[0].Draw()
    for profileCounter in range(1, len(listOfProfiles)):
        listOfProfiles[profileCounter].Draw("SAME")
    legend.Draw("SAME")
    if (setLogX): ROOT.gPad.SetLogx()
    outputCanvas.SaveAs("{outputDirectory}/{outputFileName}.png".format(**locals()))

def addEventToProfile(rechitvec, histogramToFill, startLayer, stopLayer, printDebug):
    for rechit in rechitvec:
        energy = rechit.energy()
        xpos = rechit.get_x()
        ypos = rechit.get_y()
        cellid = rechit.cellid()
        layer = rechit.layer()
        if (layer >= startLayer and layer <= stopLayer and energy > 0.):
            if (printDebug): print ("Pushing back histToFill by energy={energy} at layer={layer}".format(**locals()))
            histogramToFill.Fill(energy)

def fillHistogramWithData(inputTree, histogramToFill, startLayer, stopLayer, printDebug):
    totNEvents = inputTree.GetEntries()
    if (totNEvents == 0):
        sys.exit("No data available!")
    print ("Number of available events = " + str(totNEvents))
    rechitvec = ROOT.vector("HGCSSRecoHit")()
    inputTree.SetBranchAddress("HGCSSRecoHitVec", rechitvec)
    for event in range(0, totNEvents):
        inputTree.GetEntry(event)
        addEventToProfile(rechitvec, histogramToFill, startLayer, stopLayer, printDebug)
    
def addEnergyProfileToList(inputTree, etaRangeTitle, minHistEnergy, maxHistEnergy, nEnergyBins, startLayer, stopLayer, etaLegendItem, lineColor, listOfProfiles, legend, printDebug):
    histogramToFill = ROOT.TH1F("energyProfile_{etaRangeTitle}".format(**locals()), "Energy profile", nEnergyBins, minHistEnergy, maxHistEnergy)
    fillHistogramWithData(inputTree, histogramToFill, startLayer, stopLayer, printDebug)
    histogramToFill.Scale(1./histogramToFill.Integral("width"))
    histogramToFill.SetLineColor(lineColor)
    legendEntry = legend.AddEntry(histogramToFill, etaLegendItem)
    legendEntry.SetTextColor(lineColor)
    listOfProfiles += [histogramToFill]

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description = "Plot rechit distributions.")
    parser.add_argument('--inputDigiDirectory', action='store', help='Input Digi directory', default="store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/", type=str)
    parser.add_argument('--inputDigiFileNamePrefix', action='store', help='Digi file name prefix', type=str, required=True)
    parser.add_argument('--outputFileName', action='store', help='Output file name', type=str, required=True)
    parser.add_argument('--outputDirectory', action='store', help='Path to output directory', default='rootpy_energyProfiles', type=str)
    parser.add_argument('--minHistEnergy', action='store', help='minimum value of energy on x-axis', type=float, required=True)
    parser.add_argument('--maxHistEnergy', action='store', help='maximum value of energy on x-axis', type=float, required=True)
    parser.add_argument('--nEnergyBins', action='store', help='number of energy bins', default=150, type=int)
    parser.add_argument('--resetOutputDir', action='store_true', help='Whether or not to clean up output directory')
    parser.add_argument('--printDebug', action='store_true', help='Whether or not to print debug info')
    parser.add_argument('--setLogX', action='store_true', help='Whether or not to set log scale along x')
    inputArguments = parser.parse_args()

    inputDigiDirectory = inputArguments.inputDigiDirectory
    inputDigiFileNamePrefix = inputArguments.inputDigiFileNamePrefix
    outputFileName = inputArguments.outputFileName
    outputDirectory = inputArguments.outputDirectory
    minHistEnergy = inputArguments.minHistEnergy
    maxHistEnergy = inputArguments.maxHistEnergy
    nEnergyBins = inputArguments.nEnergyBins
    resetOutputDir = inputArguments.resetOutputDir
    printDebug = inputArguments.printDebug
    setLogX = inputArguments.setLogX
    deltaEta = 0.01
    startLayer = 36
    stopLayer = 51

    listOfEtaRanges = [[1.5, 1.7], [1.7, 1.9], [1.9, 2.1]]
    listOfEtaRangeTitles = ["outer", "intermediate", "inner"]
    listOfEtaRangeLegendItems = ["1.5 #leq #cbar#kern[0.3]{#eta}#cbar < 1.7", "1.7 #leq #cbar#kern[0.3]{#eta}#cbar < 1.9", "1.9 #leq #cbar#kern[0.3]{#eta}#cbar < 2.1"]
    listOfEtaRangeColors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue]
    
    # listOfEtaRanges = [[1.5, 2.5]]
    # listOfEtaRangeTitles = ["base"]
    # listOfEtaRangeLegendItems = ["1.5 #leq #cbar#kern[0.3]{#eta}#cbar < 2.5"]
    # listOfEtaRangeColors = [ROOT.kBlack]

    if (len(outputDirectory) > 0 and resetOutputDir):
        if (os.path.isdir("{outputDirectory}".format(**locals()))):
            if (os.system("rm -f {outputDirectory}/*".format(**locals())) != 0): sys.exit("Unable to clean up outputDirectory: {outputDirectory}".format(**locals()))
        else:
            if (os.system("mkdir {outputDirectory}".format(**locals())) != 0): sys.exit("Unable to create outputDirectory: {outputDirectory}".format(**locals()))
    elif (len(outputDirectory) == 0):
        sys.exit("Please enter nonempty string for outputDirectory")

    ROOT.gSystem.Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so")
    listOfProfiles = []
    legend = ROOT.TLegend(0.55, 0.6, 0.9, 0.9)
    legend.SetHeader("Colors for eta ranges:")
    for inputEtaRangesCounter in range(len(listOfEtaRanges)):
        etaRange = listOfEtaRanges[inputEtaRangesCounter]
        minEta = etaRange[0]
        maxEta = etaRange[1]
        inputTree = ROOT.TChain("RecoTree")
        addInputFilesToTree(inputTree, inputDigiDirectory, minEta, maxEta, deltaEta, inputDigiFileNamePrefix)
        etaRangeTitle = listOfEtaRangeTitles[inputEtaRangesCounter]
        etaRangeLegendItem = listOfEtaRangeLegendItems[inputEtaRangesCounter]
        lineColor = listOfEtaRangeColors[inputEtaRangesCounter]
        addEnergyProfileToList(inputTree, etaRangeTitle, minHistEnergy, maxHistEnergy, nEnergyBins, startLayer, stopLayer, etaRangeLegendItem, lineColor, listOfProfiles, legend, printDebug)
        inputTree.Reset()

    ROOT.gStyle.SetOptStat(0)
    plotListOfProfiles(outputDirectory, outputFileName, listOfProfiles, legend, setLogX)
