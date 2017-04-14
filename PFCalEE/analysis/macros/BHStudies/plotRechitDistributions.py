#!/usr/bin/env python

from __future__ import division, print_function

import ROOT

import sys, os, math, argparse

parser=argparse.ArgumentParser(description = "Plot rechit distributions.")
parser.add_argument('--eosPrefix', action='store', help='EOS Prefix', default='root://eoscms//eos/cms/', type=str)
parser.add_argument('--inputDigiDirectory', action='store', help='Input Digi directory', default="store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/", type=str)
parser.add_argument('--hardcodedMaxNFilesCounter', action='store', help='Maximum number of files to get data from', default=0, type=int)
parser.add_argument('--inputDigiFileNamePrefix', action='store', help='Digi file name prefix', type=str, required=True)
parser.add_argument('--outputFileName', action='store', help='Output file name', type=str, required=True)
parser.add_argument('--outputDirectory', action='store', help='Path to output directory', default='rootpy_rechitDistributions', type=str)
parser.add_argument('--resetOutputDir', action='store_true', help='Whether or not to clean up output directory')
inputArguments = parser.parse_args()

# eosPrefix = "root://eoscms//eos/cms/"
eosPrefix = inputArguments.eosPrefix
inputDigiDirectory = inputArguments.inputDigiDirectory
hardcodedMaxNFilesCounter = inputArguments.hardcodedMaxNFilesCounter
inputDigiFileNamePrefix = inputArguments.inputDigiFileNamePrefix
outputFileName = inputArguments.outputFileName
outputDirectory = inputArguments.outputDirectory
resetOutputDir = inputArguments.resetOutputDir
# inputDigiDirectory = "store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/"
# hardcodedMaxNFilesCounter = 0
# inputDigiFileNamePrefix = "Digi_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_"
# # inputDigiFileNamePrefix = "DigivariableNoise_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_"
# # inputDigiFileNamePrefix = "Digi_Pu200_IC3_bhnoise0.20mips_version33_model2_Minbias_14TeV_10Events_"
# outputFileName = "recHitDistributions_zeroNoise"
# # outputFileName = "recHitDistributions_variableNoise"
# # outputFileName = "recHitDistributions_constNoise"
# outputDirectory = "rootpy_rechitDistributions"
# layerToPlot = 42
# listOfInputCellRanges = [[1450., 1650.]]
# listOfInputCellRanges = [[900., 1000.], [1370., 1550.], [2050., 2650.]]
# resetOutputDir = False
listOfInputCellRanges = [[900., 1000.], [1370., 1550.], [2050., 2650.]]
listOfCellSizeTitles = ["Small", "Intermediate", "Large"]
listOfCellSizeLegendItems = ["Small (#approx 4.5 cm^{2})","Intermediate (#approx 10 cm^{2})","Large (#approx 25 cm^{2})"]
listOfColors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kYellow]
overallTitle = "Rechit Distributions for various cell sizes"
nEnergyBins = 130
energyMin = 0.4
energyMax = 3.
startLayer = 40

# if (len(sys.argv) > 2):
#     sys.exit("Unexpected number of arguments")
# elif (len(sys.argv) == 2):
#     if (sys.argv[1] == "resetOutputDir=True"): resetOutputDir = True
#     else: sys.exit("Argument 1 not in expected format")
# else:
#     print ("Not cleaning up output directory.")
        
if (len(outputDirectory) > 0 and resetOutputDir):
    if (os.path.isdir("{outputDirectory}".format(**locals()))):
        if (os.system("rm -f {outputDirectory}/*".format(**locals())) != 0): sys.exit("Unable to clean up outputDirectory: {outputDirectory}".format(**locals()))
    else:
        if (os.system("mkdir {outputDirectory}".format(**locals())) != 0): sys.exit("Unable to create outputDirectory: {outputDirectory}".format(**locals()))
elif (len(outputDirectory) == 0):
    sys.exit("Please enter nonempty string for outputDirectory")

def addInputFilesToTree(inputTree):
    print ("Adding input files to tree...")
    consecutiveNonexistentFiles = 0
    fileCounter = 1
    proceedFlag = True
    while (proceedFlag):
        inputFileName = inputDigiDirectory + inputDigiFileNamePrefix + "%04d"%(fileCounter) + ".root"
        # print ("Checking for file..." + inputFileName)
        if (os.system("eos ls {inputFileName} &> /dev/null".format(**locals())) != 0):
            print ("Unable to find file: {inputFileName}".format(**locals()))
            consecutiveNonexistentFiles += 1
        else:
            consecutiveNonexistentFiles = 0
            inputFileName = eosPrefix + inputFileName
            inputTree.AddFile(inputFileName)
        fileCounter += 1
        proceedFlag = (consecutiveNonexistentFiles < 4)
        if (hardcodedMaxNFilesCounter > 0): proceedFlag = (proceedFlag and fileCounter <= hardcodedMaxNFilesCounter)
    print ("Finished adding input files to tree.")

def plotListOfHistograms(outputDirectory, outputFileName, histogramsStack, legend):
    outputCanvas = ROOT.TCanvas("c_rechitDistributionsComparison", "Comparison of rechit distributions", 1024, 768)
    outputCanvas.cd()
    # for histogramCounter in range(len(histogramsList)):
    # ROOT.gPad.SetLogy()
    histogramsStack.Draw("nostack")
    legend.Draw("SAME")
    outputCanvas.SaveAs("{outputDirectory}/{outputFileName}.png".format(**locals()))
    
def addRechitHistogramToStack(inputTree, nEnergyBins, energyMin, energyMax, inputCellRange, cellSizeTitle, cellSizeLegendItem, lineColor, histogramsStack, legend):
    # inputCellRange = listOfInputCellRanges[inputCellRangesCounter]
    rangeLow = inputCellRange[0]
    rangeHigh = inputCellRange[1]
    histogramName = "rechitDistribution_{cellSizeTitle}({nEnergyBins}, {energyMin}, {energyMax})".format(**locals())
    histogramToPlot = "HGCSSRecoHitVec.energy_ >> %s"%(histogramName)
    radialDistanceString = "(sqrt(HGCSSRecoHitVec.xpos_*HGCSSRecoHitVec.xpos_ + HGCSSRecoHitVec.ypos_*HGCSSRecoHitVec.ypos_))"
    # selectionCondition = "HGCSSRecoHitVec.layer_ == %d && %s > %.2f && %s < %.2f"%(layerToPlot, radialDistanceString, rangeLow, radialDistanceString, rangeHigh)
    selectionCondition = "HGCSSRecoHitVec.layer_ >= %d && %s > %.2f && %s < %.2f"%(startLayer, radialDistanceString, rangeLow, radialDistanceString, rangeHigh)
    # selectionCondition = "%s > %.2f && %s < %.2f"%(radialDistanceString, rangeLow, radialDistanceString, rangeHigh)
    print ("Histogram to retrieve: %s, Selection condition: %s"%(histogramToPlot, selectionCondition))
    inputTree.Draw(histogramToPlot, selectionCondition)
    histogramToAdd = ROOT.gDirectory.Get("rechitDistribution_{cellSizeTitle}".format(**locals()))
    # histogramToAdd.ClearUnderflowAndOverflow()
    normalizationFactor = histogramToAdd.Integral(1, nEnergyBins)
    histogramToAdd.Scale(1./normalizationFactor)
    histogramToAdd.SetLineColor(lineColor)
    legendEntry = legend.AddEntry(histogramToAdd, cellSizeLegendItem)
    legendEntry.SetTextColor(lineColor)
    histogramsStack.Add(histogramToAdd)

    # inputTree.Draw("HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_ >> xyCheck_%d"%(inputCellRangesCounter), selectionCondition)
    # xyHist = ROOT.gDirectory.Get("xyCheck_{inputCellRangesCounter}".format(**locals()))
    # tempCanvas = ROOT.TCanvas("c_xyCheck", "x,y check", 1024, 768)
    # xyHist.Draw()
    # tempCanvas.SaveAs("rootpy_rechitDistributions/xyCheck_{inputCellRangesCounter}".format(**locals()))
    
ROOT.gSystem.Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so")
histogramsStack = ROOT.THStack("histogramsStack", "Comparison of rechit distributions")
legend = ROOT.TLegend(0.55, 0.6, 0.9, 0.9)
legend.SetHeader("Colors for cell sizes:")
inputTree = ROOT.TChain("RecoTree")
addInputFilesToTree(inputTree)
totNEvents = inputTree.GetEntries()
if (totNEvents == 0):
    sys.exit("No data available!")
    print ("Number of available events = " + str(totNEvents))
for inputCellRangesCounter in range(len(listOfInputCellRanges)):
    inputCellRange = listOfInputCellRanges[inputCellRangesCounter]
    cellSizeTitle = listOfCellSizeTitles[inputCellRangesCounter]
    cellSizeLegendItem = listOfCellSizeLegendItems[inputCellRangesCounter]
    lineColor = listOfColors[inputCellRangesCounter]
    # addRechitHistogramToStack(inputTree, nEnergyBins, energyMin, energyMax, listOfInputCellRanges, listOfCellSizeNames, listOfColors, inputCellRangesCounter, histogramsStack, legend)
    addRechitHistogramToStack(inputTree, nEnergyBins, energyMin, energyMax, inputCellRange, cellSizeTitle, cellSizeLegendItem, lineColor, histogramsStack, legend)

ROOT.gStyle.SetOptStat(0)
plotListOfHistograms(outputDirectory, outputFileName, histogramsStack, legend)
