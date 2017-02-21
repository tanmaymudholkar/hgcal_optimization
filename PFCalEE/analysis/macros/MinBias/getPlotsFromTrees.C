#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include<ctime>
#include<algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "Rtypes.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

#include "HGCSSInfo.hh"
#include "HGCSSDetector.hh"
// #include "HGCSSGeometryConversion.hh"

/// std::vector<Int_t> etaBinsToPlot(1,3,5,7,9);
const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");
// const TString inputFileName = Form("/afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias/version_33/PythiaTest/occupiedCellsTreesMaxEta3.root");
const TString inputFileName = Form("/afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias/version_33/PythiaTest/occupiedCellsTrees.root");
const TString plotsDir = Form("plots/");
const TString outputFileName = Form("plots/plotsCollection.root");
const std::string inputEtaRangesFileName = "etaBinBoundaries.dat";
const TString totalNumberOfCellsInEtaBinsDir = Form("totalNumberOfCellsInEtaBinsInfo");
// const Double_t listOfThresholds[] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 5.5}; // assumed to be in ascending order of threshold in mips
// const Double_t listOfThresholds[] = {0.5, 2.0, 5.0}; // assumed to be in ascending order of threshold in mips
// const Double_t listOfThresholds[] = {0.5, 5.0, 5.5}; // assumed to be in ascending order of threshold in mips
const Double_t listOfThresholds[] = {0.5, 5.0}; // assumed to be in ascending order of threshold in mips
// const Double_t listOfThresholds[] = {5.0, 5.5}; // assumed to be in ascending order of threshold in mips
Int_t listOfEtaBinsToPlot[] = {9,7,5,3,1}; // assumed to be in descending order of eta
Int_t startEtaBinToCountInTotalOccupancy = 1;
Int_t stopEtaBinToCountInTotalOccupancy = 10;
unsigned startLayers[] = {0,28,40};
std::vector<Int_t> etaBinsToPlot(listOfEtaBinsToPlot, listOfEtaBinsToPlot+sizeof(listOfEtaBinsToPlot)/sizeof(listOfEtaBinsToPlot[0]));
std::map<Int_t, Int_t> colorsForEtaBins = {{1, kRed}, {3, kGreen}, {5, kBlue}, {7, kYellow}, {9, kBlack}}; // number of colors should be equal to size of list of eta bins to plot
std::map<Int_t, TString> namesForEtaBins = {{1, Form("1.5 < #cbar#kern[0.3]{#eta}#cbar < 1.6")}, {3, Form("1.8 < #cbar#kern[0.3]{#eta}#cbar < 2.0")}, {5, Form("2.1 < #cbar#kern[0.3]{#eta}#cbar < 2.2")}, {7, Form("2.4 < #cbar#kern[0.3]{#eta}#cbar < 2.5")}, {9, Form("2.7 < #cbar#kern[0.3]{#eta}#cbar < 2.8")}};
const TString versionName = Form("withMedian");

HGCSSInfo *info;
double calorSizeXY = 0;
double cellSize = 0;
unsigned versionNumber = 0;
unsigned model = 0;
unsigned nLayers;
unsigned totalNEvents;

HGCSSDetector & myDetector = theDetector();

// HGCSSGeometryConversion *geomConv;

Double_t *etaBinEdges;
// Double_t *layerZPositions;
unsigned nEtaBins;

TAxis *axisWithEtaBinning;

std::map<unsigned, std::map<unsigned, Int_t> > totalNumberOfCellsMap; // first index: layer number, second index: eta bin number

std::map<Int_t, TProfile*> occupiedCellsMap;
std::map<Int_t, std::map<Int_t, std::vector<Int_t> > > occupiedCellsVectorMap;
std::map<Int_t, TGraphErrors*> occupancyMap;
std::map<Int_t, TGraphErrors*> occupancyMedianMap;
std::map<Double_t, Int_t> occupiedCellsVsThresholdMap;

std::map<unsigned, Double_t> occupancyMaxval;

void initializeOccupancyMaxval() {
  for (unsigned startLayerCounter = 0; startLayerCounter < (sizeof(startLayers)/sizeof(startLayers[0])); ++startLayerCounter) {
    unsigned startLayer = startLayers[startLayerCounter];
    occupancyMaxval[startLayer] = -1.;
  }
}

void fillEtaBinEdges() {
  std::vector<Double_t> binEdges;
  std::string line;
  std::ifstream inputFile(inputEtaRangesFileName.c_str());
  if (inputFile.is_open()) {
    while(std::getline(inputFile, line)) {
      binEdges.push_back(atof(line.c_str()));
    }
    inputFile.close();
  }
  else {
    std::cout << "ERROR: Unable to open input file to read in eta ranges" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  etaBinEdges = new Double_t[binEdges.size()];
  nEtaBins = binEdges.size()-1;
  for (unsigned edgeCounter = 0; edgeCounter < binEdges.size(); ++edgeCounter) {
    etaBinEdges[edgeCounter] = binEdges[edgeCounter];
  }
  binEdges.clear();
}

void initializeAxisWithEtaBinning() {
  axisWithEtaBinning = new TAxis(nEtaBins, etaBinEdges);
}

void printAxisBinning(TAxis *inputAxis) {
  unsigned nBins = inputAxis->GetNbins();
  for (unsigned binCounter = 0; binCounter <= 1+nBins; ++binCounter) {
    std::cout << "Bin number " << binCounter << " ranges from " << inputAxis->GetBinLowEdge(binCounter) << " to " << inputAxis->GetBinUpEdge(binCounter) << " and its center is " << inputAxis->GetBinCenter(binCounter) << std::endl;
  }
}

bool testFile(TString inputPath) {
  TFile *testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    return false;
  }
  return true;
}

void initializeCalorimeterProperties() {
  std::cout << "Initializing calorimeter properties..." << std::endl;
  // TString digiFileName = dataDir+Form("DiginoXTalk_Pu200_IC3_version33_000001.root");
  TString digiFileName = dataDir+Form("DigiPu200_IC3_version33_000001.root");
  bool dataExists = testFile(digiFileName);
  if (dataExists) {
    TFile *inputFile = TFile::Open(digiFileName);
    info=(HGCSSInfo*)inputFile->Get("Info");
    calorSizeXY = info->calorSizeXY();
    cellSize = info->cellSize();
    versionNumber = info->version();
    model = info->model();
    inputFile->Close();

    std::cout << " -- Calor size XY = " << calorSizeXY
              << ", version number = " << versionNumber 
              << ", model = " << model
              << ", cellSize = " << cellSize
              << std::endl;
  }
  else {
    std::cout << "Input digi files not found!" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  myDetector.buildDetector(versionNumber,true,false,false);
  nLayers = myDetector.nLayers();
}

void readTotalNumberOfCellsFromFiles() {
  std::cout << "Reading total number of cells from files..." << std::endl;
  for (unsigned layerCounter = 0; layerCounter < nLayers; ++layerCounter) {
    std::map<unsigned, Int_t> totalNumberOfCellsInEtaBinsByLayer;
    std::ifstream inputFile;
    std::string line;
    unsigned binNumber;
    Int_t totalNumberOfCells;
    inputFile.open((totalNumberOfCellsInEtaBinsDir+Form("/totalNumberOfCellsInEtaBins_layer%i.dat", layerCounter)).Data());
    if (inputFile.is_open()) {
      while(std::getline(inputFile, line)) {
        std::istringstream issline(line);
        issline >> binNumber >> totalNumberOfCells;
        totalNumberOfCellsInEtaBinsByLayer[binNumber] = totalNumberOfCells;
      }
      inputFile.close();
    }
    else {
      std::cout << "ERROR: Unable to open input file to read in total number of cells in eta bins for layer index " << layerCounter << std::endl;
      std::exit(EXIT_FAILURE);
    }
    totalNumberOfCellsMap[layerCounter] = totalNumberOfCellsInEtaBinsByLayer;
  }
}

void initializeOccupiedCellsMaps(Double_t threshold, unsigned startLayer) {
  for (unsigned etaBinCounter = 0; etaBinCounter <= (1+nEtaBins); ++etaBinCounter) {
    TString profileName = Form("occupiedCells_etaBinCounter%i_startLayer%i_threshold%.1f", etaBinCounter, startLayer, threshold);
    occupiedCellsMap[etaBinCounter] = new TProfile(profileName, profileName, nLayers-startLayer, startLayer-0.5, nLayers-0.5);
    for (unsigned layerCounter = startLayer; layerCounter < nLayers; ++layerCounter) {
      occupiedCellsVectorMap[etaBinCounter][layerCounter].push_back(0);
      occupiedCellsVectorMap[etaBinCounter][layerCounter].clear();
    }
  }
  occupiedCellsVsThresholdMap[threshold] = 0;
}

void printETA(Double_t timeRemainingInSec) {
  Int_t hoursRemaining = static_cast<Int_t>(timeRemainingInSec/3600.);
  Int_t minutesRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining)/60.);
  Int_t secondsRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining - 60*minutesRemaining) + 0.5);
  
  std::cout << Form("\rApproximately %02d h: %02d m: %02d s remaining...", hoursRemaining, minutesRemaining, secondsRemaining) << std::flush;
}

void fillOccupiedCellsMapsFromData(TFile *inputFile, Double_t threshold, unsigned startLayer) {
  std::cout << "Filling occupied cells map from data..." << std::endl;
  TTree *inputTree = (TTree*)inputFile->Get(Form("occupiedCellsTree_threshold_%.1f", threshold));
  unsigned nEntries = inputTree->GetEntries();
  totalNEvents = nEntries/(nLayers*(nEtaBins+2));
  std::cout << "Checking: Total available number of events = " << totalNEvents << std::endl;
  Int_t etaBin = -1;
  unsigned layer = 1000;
  Int_t occupiedCells = -1;
  inputTree->SetBranchAddress("etaBin",&etaBin);
  inputTree->SetBranchAddress("layer",&layer);
  inputTree->SetBranchAddress("occupiedCells",&occupiedCells);

  Double_t fractionCompleted;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  std::time_t initialTimeClocked = std::time(nullptr);
  inputTree->GetEntry(0);
  
  if (layer >= startLayer) {
    occupiedCellsMap[etaBin]->Fill(layer, occupiedCells);
    occupiedCellsVectorMap[etaBin][layer].push_back(occupiedCells);
    if (etaBin >= startEtaBinToCountInTotalOccupancy && etaBin <= stopEtaBinToCountInTotalOccupancy) {
      occupiedCellsVsThresholdMap[threshold] += occupiedCells;
    }
  }
    
  for (unsigned entryCounter(1); entryCounter<nEntries; ++entryCounter){// loop on events
    fractionCompleted = float(entryCounter*1./nEntries);
    timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
    guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
    printETA(guess_timeRemainingInSec);

    inputTree->GetEntry(entryCounter);
    if (layer >= startLayer) {
      occupiedCellsMap[etaBin]->Fill(layer, occupiedCells);
      occupiedCellsVectorMap[etaBin][layer].push_back(occupiedCells);
      if (etaBin >= startEtaBinToCountInTotalOccupancy && etaBin <= stopEtaBinToCountInTotalOccupancy) {
        occupiedCellsVsThresholdMap[threshold] += occupiedCells;
      }
    }
  }
  std::cout << "Completed filling occupied cells map from data!" << std::endl;
}

void plotVerticalLayerBoundaries(TCanvas *outputCanvas, Double_t ymax, Double_t ymin) {
  outputCanvas->cd();
  // TLine *l_gapEvenBegin = new TLine(gapEvenBegin, ymin - 0.085*(ymax - ymin), gapEvenBegin, ymax + 0.085*(ymax - ymin));
  // TLine *l_gapEvenEnd = new TLine(gapEvenEnd, ymin - 0.085*(ymax - ymin), gapEvenEnd, ymax + 0.085*(ymax - ymin));
  // TLine *l_gapOddBegin = new TLine(gapOddBegin, ymin - 0.085*(ymax - ymin), gapOddBegin, ymax + 0.085*(ymax - ymin));
  // TLine *l_gapOddEnd = new TLine(gapOddEnd, ymin - 0.085*(ymax - ymin), gapOddEnd, ymax + 0.085*(ymax - ymin));
  // TLine *l_eefe = new TLine(27.5, ymin - 0.085*(ymax - ymin), 27.5, ymax + 0.085*(ymax - ymin));
  TLine *l_eefe = new TLine(27.5, ymin, 27.5, 0.8*ymax);
  TLine *l_febe = new TLine(39.5, ymin, 39.5, 0.8*ymax);

  // l_gapEvenBegin->SetLineColor(kRed);
  // l_gapEvenEnd->SetLineColor(kRed);
  // l_gapOddBegin->SetLineColor(kBlue);
  // l_gapOddEnd->SetLineColor(kBlue);
  
  // l_gapEvenBegin->Draw("SAME");
  // l_gapEvenEnd->Draw("SAME");
  // l_gapOddBegin->Draw("SAME");
  // l_gapOddEnd->Draw("SAME");
  l_eefe->Draw("SAME");
  l_febe->Draw("SAME");
}

void saveOccupiedCellsMapToFiles(TFile *outputFile, Double_t threshold, unsigned startLayer) {
  TString outputCanvasName = Form("occupiedCells_startLayer%i_threshold%.1f", startLayer, threshold);
  TCanvas *outputCanvas = new TCanvas(outputCanvasName, outputCanvasName);
  gPad->SetLogy(0);
  TLegend *legend = new TLegend(0.1, 0.1, 0.4, 0.5, "Colors for #eta-ranges:");
  for (std::vector<Int_t>::iterator etaBinsToPlotIterator = etaBinsToPlot.begin(); etaBinsToPlotIterator != etaBinsToPlot.end(); ++etaBinsToPlotIterator) {
    Int_t etaBin = *etaBinsToPlotIterator;
    // outputFile->WriteTObject(occupiedCellsMap[etaBin]);
    if (etaBin == etaBinsToPlot[0]) {
      gPad->SetLogy();
      TLegendEntry *legendEntry = legend->AddEntry(occupiedCellsMap[etaBin], namesForEtaBins[etaBin], "");
      legendEntry->SetTextColor(colorsForEtaBins[etaBin]);
      occupiedCellsMap[etaBin]->SetTitle(Form("Number of cells with deposited rechit energy above %.1f mips", threshold));
      occupiedCellsMap[etaBin]->GetXaxis()->SetTitle("Layer");
      occupiedCellsMap[etaBin]->SetLineColor(colorsForEtaBins[etaBin]);
      occupiedCellsMap[etaBin]->GetYaxis()->SetRangeUser(0.5, 1000);
      occupiedCellsMap[etaBin]->Draw("HIST LE");
      plotVerticalLayerBoundaries(outputCanvas, 1000, 0.5);
    }
    else {
      TLegendEntry *legendEntry = legend->AddEntry(occupiedCellsMap[etaBin], namesForEtaBins[etaBin], "");
      legendEntry->SetTextColor(colorsForEtaBins[etaBin]);
      occupiedCellsMap[etaBin]->SetLineColor(colorsForEtaBins[etaBin]);
      occupiedCellsMap[etaBin]->Draw("HIST SAME LE");
    }
  }
  legend->SetTextSize(0.05);
  legend->Draw();
  outputCanvas->SaveAs(plotsDir+versionName+Form("_")+outputCanvasName+Form(".png"));
  // outputCanvas->SaveAs(plotsDir+outputCanvasName+Form(".png"));
  outputFile->WriteTObject(outputCanvas);
}

void initializeOccupancyMaps(Double_t threshold, unsigned startLayer) {
  (void)(threshold);
  for (unsigned etaBinCounter = 0; etaBinCounter <= (1+nEtaBins); ++etaBinCounter) {
    // TString graphName = Form("occupancy_etaBinCounter%i_startLayer%i_threshold%.1f", etaBinCounter, startLayer, threshold);
    occupancyMap[etaBinCounter] = new TGraphErrors(nLayers-startLayer);
    // graphName = Form("occupancyMedian_etaBinCounter%i_startLayer%i_threshold%.1f", etaBinCounter, startLayer, threshold);
    occupancyMedianMap[etaBinCounter] = new TGraphErrors(nLayers-startLayer);
  }
}

void calculateOccupancies(Double_t threshold, unsigned startLayer) {
  std::cout << "Calculating occupancies..." << std::endl;
  for (Int_t etaBinCounter = 0; etaBinCounter <= (1+nEtaBins); ++etaBinCounter) {
    TProfile *occupiedCellsProfile = occupiedCellsMap[etaBinCounter];
    for (unsigned layerCounter = startLayer; layerCounter < nLayers; ++layerCounter) {
      Double_t occupancy = (occupiedCellsProfile->GetBinContent(1+layerCounter-startLayer))/(totalNumberOfCellsMap[layerCounter][etaBinCounter]);
      if (threshold == listOfThresholds[0] && etaBinCounter == listOfEtaBinsToPlot[0]) {
        if (occupancy > occupancyMaxval[startLayer]) occupancyMaxval[startLayer] = occupancy;
      }
      Double_t occupancyError = (occupiedCellsProfile->GetBinError(1+layerCounter-startLayer))/(totalNumberOfCellsMap[layerCounter][etaBinCounter]);
      occupancyMap[etaBinCounter]->SetPoint(layerCounter-startLayer, layerCounter, occupancy);
      occupancyMap[etaBinCounter]->SetPointError(layerCounter-startLayer, 0, occupancyError);

      std::vector<Int_t> occupiedCellsVector = occupiedCellsVectorMap[etaBinCounter][layerCounter];
      // Double_t *occupiedCellsArray = &occupiedCellsVector[0];
      Double_t occupancyMedian = TMath::Median(occupiedCellsVector.size(), &occupiedCellsVector[0])/(totalNumberOfCellsMap[layerCounter][etaBinCounter]);
      occupancyMedianMap[etaBinCounter]->SetPoint(layerCounter-startLayer, layerCounter, occupancyMedian);
      occupancyMedianMap[etaBinCounter]->SetPointError(layerCounter-startLayer, 0, occupancyError);
    }
  }
}

void saveOccupancyMapsToFiles(TFile *outputFile, Double_t threshold, unsigned startLayer, Bool_t plotMedian) {
  TString typePrefix = plotMedian ? Form("median_") : Form("mean_");
  TString outputCanvasName = typePrefix+Form("occupancyPlots_startLayer%i_threshold%.1f", startLayer, threshold);
  TCanvas *outputCanvas = new TCanvas(outputCanvasName, outputCanvasName);
  gPad->SetLogy(0);
  // TLegend *legend = new TLegend(0.7, 0.6, 0.9, 0.9, "Colors for #eta-ranges:");
  // TLegend *legend = new TLegend(0.8, 0.7, 1.0, 1.0, "Colors for #eta-ranges:");
  TLegend *legend = new TLegend(0.1, 0.1, 0.3, 0.4, "Colors for #eta-ranges:");
  for (std::vector<Int_t>::iterator etaBinsToPlotIterator = etaBinsToPlot.begin(); etaBinsToPlotIterator != etaBinsToPlot.end(); ++etaBinsToPlotIterator) {
    Int_t etaBin = *etaBinsToPlotIterator;
    TGraphErrors *occupancyMapToPlot = plotMedian ? occupancyMedianMap[etaBin] : occupancyMap[etaBin];
    // outputFile->WriteTObject(occupancyMapToPlot);
    if (etaBin == etaBinsToPlot[0] && startLayer == startLayers[0]) {
      gPad->SetLogy();
      TLegendEntry *legendEntry = legend->AddEntry(occupancyMapToPlot, namesForEtaBins[etaBin], "");
      legendEntry->SetTextColor(colorsForEtaBins[etaBin]);
      occupancyMapToPlot->SetMaximum(1.1*occupancyMaxval[startLayer]);
      TString altTypePrefix = plotMedian ? Form("Median "): Form("Mean ");
      occupancyMapToPlot->SetTitle(altTypePrefix + Form("Fraction of cells with deposited rechit energy above %.1f mips", threshold));
      occupancyMapToPlot->GetXaxis()->SetTitle("Layer");
      occupancyMapToPlot->SetLineColor(colorsForEtaBins[etaBin]);
      occupancyMapToPlot->GetYaxis()->SetRangeUser(0.00001, 1.1);
      occupancyMapToPlot->Draw();
      plotVerticalLayerBoundaries(outputCanvas, 1.1, 0.00001);
    }
    else {
      TLegendEntry *legendEntry = legend->AddEntry(occupancyMapToPlot, namesForEtaBins[etaBin], "");
      legendEntry->SetTextColor(colorsForEtaBins[etaBin]);
      occupancyMapToPlot->SetLineColor(colorsForEtaBins[etaBin]);
      occupancyMapToPlot->Draw("SAME");
    }
  }
  legend->SetTextSize(0.03);
  legend->Draw();
  outputCanvas->SaveAs(plotsDir+versionName+Form("_")+outputCanvasName+Form(".png"));
  // outputCanvas->SaveAs(plotsDir+typePrefix+outputCanvasName+Form(".png"));
  outputFile->WriteTObject(outputCanvas);
}

void calculateAndSaveOccupiedCellsVsThresholdToFile(TFile *outputFile, unsigned startLayer) {
  unsigned numberOfThresholds = sizeof(listOfThresholds)/sizeof(listOfThresholds[0]);
  TGraph *occupiedCellsVsThreshold = new TGraph(numberOfThresholds);
  for (unsigned thresholdCounter = 0; thresholdCounter < numberOfThresholds; ++thresholdCounter) {
    std::cout << "At threshold = " << listOfThresholds[thresholdCounter] << ", total number of occupied cells per event = " << (1./totalNEvents)*occupiedCellsVsThresholdMap[listOfThresholds[thresholdCounter]] << std::endl;
    // std::cout << "Checking... totalNEvents = " << totalNEvents << std::endl;
    occupiedCellsVsThreshold->SetPoint(thresholdCounter, listOfThresholds[thresholdCounter], (1./totalNEvents)*occupiedCellsVsThresholdMap[listOfThresholds[thresholdCounter]]);
  }
  
  TString outputCanvasName = Form("totalOccupiedCells_startLayer%i", startLayer);
  TCanvas *outputCanvas = new TCanvas(outputCanvasName, outputCanvasName);
  gPad->SetLogy(0);
  occupiedCellsVsThreshold->SetTitle(Form("Total number of occupied cells per event for layer number >= %i", startLayer));
  occupiedCellsVsThreshold->GetXaxis()->SetTitle("Threshold");
  occupiedCellsVsThreshold->Draw("AC*");
  outputCanvas->SaveAs(plotsDir+versionName+Form("_")+outputCanvasName+Form(".png"));
  outputFile->WriteTObject(outputCanvas);
}

void getPlotsFromTrees() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  gStyle->SetOptStat(0);

  initializeOccupancyMaxval();
  fillEtaBinEdges();
  initializeAxisWithEtaBinning();
  printAxisBinning(axisWithEtaBinning);
  initializeCalorimeterProperties();
  readTotalNumberOfCellsFromFiles();

  TFile *inputFile = TFile::Open(inputFileName);
  TFile *outputFile = new TFile(outputFileName, "RECREATE");
  unsigned numberOfThresholds = sizeof(listOfThresholds)/sizeof(listOfThresholds[0]);
  for (unsigned startLayerCounter = 0; startLayerCounter < (sizeof(startLayers)/sizeof(startLayers[0])); ++startLayerCounter) {
    unsigned startLayer = startLayers[startLayerCounter];
    for (unsigned thresholdCounter = 0; thresholdCounter < numberOfThresholds; ++thresholdCounter) {
      Double_t threshold = listOfThresholds[thresholdCounter];
      initializeOccupiedCellsMaps(threshold, startLayer);
      fillOccupiedCellsMapsFromData(inputFile, threshold, startLayer);
      // saveOccupiedCellsMapToFile(outputFile);
      saveOccupiedCellsMapToFiles(outputFile, threshold, startLayer);
      initializeOccupancyMaps(threshold, startLayer);
      calculateOccupancies(threshold, startLayer);
      saveOccupancyMapsToFiles(outputFile, threshold, startLayer, false);
      saveOccupancyMapsToFiles(outputFile, threshold, startLayer, true);
      occupiedCellsMap.clear();
      occupiedCellsVectorMap.clear();
      occupancyMap.clear();
      occupancyMedianMap.clear();
    }
    calculateAndSaveOccupiedCellsVsThresholdToFile(outputFile, startLayer);
    occupiedCellsVsThresholdMap.clear();
  }
  outputFile->Close();
  inputFile->Close();
}

# ifndef __CINT__
int main() {
  getPlotsFromTrees();
  return 0;
}
# endif
