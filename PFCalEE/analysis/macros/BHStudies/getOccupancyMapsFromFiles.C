#include <fstream>
#include <cstdlib>
#include <iostream>
#include <map>
#include <ctime>

#include "TAxis.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"

#include "HGCSSHardcodedConstants.hh"

void initializeHistos(TH2F *h_occupancy, TH2F *h_occupancy_errors) {
  int nbins_x = h_occupancy->GetNbinsX();
  int nbins_y = h_occupancy->GetNbinsY();
  for (int binxcounter = 0; binxcounter <= (1+nbins_x); ++binxcounter) {
    for (int binycounter = 0; binycounter <= (1+nbins_y); ++binycounter) {
      h_occupancy->SetBinContent(binxcounter, binycounter, 0.);
      h_occupancy_errors->SetBinContent(binxcounter, binycounter, 0.);
    }
  }
}

template<typename typeOfData>
std::vector<typeOfData> getMapFromFile(std::string inputFileName) {
  std::vector<typeOfData> dataFromFile;
  std::ifstream inputFile;
  inputFile.open(inputFileName);
  if (inputFile.is_open()) {
    // std::cout << "Opened file to read in data!" << std::endl;
    int radialBin;
    typeOfData dataValue;
    while (inputFile >> radialBin >> dataValue) {
      dataFromFile.push_back(dataValue);
    }
  }
  else {
    std::cout << "Could not open " << inputFileName << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return dataFromFile;
}

template<typename typeOfData>
void printVector(std::vector<typeOfData> inputVector, std::string vectorName) {
  std::cout << "Vector name: " << vectorName << std::endl;
  typename std::vector<typeOfData>::iterator inputVectorIterator = inputVector.begin();
  for (; inputVectorIterator != inputVector.end(); ++inputVectorIterator) {
    std::cout << *inputVectorIterator << std::endl;
  }
  std::cout << "Size of vector: " << inputVector.size() << std::endl;
}

void fillHistosFromFiles(TH2F *h_occupancy, TH2F *h_occupancy_errors, std::string totalNumberOfCellsDir, std::string occupiedNumberOfCellsDir, double threshold, int nRadialBins, int startLayer, int stopLayer) {
  std::cout << "Starting to fill histograms from files..." << std::endl;
  for (int layer = startLayer; layer <= stopLayer; ++layer) {
    Int_t layerBinNumber = h_occupancy->GetXaxis()->FindFixBin(double(layer));
    std::string totalNumberOfCellsInputFileName = totalNumberOfCellsDir + std::string(Form("/totalNumberOfCellsInRadialBins_layer%d.dat", layer));
    std::vector<int> totalNumberOfCellsList = getMapFromFile<int>(totalNumberOfCellsInputFileName);
    std::string occupiedCellsInputFileName = occupiedNumberOfCellsDir + std::string(Form("/threshold_%.1f/occupiedCellsInRadialBins_layer%d.dat", threshold, layer));
    std::vector<double> occupiedCellsList = getMapFromFile<double>(occupiedCellsInputFileName);
    std::string occupiedCellsErrorsInputFileName = occupiedNumberOfCellsDir + std::string(Form("/threshold_%.1f/occupiedCellsErrorsInRadialBins_layer%d.dat", threshold, layer));
    std::vector<double> occupiedCellsErrorsList = getMapFromFile<double>(occupiedCellsErrorsInputFileName);
    if ((int(totalNumberOfCellsList.size()) != 2+nRadialBins) || (int(occupiedCellsList.size()) != 2+nRadialBins) || (int(occupiedCellsErrorsList.size()) != 2+nRadialBins)) {
      std::cout << "Something has gone terribly wrong... an unexpected number of radial bins found !" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for (int radialBinCounter = 1; radialBinCounter <= (nRadialBins); ++radialBinCounter) {
      if (totalNumberOfCellsList[radialBinCounter] > 0) {
        h_occupancy->SetBinContent(layerBinNumber, radialBinCounter, occupiedCellsList[radialBinCounter]/(1.0*totalNumberOfCellsList[radialBinCounter]));
        h_occupancy_errors->SetBinContent(layerBinNumber, radialBinCounter, occupiedCellsErrorsList[radialBinCounter]/(1.0*totalNumberOfCellsList[radialBinCounter]));
      }
    }
  }
}

void plotHistogram(TH2F *histogramToPlot, std::string outputDir, std::string histogramName, double threshold) {
  std::cout << "Saving histogram " << histogramName << std::endl;
  TCanvas *c_hist = new TCanvas(histogramToPlot->GetName(), histogramToPlot->GetTitle(), 1024, 768);
  gStyle->SetOptStat(0);
  gPad->SetLogz();
  histogramToPlot->Draw("colz");
  c_hist->SaveAs((outputDir + "/" + histogramName + std::string(Form("_threshold%.1f.png", threshold))).c_str());
  delete c_hist;
}

void getOccupancyMapsFromFiles(double minR, double maxR, int nRadialBins, std::string totalNumberOfCellsDir, std::string occupiedNumberOfCellsDir, std::string outputDir, double threshold) {
  std::cout << "Fetching occupancy statistics from files..." << std::endl;
  
  TH2F *h_occupancy = new TH2F("occupancy", "Occupancy;Layer;r(mm)", NFHBHLAYERS, -0.5+ANNULARGEOMETRYFIRSTLAYER,0.5+BHLASTLAYER, nRadialBins, minR, maxR);
  TH2F *h_occupancy_errors = new TH2F("occupancyErrors", "Occupancy Error;Layer;r(mm)", NFHBHLAYERS, -0.5+ANNULARGEOMETRYFIRSTLAYER,0.5+BHLASTLAYER, nRadialBins, minR, maxR);
  initializeHistos(h_occupancy, h_occupancy_errors);
  fillHistosFromFiles(h_occupancy, h_occupancy_errors, totalNumberOfCellsDir, occupiedNumberOfCellsDir, threshold, nRadialBins, ANNULARGEOMETRYFIRSTLAYER, BHLASTLAYER);

  if (system(("mkdir -p " + outputDir).c_str())) {
    std::cout << "Unable to create output directory: " << outputDir << std::endl;
    std::exit(EXIT_FAILURE);
  }
  plotHistogram(h_occupancy, outputDir, "occupancyProfile", threshold);
  delete h_occupancy;
  plotHistogram(h_occupancy_errors, outputDir, "occupancyErrorsProfile", threshold);
  delete h_occupancy_errors;
}

# ifndef __CINT__
int main(int argc, char ** argv) {
  int nExpectedArguments = 7;
  if (argc != 1+nExpectedArguments) {
    std::cout << "Error: number of arguments expected: " << nExpectedArguments << "; number of arguments supplied: " << argc-1 << std::endl;
    std::exit(EXIT_FAILURE);
  }
  double minR = std::atof(argv[1]);
  double maxR = std::atof(argv[2]);
  int nRadialBins = std::atoi(argv[3]);
  std::string totalNumberOfCellsDir = std::string(argv[4]);
  std::string occupiedNumberOfCellsDir = std::string(argv[5]);
  std::string outputDir = std::string(argv[6]);
  double threshold = std::atof(argv[7]);
  std::cout << "Arguments supplied: " << std::endl
            << "minR = " << minR << std::endl
            << "maxR = " << maxR << std::endl
            << "nRadialBins = " << nRadialBins << std::endl
            << "totalNumberOfCellsDir = " << totalNumberOfCellsDir << std::endl
            << "occupiedNumberOfCellsDir = " << occupiedNumberOfCellsDir << std::endl
            << "outputDir = " << outputDir << std::endl
            << "threshold = " << threshold << std::endl
            << std::endl;
  getOccupancyMapsFromFiles(minR, maxR, nRadialBins, totalNumberOfCellsDir, occupiedNumberOfCellsDir, outputDir, threshold);
  return 0;
}
# endif
