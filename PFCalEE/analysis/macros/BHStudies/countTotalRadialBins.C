#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TAxis.h"
#include "Rtypes.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

#include "HGCSSInfo.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"
#include "HGCSSHardcodedConstants.hh"

void printAxisBinning(TAxis *inputAxis) {
  unsigned nBins = inputAxis->GetNbins();
  for (unsigned binCounter = 0; binCounter <= 1+nBins; ++binCounter) {
    std::cout << "Bin number " << binCounter << " ranges from " << inputAxis->GetBinLowEdge(binCounter) << " to " << inputAxis->GetBinUpEdge(binCounter) << " and its center is " << inputAxis->GetBinCenter(binCounter) << std::endl;
  }
}

bool testFile(std::string inputPath) {
  TFile *testFile = TFile::Open(inputPath.c_str());
  if ( !testFile ) {
    return false;
  }
  return true;
}

HGCSSGeometryConversion* getGeometryConversion(HGCSSDetector &myDetector, std::string inputFileName) {
  HGCSSInfo *info;
  double calorSizeXY = 0;
  double cellSize = 0;
  unsigned versionNumber = 0;
  unsigned model = 0;
  std::cout << "Initializing calorimeter properties..." << std::endl;
  // TString inputFileName = dataDir+Form("");
  bool dataExists = testFile(inputFileName);
  if (dataExists) {
    TFile *inputFile = TFile::Open(inputFileName.c_str());
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
  
  std::cout << "Initializing geometry conversion..." << std::endl;
  HGCSSGeometryConversion *geomConv = new HGCSSGeometryConversion(model,cellSize,false,2);
  geomConv->setXYwidth(calorSizeXY);
  geomConv->initialiseHoneyComb(calorSizeXY,cellSize);
  // geomConv->initialiseSquareMap(calorSizeXY,10.);
  std::string geometryInputFolder = "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/test/FH_BH_Geometry";
  geomConv->initialiseFHBHMaps(geometryInputFolder);
  return geomConv;
}

int getNLayers(HGCSSDetector &myDetector) {
  return myDetector.nLayers();
}

std::map<Int_t, Int_t> getNumberOfCellsMap(unsigned layerCounter, int nRadialBins, HGCSSGeometryConversion* geomConv, TAxis *axisWithRadialBinning) {
  std::cout << "Getting number of cells for layer index " << layerCounter << " ..." << std::endl;
  // Double_t zPosition = layerZPositions[layerCounter];

  // const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layerCounter);
  // bool isScint = subdet.isScint;
  
  std::map<Int_t, Int_t> numberOfCellsMap;
  for (int radialBinCounter = 0; radialBinCounter <= 1+nRadialBins; ++radialBinCounter) {
    numberOfCellsMap[radialBinCounter] = 0;
  }

  bool radialMapToBeUsed = (layerCounter >= ANNULARGEOMETRYFIRSTLAYER);

  std::map<int,std::pair<double,double> >& mapForLayer = radialMapToBeUsed?geomConv->fhbhGeoms[layerCounter-ANNULARGEOMETRYFIRSTLAYER]:geomConv->hexaGeom;
  unsigned nBins = mapForLayer.size();
  std::cout << "For this layer, nBins = " << nBins << std::endl;

  for (unsigned cellid = 1; cellid <= nBins; ++cellid){
    std::pair<double,double> xyPosition = mapForLayer[cellid];
    Double_t xPosition = xyPosition.first;
    Double_t yPosition = xyPosition.second;
    Double_t radialValue = sqrt(xPosition*xPosition + yPosition*yPosition);
    int correspondingRadialBin = axisWithRadialBinning->FindFixBin(radialValue);
    if (layerCounter >= ANNULARGEOMETRYFIRSTLAYER && (correspondingRadialBin <= 0 || correspondingRadialBin >= 1+nRadialBins)) {
      std::cout << "Getting a cell outside radial binning, with cell ID: " << cellid << ", x = " << xPosition << ", y = " << yPosition << std::endl;
      std::exit(EXIT_FAILURE);
    }
    numberOfCellsMap[correspondingRadialBin] += 1;
  }
  return numberOfCellsMap;
}

void writeNumberOfCellsToFiles(const int &nLayers, const int &nRadialBins, HGCSSGeometryConversion *geomConv, const std::string &outputDir, TAxis *axisWithRadialBinning) {
  TString TOutputDir = Form(outputDir.c_str());
  std::cout << "Writing number of cells to files..." << std::endl;
  if (system(("mkdir -p " + outputDir + " && rm " + outputDir + "/*").c_str())) {
    std::cout << "Output directory now fresh!" << std::endl;
  }
  else {
    std::cout << "Error, unable to create a new directory or empty the existing directory: " << outputDir << std::endl;
  }
  for (int layerCounter = 0; layerCounter < nLayers; ++layerCounter) {
    ofstream outputFile;
    outputFile.open((TOutputDir+Form("/totalNumberOfCellsInRadialBins_layer%i.dat", layerCounter)).Data());
    std::map<Int_t, Int_t> numberOfCellsMap = getNumberOfCellsMap(layerCounter, nRadialBins, geomConv, axisWithRadialBinning);
    for (std::map<Int_t, Int_t>::iterator numberOfCellsMapIterator = numberOfCellsMap.begin(); numberOfCellsMapIterator != numberOfCellsMap.end(); ++numberOfCellsMapIterator) {
      outputFile << numberOfCellsMapIterator->first << "    " << numberOfCellsMapIterator->second << std::endl;
    }
    outputFile.close();
  }
}

void countTotalBins(double minR, double maxR, int nRadialBins, std::string inputFileName, std::string outputDir) {
  // load the shared library for HGCSS* classes:
  // gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  TAxis *axisWithRadialBinning = new TAxis(nRadialBins, minR, maxR);
  printAxisBinning(axisWithRadialBinning);
  HGCSSDetector &myDetector = theDetector();
  HGCSSGeometryConversion *geomConv = getGeometryConversion(myDetector, inputFileName);
  int nLayers = getNLayers(myDetector);
  writeNumberOfCellsToFiles(nLayers, nRadialBins, geomConv, outputDir, axisWithRadialBinning);
  delete axisWithRadialBinning;
  delete geomConv;
}

# ifndef __CINT__
int main(int argc, char ** argv) {
  int nExpectedArguments = 5;
  if (argc != 1+nExpectedArguments) {
    std::cout << "Error: number of arguments expected: " << nExpectedArguments << "; number of arguments supplied: " << argc-1 << std::endl;
    std::exit(EXIT_FAILURE);
  }
  double minR = std::atof(argv[1]);
  double maxR = std::atof(argv[2]);
  int nRadialBins = std::atoi(argv[3]);
  std::string inputFileName = std::string(argv[4]);
  std::string outputDir = std::string(argv[5]);
  countTotalBins(minR, maxR, nRadialBins, inputFileName, outputDir);
  return 0;
}
# endif
