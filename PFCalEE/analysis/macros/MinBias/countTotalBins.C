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

// const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");
const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_SingleEvents/");
const std::string inputEtaRangesFileName = "etaBinBoundaries.dat";
const std::string inputLayerZPositionsFileName = "layerZPositions.dat";
const TString outputDir = Form("totalNumberOfCellsInEtaBinsInfo");

HGCSSInfo *info;
double calorSizeXY = 0;
double cellSize = 0;
unsigned versionNumber = 0;
unsigned model = 0;
unsigned nLayers;

HGCSSDetector & myDetector = theDetector();

HGCSSGeometryConversion *geomConv;

Double_t *etaBinEdges;
Double_t *layerZPositions;
unsigned nEtaBins;

TAxis *axisWithEtaBinning;

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
  TString digiFileName = dataDir+Form("Digi_Pu200_IC3_version33_Minbias_14TeV_SingleEvents_0001.root");
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

void fillLayerZPositions() {
  std::vector<Double_t> zPositions;
  std::string line;
  std::ifstream inputFile(inputLayerZPositionsFileName.c_str());
  if (inputFile.is_open()) {
    while(std::getline(inputFile, line)) {
      zPositions.push_back(atof(line.c_str()));
    }
    inputFile.close();
  }
  else {
    std::cout << "ERROR: Unable to open input file to read in layer z positions" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  layerZPositions = new Double_t[zPositions.size()];
  for (unsigned layerCounter = 0; layerCounter < zPositions.size(); ++layerCounter) {
    layerZPositions[layerCounter] = zPositions[layerCounter];
    std::cout << "z position [" << layerCounter << "] = " << layerZPositions[layerCounter] << std::endl;
  }
  if (zPositions.size() != nLayers) {
    std::cout << "Something has gone terribly wrong... number of layers present in input file for z positions of layers is not the same as the number of layers obtained from initializing detector with this design" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  zPositions.clear();
}

void initializeGeometryConversion() {
  std::cout << "Initializing geometry conversion..." << std::endl;
  geomConv = new HGCSSGeometryConversion(model,cellSize,false,2);
  geomConv->setXYwidth(calorSizeXY);
  geomConv->initialiseHoneyComb(calorSizeXY,cellSize);
  // geomConv->initialiseSquareMap(calorSizeXY,10.);
  std::string geometryInputFolder = "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/test/FH_BH_Geometry";
  geomConv->initialiseFHBHMaps(geometryInputFolder);
}

inline Double_t getEta(Double_t xPosition, Double_t yPosition, Double_t zPosition) {
  return (ROOT::Math::XYZPoint(xPosition/10.,yPosition/10.,zPosition/10.)).eta();
}

std::map<Int_t, Int_t> getNumberOfCellsMap(unsigned layerCounter) {
  std::cout << "Getting number of cells for layer index " << layerCounter << " ..." << std::endl;
  Double_t zPosition = layerZPositions[layerCounter];

  // const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layerCounter);
  // bool isScint = subdet.isScint;
  
  std::map<Int_t, Int_t> numberOfCellsMap;
  for (unsigned etaBinCounter = 0; etaBinCounter <= 1+nEtaBins; ++etaBinCounter) {
    numberOfCellsMap[etaBinCounter] = 0;
  }

  bool radialMapToBeUsed = (layerCounter >= ANNULARGEOMETRYFIRSTLAYER);

  std::map<int,std::pair<double,double> >& mapForLayer = radialMapToBeUsed?geomConv->fhbhGeoms[layerCounter-ANNULARGEOMETRYFIRSTLAYER]:geomConv->hexaGeom;
  unsigned nBins = mapForLayer.size();
  std::cout << "For this layer, nBins = " << nBins << std::endl;

  for (unsigned cellid = 1; cellid <= nBins; ++cellid){
    std::pair<double,double> xyPosition = mapForLayer[cellid];
    Double_t xPosition = xyPosition.first;
    Double_t yPosition = xyPosition.second;
    Double_t etaValue = getEta(xPosition, yPosition, zPosition);
    unsigned correspondingEtaBin = axisWithEtaBinning->FindFixBin(etaValue);
    numberOfCellsMap[correspondingEtaBin] += 1;
  }
  return numberOfCellsMap;
}

void writeNumberOfCellsToFiles() {
  std::cout << "Writing number of cells to files..." << std::endl;
  for (unsigned layerCounter = 0; layerCounter < nLayers; ++layerCounter) {
    ofstream outputFile;
    outputFile.open((outputDir+Form("/totalNumberOfCellsInEtaBins_layer%i.dat", layerCounter)).Data());
    std::map<Int_t, Int_t> numberOfCellsMap = getNumberOfCellsMap(layerCounter);
    for (std::map<Int_t, Int_t>::iterator numberOfCellsMapIterator = numberOfCellsMap.begin(); numberOfCellsMapIterator != numberOfCellsMap.end(); ++numberOfCellsMapIterator) {
      outputFile << numberOfCellsMapIterator->first << "    " << numberOfCellsMapIterator->second << std::endl;
    }
    outputFile.close();
  }
}

void countTotalBins() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  fillEtaBinEdges();
  initializeAxisWithEtaBinning();
  printAxisBinning(axisWithEtaBinning);
  initializeCalorimeterProperties();
  fillLayerZPositions();
  initializeGeometryConversion();
  writeNumberOfCellsToFiles();
}

# ifndef __CINT__
int main() {
  countTotalBins();
  return 0;
}
# endif
