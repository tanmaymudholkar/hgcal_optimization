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

#define CONSTNOISE 0.2

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

HGCSSGeometryConversion* getGeometryConversion(HGCSSDetector &myDetector, std::string inputFileName, std::string geometryInputFolder) {
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
  // std::string geometryInputFolder = "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/test/FH_BH_Geometry";
  geomConv->initialiseFHBHMaps(geometryInputFolder);
  return geomConv;
}

int getNLayers(HGCSSDetector &myDetector) {
  return myDetector.nLayers();
}

std::map<Int_t, Double_t> getNoiseEstimateMap(unsigned layerCounter, int nRadialBins, HGCSSGeometryConversion* geomConv, TAxis *axisWithRadialBinning, std::map<int, double> variableNoiseForLayer, TAxis *variableBinningAxisForLayer, double threshold) {
  std::cout << "Getting number of cells for layer index " << layerCounter << " ..." << std::endl;
  std::map<Int_t, Double_t> noiseEstimateMap;
  for (int radialBinCounter = 0; radialBinCounter <= 1+nRadialBins; ++radialBinCounter) {
    noiseEstimateMap[radialBinCounter] = 0.;
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
    int radialBinNumberNoiseBinning = (variableBinningAxisForLayer)->FindFixBin(radialValue);
    double noiseSigmaVal = variableNoiseForLayer[radialBinNumberNoiseBinning];
    if (layerCounter >= ANNULARGEOMETRYFIRSTLAYER && (correspondingRadialBin <= 0 || correspondingRadialBin >= 1+nRadialBins)) {
      std::cout << "Getting a cell outside radial binning, with cell ID: " << cellid << ", x = " << xPosition << ", y = " << yPosition << std::endl;
      std::exit(EXIT_FAILURE);
    }
    noiseEstimateMap[correspondingRadialBin] += 0.5*erfc(threshold/(noiseSigmaVal*sqrt(2)));
  }
  return noiseEstimateMap;
}

void writeNoiseEstimatesToFiles(const int &nLayers, const int &nRadialBins, HGCSSGeometryConversion *geomConv, const std::string &outputDir, TAxis *axisWithRadialBinning, std::map<unsigned, std::map<int, double> > variableNoise, std::map<unsigned, TAxis*> variableBinningAxis, double threshold) {
  std::cout << "Writing number of cells to files..." << std::endl;
  std::string outputDirWithThreshold = outputDir + std::string(Form("/threshold_%.1f", threshold));
  if (system(("mkdir -p " + outputDirWithThreshold + " && rm -rf " + outputDirWithThreshold + "/*").c_str()) != 0) {
    std::cout << "Error, unable to create a new directory or empty the existing directory: " << outputDirWithThreshold << std::endl;
  }
  else {
    std::cout << "Output directory now fresh!" << std::endl;
  }
  for (int layerCounter = 0; layerCounter < nLayers; ++layerCounter) {
    ofstream outputFile;
    ofstream errorsOutputFile;
    outputFile.open(outputDirWithThreshold+std::string(Form("/occupiedCellsInRadialBins_layer%i.dat", layerCounter)));
    errorsOutputFile.open(outputDirWithThreshold+std::string(Form("/occupiedCellsErrorsInRadialBins_layer%i.dat", layerCounter)));
    if (layerCounter >= ANNULARGEOMETRYFIRSTLAYER) {
      std::map<Int_t, Double_t> noiseEstimateMap = getNoiseEstimateMap(layerCounter, nRadialBins, geomConv, axisWithRadialBinning, variableNoise[layerCounter], variableBinningAxis[layerCounter], threshold);
      for (std::map<Int_t, Double_t>::iterator noiseEstimateMapIterator = noiseEstimateMap.begin(); noiseEstimateMapIterator != noiseEstimateMap.end(); ++noiseEstimateMapIterator) {
        outputFile << noiseEstimateMapIterator->first << "    " << noiseEstimateMapIterator->second << std::endl;
        errorsOutputFile << noiseEstimateMapIterator->first << "    0.0" << std::endl;
      }
    }
    else {
      for (unsigned radialBinCounter = 0; radialBinCounter <= (1+nRadialBins); ++radialBinCounter) {
        outputFile << radialBinCounter << "    0.0" << std::endl;
        errorsOutputFile << radialBinCounter << "    0.0" << std::endl;
      }
    }
    outputFile.close();
    errorsOutputFile.close();
  }
}

void setVariableNoise(const unsigned & alay, const std::string inputGeometryFilePath, std::map<unsigned, std::map<int, double> > & variableNoise, std::map<unsigned, TAxis*> & variableBinningAxis, bool constantNoiseSwitch) {
  std::cout << "Starting to set noise for layer " << alay << "... reading in geometry..." << std::endl;
  std::vector<Double_t> radialBinEdgesVector;
  std::vector<Double_t> mipNoiseVector;
  std::ifstream inputFile;
  inputFile.open(inputGeometryFilePath);
  if (inputFile.is_open()) {
    std::cout << "Opened file! Now reading in geometry..." << std::endl;
    Int_t layerNumber;
    Double_t centerR, area, mipsig, sipm_noise, sig_over_noise, power, fluence, dose, outerR, innerR;
    bool isFirstLine = true;
    while (inputFile >> layerNumber >> centerR >> area >> mipsig >> sipm_noise >> sig_over_noise >> power >> fluence >> dose >> outerR >> innerR) {
      std::cout << "layerNumber: " << layerNumber
                << "; centerR = " << centerR
                << "; area = " << area
                << "; mipsig = " << mipsig
                << "; sipm_noise = " << sipm_noise
                << "; sig_over_noise = " << sig_over_noise
                << "; power = " << power
                << "; fluence = " << fluence
                << "; dose = " << dose
                << "; outerR = " << outerR
                << "; innerR = " << innerR << std::endl;
      
      radialBinEdgesVector.push_back(innerR);
      mipNoiseVector.push_back(1./sig_over_noise);
      if (isFirstLine) {
        isFirstLine = false;
      }
    }
    radialBinEdgesVector.push_back(outerR);
  }
  else {
    std::cout << "Could not open file " << inputGeometryFilePath << " for reading!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  inputFile.close();
  std::cout << "Successfully read geometry from file! Now setting binning..." << std::endl;
  int nAxisBins = radialBinEdgesVector.size()-1;
  Double_t *radialBinEdges = &radialBinEdgesVector[0];
  variableBinningAxis[alay] = new TAxis(nAxisBins, radialBinEdges);
  std::map<int, double> noiseMap;
  noiseMap[0] = 0.;
  for (int binCounter = 1; binCounter <= nAxisBins; ++binCounter) {
    noiseMap[binCounter] = (constantNoiseSwitch? CONSTNOISE: mipNoiseVector[binCounter-1]);
  }
  noiseMap[1+nAxisBins] = 0.;
  variableNoise[alay] = noiseMap;
  std::cout << "Checking..." << std::endl;
  std::cout << "At layer = " << alay << ":" << std::endl;
  for (std::map<int, double>::iterator noiseMapIterator = noiseMap.begin(); noiseMapIterator != noiseMap.end(); ++noiseMapIterator) {
    std::cout << "noise[" << noiseMapIterator->first << "] = " << noiseMapIterator->second << std::endl;
  }
}

void setVariableNoiseForAllLayers(std::string geometryInputFolder, std::map<unsigned, std::map<int, double> > & variableNoise, std::map<unsigned, TAxis*> & variableBinningAxis, bool constantNoiseSwitch) {
  for (unsigned layerCounter = ANNULARGEOMETRYFIRSTLAYER; layerCounter <= BHLASTLAYER; ++layerCounter) {
    std::string inputGeometryFilePath;
    unsigned offset;
    if (layerCounter >= ANNULARGEOMETRYFIRSTLAYER && layerCounter <= FHLASTLAYER) {
      offset = FHFIRSTLAYER-1;
      inputGeometryFilePath = geometryInputFolder + "/" + "geometry_FH" + std::to_string(layerCounter-offset) + ".txt";
    }
    else if (layerCounter >= BHFIRSTLAYER && layerCounter <= BHLASTLAYER) {
      offset = BHFIRSTLAYER-1;
      inputGeometryFilePath = geometryInputFolder + "/" + "geometry_BH" + std::to_string(layerCounter-offset) + ".txt";
    }
    setVariableNoise(layerCounter, inputGeometryFilePath, variableNoise, variableBinningAxis, constantNoiseSwitch);
  }
}

void getNoiseEstimate(double minR, double maxR, int nRadialBins, std::string inputFileName, std::string outputDir, bool constantNoiseSwitch, double threshold) {
  TAxis *axisWithRadialBinning = new TAxis(nRadialBins, minR, maxR);
  printAxisBinning(axisWithRadialBinning);
  HGCSSDetector &myDetector = theDetector();
  std::string geometryInputFolder = "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/test/FH_BH_Geometry";
  HGCSSGeometryConversion *geomConv = getGeometryConversion(myDetector, inputFileName, geometryInputFolder);
  int nLayers = getNLayers(myDetector);
  std::map<unsigned, std::map<int, double> > variableNoise;
  std::map<unsigned, TAxis*> variableBinningAxis;
  setVariableNoiseForAllLayers(geometryInputFolder, variableNoise, variableBinningAxis, constantNoiseSwitch);
  writeNoiseEstimatesToFiles(nLayers, nRadialBins, geomConv, outputDir, axisWithRadialBinning, variableNoise, variableBinningAxis, threshold);
  delete axisWithRadialBinning;
  delete geomConv;
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
  std::string inputFileName = std::string(argv[4]);
  std::string outputDir = std::string(argv[5]);
  std::string noiseType = std::string(argv[6]);
  bool constantNoiseSwitch;
  if (noiseType == "variable") {
    constantNoiseSwitch = false;
  }
  else if (noiseType == "constant") {
    constantNoiseSwitch = true;
  }
  else {
    std::cout << "Only variable or constant noise profiles accepted!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  double threshold = std::atof(argv[7]);
  getNoiseEstimate(minR, maxR, nRadialBins, inputFileName, outputDir, constantNoiseSwitch, threshold);
  return 0;
}
# endif
