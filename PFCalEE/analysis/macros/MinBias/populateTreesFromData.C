#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TStyle.h"
#include "TAxis.h"
#include "Rtypes.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "Math/Point3D.h"
#include "Math/Point3Dfwd.h"

// #include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/Minbias_14TeV_10Events/");
const TString inputFileNamePrefix = Form("Digi_Pu200_IC3_version33_model2_Minbias_14TeV_10Events_");
const std::string inputEtaRangesFileName = "etaBinBoundaries.dat";
const TString totalNumberOfCellsInEtaBinsDir = Form("totalNumberOfCellsInEtaBinsInfo");
const TString outputDir = Form("/afs/cern.ch/work/t/tmudholk/public/bhStudies/occupancyRootFiles/");
const TString outputFileNamePrefix = Form("occupiedCellsTrees_Minbias_14TeV_10Events");
// const TString versionName = Form("noXTalk");

// const Double_t listOfThresholds[] = {0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 5.5};
const Double_t listOfThresholds[] = {0.5, 5.0};
// const Double_t listOfThresholds[] = {0.5, 5.0, 5.5};
// const Double_t listOfThresholds[] = {0.5, 2.0, 5.0}; // assumed to be in ascending order of threshold in mips
unsigned hardcodedMaxDataCounter = 0;

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

std::map<unsigned, std::map<unsigned, Int_t> > totalNumberOfCellsMap; // first index: layer number, second index: eta bin number

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
  TString digiFileName = dataDir + inputFileNamePrefix + Form("%04d", 1) + Form(".root");
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

void initializeGeometryConversion() {
  std::cout << "Initializing geometry conversion..." << std::endl;
  geomConv = new HGCSSGeometryConversion(model,cellSize,false,2);
  geomConv->setXYwidth(calorSizeXY);
  geomConv->initialiseHoneyComb(calorSizeXY,cellSize);
  // geomConv->initialiseSquareMap(calorSizeXY,10.);
  std::string geometryInputFolder = "/afs/cern.ch/user/t/tmudholk/public/research/hgcal_BHStudies/test/FH_BH_Geometry";
  geomConv->initialiseFHBHMaps(geometryInputFolder);
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

void readInputFiles(TChain *lRecTree) {
  // TString digiPrefix = Form("DigiPu200_IC3_version33_000");
  TString suffix = Form(".root");
  unsigned dataCounter = 0;
  unsigned consecutiveNonexistentFiles = 0;
  // while (consecutiveNonexistentFiles < 4 && (hardcodedMaxDataCounter == 0 || dataCounter < hardcodedMaxDataCounter)) {
  bool continueReading = true;
  // while (consecutiveNonexistentFiles < 4) {
  while (continueReading) {
    std::cout << "Reading input files for data counter " << dataCounter << std::endl;
    // TFile *inputFile;
    TString fileName = dataDir + inputFileNamePrefix + Form("%04d", dataCounter) + suffix;
    bool dataExists = testFile(fileName);
    if (dataExists) {
      lRecTree->AddFile(fileName);
      consecutiveNonexistentFiles = 0;
      if (dataCounter == 0) {
        TFile *inputFile = TFile::Open(fileName);
        info=(HGCSSInfo*)inputFile->Get("Info");
        calorSizeXY = info->calorSizeXY();
        cellSize = info->cellSize();
        versionNumber = info->version();
        model = info->model();
        // inputFile->Close();

        std::cout << " -- Calor size XY = " << calorSizeXY
                  << ", version number = " << versionNumber 
                  << ", model = " << model
                  << ", cellSize = " << cellSize
                  << std::endl;
      }
    }
    else {
      ++consecutiveNonexistentFiles;
    }
    ++dataCounter;
    continueReading = (consecutiveNonexistentFiles < 4);
    if (hardcodedMaxDataCounter != 0) continueReading = (continueReading && dataCounter < hardcodedMaxDataCounter);
  }
}

void printETA(Double_t timeRemainingInSec, unsigned ievt, unsigned nEvts_total) {
  Int_t hoursRemaining = static_cast<Int_t>(timeRemainingInSec/3600.);
  Int_t minutesRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining)/60.);
  Int_t secondsRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining - 60*minutesRemaining) + 0.5);
  
  std::cout << Form("\rApproximately %02d h: %02d m: %02d s remaining. Analyzing %04d/%04d", hoursRemaining, minutesRemaining, secondsRemaining, ievt, nEvts_total) << std::flush;
}

// void fillTree(std::vector<HGCSSRecoHit> *rechitvec, TTree *outputTree, Double_t &energy, unsigned &cellid, Double_t &eta, Double_t &xpos, Double_t &ypos, unsigned &layer) {
void fillTree(std::vector<HGCSSRecoHit> *rechitvec, TTree *outputTree, Int_t &etaBin, unsigned &layer, Int_t &occupiedCells, Double_t threshold) {
  std::map<unsigned, std::map<unsigned, Int_t> > occupiedCellsMap;
  for (unsigned layerCounter = 0; layerCounter < nLayers; ++layerCounter) {
    for (unsigned etaBinCounter = 0; etaBinCounter <= (1+nEtaBins); ++etaBinCounter) {
      occupiedCellsMap[layerCounter][etaBinCounter] = 0;
    }
  }
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
    Double_t energy = ((*rechitvec)[iH]).energy();
    if (energy > threshold) {
      // cellid = ((*rechitvec)[iH]).cellid();
      Double_t eta = ((*rechitvec)[iH]).eta();
      unsigned etaBinForRecHit = axisWithEtaBinning->FindFixBin(eta);
      // xpos = ((*rechitvec)[iH]).get_x();
      // ypos = ((*rechitvec)[iH]).get_y();
      unsigned layerForRecHit = ((*rechitvec)[iH]).layer();
      occupiedCellsMap[layerForRecHit][etaBinForRecHit] += 1;
    }
  }
  for (unsigned layerCounter = 0; layerCounter < nLayers; ++layerCounter) {
    layer = layerCounter;
    for (unsigned etaBinCounter = 0; etaBinCounter <= (1+nEtaBins); ++etaBinCounter) {
      etaBin = etaBinCounter;
      occupiedCells = occupiedCellsMap[layerCounter][etaBinCounter];
      outputTree->Fill();
    }
  }
}

// void readDataIntoTree(TTree *outputTree, Double_t &energy, unsigned &cellid, Double_t &eta, Double_t &xpos, Double_t &ypos, unsigned &layer) {
void readDataIntoTree(TTree *outputTree, Int_t &etaBin, unsigned &layer, Int_t &occupiedCells, Double_t threshold) {
  TChain  *lRecTree = new TChain("RecoTree");
  readInputFiles(lRecTree);
  std::vector<HGCSSRecoHit> *rechitvec = 0;
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  unsigned nEvts_total = lRecTree->GetEntries();
  std::cout << "Total number of available events = " << nEvts_total << std::endl;
  Double_t fractionCompleted;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  std::time_t initialTimeClocked = std::time(nullptr);
  lRecTree->GetEntry(0);
  fillTree(rechitvec, outputTree, etaBin, layer, occupiedCells, threshold);
  for (unsigned ievt(1); ievt<nEvts_total; ++ievt){// loop on events
    fractionCompleted = float(ievt*1./nEvts_total);
    timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
    guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
    printETA(guess_timeRemainingInSec, ievt, nEvts_total);
        
    lRecTree->GetEntry(ievt);
    fillTree(rechitvec, outputTree, etaBin, layer, occupiedCells, threshold);
  }
  std::cout << std::endl;
}

// void createBranches(TTree *outputTree, Double_t &energy, unsigned &cellid, Double_t &eta, Double_t &xpos, Double_t &ypos, unsigned &layer) {
void createBranches(TTree *outputTree, Int_t &etaBin, unsigned &layer, Int_t &occupiedCells) {
  // TBranch *energyBranch = outputTree->Branch("energy", &energy);
  // TBranch *cellidBranch = outputTree->Branch("cellid", &cellid);
  // TBranch *etaBranch = outputTree->Branch("eta", &eta);
  // TBranch *xposBranch = outputTree->Branch("xpos", &xpos);
  // TBranch *yposBranch = outputTree->Branch("ypos", &ypos);
  // TBranch *layerBranch = outputTree->Branch("layer", &layer);

  outputTree->Branch("etaBin", &etaBin);
  outputTree->Branch("layer", &layer);
  outputTree->Branch("occupiedCells", &occupiedCells);

  // (void)(energyBranch);
  // (void)(cellidBranch);
  // (void)(etaBranch);
  // (void)(xposBranch);
  // (void)(yposBranch);
  // (void)(layerBranch);
}

// void populateTreesFromData(Double_t &threshold) {
void populateTreesFromData() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");
  
  // Create addresses for variables to put in corresponding branches
  // Double_t energy = 0;
  // unsigned cellid = 0;
  // Double_t eta = 0;
  // Double_t xpos = 0;
  // Double_t ypos = 0;
  // unsigned layer = 0;

  fillEtaBinEdges();
  initializeAxisWithEtaBinning();
  printAxisBinning(axisWithEtaBinning);
  initializeCalorimeterProperties();
  // initializeGeometryConversion();
  
  Int_t etaBin = 0;
  unsigned layer = 0;
  Int_t occupiedCells = 0;

  // TFile *outputFile = TFile::Open(outputDir+versionName+Form("_")+outputFileNamePrefix+Form(".root"), "RECREATE");
  TFile *outputFileInit = TFile::Open(outputDir+outputFileNamePrefix+Form(".root"), "RECREATE");
  outputFileInit->Close();
  unsigned numberOfThresholds = sizeof(listOfThresholds)/sizeof(listOfThresholds[0]);
  for (unsigned thresholdCounter = 0; thresholdCounter < numberOfThresholds; ++thresholdCounter) {
    TFile *outputFile = TFile::Open(outputDir+outputFileNamePrefix+Form(".root"), "UPDATE");
    Double_t threshold = listOfThresholds[thresholdCounter];
    TTree *outputTree = new TTree(Form("occupiedCellsTree_threshold_%.1f", threshold), Form("occupiedCellsTree_threshold_%.1f", threshold));
    createBranches(outputTree, etaBin, layer, occupiedCells);
    readDataIntoTree(outputTree, etaBin, layer, occupiedCells, threshold);
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
    // delete outputTree;
  }
  // outputFile->Close();
}

# ifndef __CINT__
int main() {
  populateTreesFromData();
  return 0;
}
# endif
