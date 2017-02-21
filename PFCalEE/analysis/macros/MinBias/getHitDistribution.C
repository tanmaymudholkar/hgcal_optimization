#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include<ctime>
#include<cmath>
#include<boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSDetector.hh"
// #include "HGCSSGeometryConversion.hh"

// const TString dataDir = Form("/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");

const bool useRaw = true;

const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");
const TString plotsDir = Form("plots/");
const TString outputDir = Form("/afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias/version_33/PythiaTest/");
const std::string inputLayerZPositionsFileName = "layerZPositions.dat";
const TString outputFileName = Form("rawHitDistribution_limitedStats_timing.root");
unsigned layerBoundaries[] = {0,28,40};
// const std::string granulStr = "0-27:4,28-39:4,40-51:4";

const unsigned hardcodedMaxDataCounter = 10;

HGCSSInfo *info;
double calorSizeXY = 0;
double cellSize = 0;
unsigned versionNumber = 0;
unsigned model = 0;
unsigned nLayers;

HGCSSDetector & myDetector = theDetector();

Double_t *layerZPositions;

const std::string granulStr = "0-27:4,28-39:4,40-51:4";

template <class T>
void extractParameterFromStr(std::string aStr,T & vec){ 
  if (aStr == "") return;
  std::vector<std::string> layVec;
  boost::split( layVec, aStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for parameter given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      exit(1);
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      exit(1);
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      if (iL < vec.size())
	std::istringstream(lPair[1])>>vec[iL];
      else {
	std::cout << " -- WARNING! Input parameter has more layers: " << endIdx << " than detector : " 
		  << vec.size()
		  << ". Ignoring additional layer #" << iL << "... PLEASE CHECK SETTINGS ARE CORRECT FOR EXISTING LAYERS!!"
		  << std::endl;
      }
    }
  }//loop on elements
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


bool testFile(TString inputPath) {
  TFile *testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    return false;
  }
  return true;
}

void readInputFiles(TChain *inputTree) {
  TString namePrefix = Form("DiginoXTalk_Test_noNoise_ICoff_Pu5_IC3_version33_000");
  if (useRaw) namePrefix = Form("HGcal_version33_000");
  TString suffix = Form(".root");
  unsigned dataCounter = 0;
  unsigned consecutiveNonexistentFiles = 0;
  while (consecutiveNonexistentFiles < 4 && (hardcodedMaxDataCounter == 0 || dataCounter < hardcodedMaxDataCounter)) {
    // while (consecutiveNonexistentFiles < 4) {
    std::cout << "Reading input files for data counter " << dataCounter << std::endl;
    // TFile *inputFile;
    TString fileName = dataDir + namePrefix + Form("%03d", dataCounter) + suffix;
    bool dataExists = testFile(fileName);
    if (dataExists) {
      inputTree->AddFile(fileName);
      consecutiveNonexistentFiles = 0;
      if (dataCounter == 0) {
        TFile *inputFile = TFile::Open(fileName);
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
    }
    else {
      ++consecutiveNonexistentFiles;
    }
    ++dataCounter;
  }
}

void printETA(Double_t timeRemainingInSec, unsigned ievt, unsigned nEvts_total) {
  Int_t hoursRemaining = static_cast<Int_t>(timeRemainingInSec/3600.);
  Int_t minutesRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining)/60.);
  Int_t secondsRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining - 60*minutesRemaining) + 0.5);
  
  std::cout << Form("\rApproximately %02d h: %02d m: %02d s remaining. Analyzing %04d/%04d", hoursRemaining, minutesRemaining, secondsRemaining, ievt, nEvts_total) << std::flush;
}

// template <typename hitType>
// void fillTree(std::vector<hitType> *hitVector, TTree *outputTree, Double_t &energy, unsigned &layer) {
//   for (unsigned iH(0); iH<(*hitVector).size(); ++iH){ // loop over rechits
//     hitType lHit = (*hitVector)[iH];
//     energy = lHit.energy();
//     if (energy > 0) {
//       layer = lHit.layer();
//       outputTree->Fill();
//     }
//   }
// }

inline Double_t getEta(Double_t xPosition, Double_t yPosition, Double_t zPosition) {
  return (ROOT::Math::XYZPoint(xPosition/10.,yPosition/10.,zPosition/10.)).eta();
}

void fillTreeRaw(std::vector<HGCSSSimHit> *simhitvec, TTree *outputTree, Double_t &energy, unsigned &layer, Double_t &xpos, Double_t &ypos, Double_t &eta, Double_t &radius, Double_t &timing, HGCSSGeometryConversion &geomConv) {
  unsigned prevLayer = 10000;
  bool isScint = false;

  for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){ // loop over simhits
    HGCSSSimHit lHit = (*simhitvec)[iH];
    energy = lHit.energy();
    if (energy > 0) {
      layer = lHit.layer();
      timing = lHit.time();
      if (layer != prevLayer) {
        const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
        isScint = subdet.isScint;
        prevLayer = layer;
      }
      std::pair<double,double> xy = lHit.get_xy(isScint,geomConv);
      xpos = xy.first;
      ypos = xy.second;
      radius = sqrt(xy.first*xy.first+xy.second*xy.second);
      Double_t zpos = layerZPositions[layer];
      eta = getEta(xpos, ypos, zpos);
      outputTree->Fill();
    }
  }
}

void fillTreeDigi(std::vector<HGCSSRecoHit> *rechitvec, TTree *outputTree, Double_t &energy, unsigned &layer, Double_t &xpos, Double_t &ypos, Double_t &eta, Double_t &radius, Double_t &timing, HGCSSGeometryConversion &geomConv) {
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
    HGCSSRecoHit lHit = (*rechitvec)[iH];
    energy = lHit.energy();
    if (energy > 0) {
      layer = lHit.layer();
      xpos = 0.;
      ypos = 0.;
      eta = 0.;
      radius = 0.;
      timing = 0.;
      outputTree->Fill();
    }
  }
  (void)geomConv;
}

void readDataIntoTree(TTree *outputTree, Double_t &energy, unsigned &layer, Double_t &xpos, Double_t &ypos, Double_t &eta, Double_t &radius, Double_t &timing) {
  TChain *lRecTree = new TChain("RecoTree");
  TChain *lSimTree = new TChain("HGCSSTree");

  if (useRaw) readInputFiles(lSimTree);
  else readInputFiles(lRecTree);

  myDetector.buildDetector(versionNumber,true,false,false);

  nLayers = myDetector.nLayers();

  HGCSSGeometryConversion geomConv(model,cellSize,false,2);
  geomConv.setXYwidth(calorSizeXY);
  geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  geomConv.initialiseSquareMap(calorSizeXY,10.);

  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1);
  extractParameterFromStr<std::vector<unsigned> >(granulStr,granularity);
  geomConv.setGranularity(granularity);

  fillLayerZPositions();

  std::vector<HGCSSSimHit> *simhitvec = 0;
  std::vector<HGCSSRecoHit> *rechitvec = 0;
  if (useRaw) lSimTree->SetBranchAddress("HGCSSSimHitVec", &simhitvec);
  else lRecTree->SetBranchAddress("HGCSSRecoHitVec", &rechitvec);
  unsigned nEvts_total = lRecTree->GetEntries();
  if (useRaw) nEvts_total = lSimTree->GetEntries();
  Double_t fractionCompleted;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  std::time_t initialTimeClocked = std::time(nullptr);
  if (useRaw) lSimTree->GetEntry(0);
  else lRecTree->GetEntry(0);
  if (useRaw) fillTreeRaw(simhitvec, outputTree, energy, layer, xpos, ypos, eta, radius, timing, geomConv);
  else fillTreeDigi(rechitvec, outputTree, energy, layer, xpos, ypos, eta, radius, timing, geomConv);

  for (unsigned ievt(1); ievt<nEvts_total; ++ievt){// loop on events
    
    fractionCompleted = float(ievt*1./nEvts_total);
    timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
    guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
    printETA(guess_timeRemainingInSec, 1+ievt, nEvts_total);
    
    if (useRaw) lSimTree->GetEntry(ievt);
    else lRecTree->GetEntry(ievt);
    if (useRaw) fillTreeRaw(simhitvec, outputTree, energy, layer, xpos, ypos, eta, radius, timing, geomConv);
    else fillTreeDigi(rechitvec, outputTree, energy, layer, xpos, ypos, eta, radius, timing, geomConv);
  }
  std::cout << std::endl;
}

void createBranches(TTree *outputTree, Double_t &energy, unsigned &layer, Double_t &xpos, Double_t &ypos, Double_t &eta, Double_t &radius, Double_t &timing) {
  outputTree->Branch("energy", &energy);
  outputTree->Branch("layer", &layer);
  outputTree->Branch("xpos", &xpos);
  outputTree->Branch("ypos", &ypos);
  outputTree->Branch("eta", &eta);
  outputTree->Branch("radius", &radius);
  outputTree->Branch("time", &timing);
}

void getHitDistribution() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  // Create addresses for variables to put in corresponding branches
  Double_t energy = 0;
  unsigned layer = 0;
  Double_t xpos = 0;
  Double_t ypos = 0;
  Double_t eta = 0;
  Double_t radius = 0;
  Double_t timing = 0;
  // Double_t zpos = 0;

  TFile *outputFile = new TFile(outputDir+outputFileName, "RECREATE");
  TTree *outputTree = new TTree("baseTree", "baseTree");
  createBranches(outputTree, energy, layer, xpos, ypos, eta, radius, timing);

  readDataIntoTree(outputTree, energy, layer, xpos, ypos, eta, radius, timing);

  outputFile->cd();
  outputTree->Write();
  delete outputTree;
  // outputFile->WriteTObject(outputTree);
  outputFile->Close();
  // delete outputFile;
  // delete outputTree;
}

# ifndef __CINT__
int main() {
  getHitDistribution();
  return 0;
}
# endif
