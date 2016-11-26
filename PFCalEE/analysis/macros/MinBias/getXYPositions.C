#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include<ctime>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
// #include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSGenParticle.hh"
#include "HGCSSDetector.hh"
#include "HGCSSGeometryConversion.hh"

// const TString dataDir = Form("/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");

const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");
const TString outputDir = Form("/afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias/version_33/PythiaTest/");
const TString outputFileName = Form("xyPositionsRaw_minimal.root");
const std::string granulStr = "0-27:4,28-39:4,40-51:4";

const unsigned hardcodedMaxDataCounter = 2;
const unsigned hardcodedMaxEvents = 100;

HGCSSInfo *info;
double calorSizeXY = 0;
double cellSize = 0;
unsigned versionNumber = 0;
unsigned model = 0;

HGCSSDetector & myDetector = theDetector();

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

bool testFile(TString inputPath) {
  TFile *testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    return false;
  }
  return true;
}

void readInputFiles(TChain *lSimTree) {
  TString hgcalPrefix = Form("HGcal_version33_000");
  TString suffix = Form(".root");
 
  std::cout << "Reading first input data file" << std::endl;
  TString fileName = dataDir + hgcalPrefix + Form("000") + suffix;
  bool dataExists = testFile(fileName);

  if (dataExists) {
    lSimTree->AddFile(fileName);
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
  else {
    std::cout << "Input files not found" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  unsigned dataCounter = 1;
  unsigned consecutiveNonexistentFiles = 0;
  while (consecutiveNonexistentFiles < 4 && (hardcodedMaxDataCounter == 0 || dataCounter < hardcodedMaxDataCounter)) {
    // while (consecutiveNonexistentFiles < 4) {
    std::cout << "Reading input files for data counter " << dataCounter << std::endl;
    // TFile *inputFile;
    fileName = dataDir + hgcalPrefix + Form("%03d", dataCounter) + suffix;
    dataExists = testFile(fileName);
    if (dataExists) {
      lSimTree->AddFile(fileName);
      consecutiveNonexistentFiles = 0;
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

void fillTree(std::vector<HGCSSSimHit> *simhitvec, TTree *outputTree, Double_t &energy, Double_t &xpos, Double_t &ypos, unsigned &layer, unsigned &cellid, HGCSSGeometryConversion &geomConv) {
  unsigned prevLayer = 10000;
  bool isScint = false;

  for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){ // loop over simhits
    HGCSSSimHit lHit = (*simhitvec)[iH];
    energy = lHit.energy();
    if (energy > 0) {
      layer = lHit.layer();
      cellid = lHit.cellid();
      if (layer != prevLayer) {
        const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
        isScint = subdet.isScint;
        prevLayer = layer;
      }
      
      std::pair<double,double> xy = lHit.get_xy(isScint,geomConv);
      xpos = xy.first;
      ypos = xy.second;
          
      // xpos = 0.;
      // ypos = 0.;
    
      // xpos = lHit.get_x();
      // ypos = lHit.get_y();
    
      // if (eta < 10) {
      outputTree->Fill();
      // if (layer > 42) std::cout << "layer number " << layer << std::endl;
      // }
      // if (layer >=39 && layer <= 41) {
      //   std::cout << "energy = " << energy << ", xpos = " << xpos << ", ypos = " << ypos << ", layer = " << layer << ", cellid = " << cellid <<  std::endl;
      //   std::cout << "isScint: " << (isScint? "true" : "false") << std::endl;
      // }
    }
  }
}

void readDataIntoTree(TTree *outputTree, Double_t &energy, Double_t &xpos, Double_t &ypos, unsigned &layer, unsigned &cellid) {

  TChain  *lSimTree = new TChain("HGCSSTree");
  readInputFiles(lSimTree);

  // HGCSSDetector & myDetector = theDetector();
  myDetector.buildDetector(versionNumber,true,false,false);

  HGCSSGeometryConversion geomConv(model,cellSize,false,2);
  geomConv.setXYwidth(calorSizeXY);
  geomConv.initialiseHoneyComb(calorSizeXY,cellSize);
  geomConv.initialiseSquareMap(calorSizeXY/2.,10.);

  const unsigned nLayers = myDetector.nLayers();
  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1);
  extractParameterFromStr<std::vector<unsigned> >(granulStr,granularity);

  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "granularity at " << iL << " : " << granularity[iL] << std::endl;
  }
  
  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  std::vector<HGCSSSimHit> *simhitvec = 0;
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  unsigned nEvts_total = lSimTree->GetEntries();
  Double_t fractionCompleted;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  std::time_t initialTimeClocked = std::time(nullptr);
  lSimTree->GetEntry(0);
  fillTree(simhitvec, outputTree, energy, xpos, ypos, layer, cellid, geomConv);

  for (unsigned ievt(1); ievt<nEvts_total && (hardcodedMaxEvents == 0 || ievt < hardcodedMaxEvents); ++ievt){// loop on events
    
    fractionCompleted = float(ievt*1./nEvts_total);
    timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
    guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
    printETA(guess_timeRemainingInSec, 1+ievt, nEvts_total);
    
    lSimTree->GetEntry(ievt);
    fillTree(simhitvec, outputTree, energy, xpos, ypos, layer, cellid, geomConv);
  }
  std::cout << std::endl;
}

void createBranches(TTree *outputTree, Double_t &energy, Double_t &xpos, Double_t &ypos, unsigned &layer, unsigned &cellid) {
  TBranch *energyBranch = outputTree->Branch("energy", &energy);
  // TBranch *etaBranch = outputTree->Branch("eta", &eta);
  TBranch *xposBranch = outputTree->Branch("xpos", &xpos);
  TBranch *yposBranch = outputTree->Branch("ypos", &ypos);
  TBranch *layerBranch = outputTree->Branch("layer", &layer);
  TBranch *cellidBranch = outputTree->Branch("cellid", &cellid);

  (void)(energyBranch);
  // (void)(etaBranch);
  (void)(xposBranch);
  (void)(yposBranch);
  (void)(layerBranch);
  (void)(cellidBranch);
}

void getXYPositions() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  // Create addresses for variables to put in corresponding branches
  Double_t energy = 0;
  // Double_t eta = 0;
  Double_t xpos = 0;
  Double_t ypos = 0;
  unsigned layer = 0;
  unsigned cellid = 0;

  TFile *outputFile = new TFile(outputDir+outputFileName, "RECREATE");

  TTree *outputTree = new TTree("baseTree", "baseTree");
  createBranches(outputTree, energy, xpos, ypos, layer, cellid);

  readDataIntoTree(outputTree, energy, xpos, ypos, layer, cellid);
  
  outputFile->WriteTObject(outputTree);
  outputFile->Close();
  // delete outputFile;
  // delete outputTree;
}

# ifndef __CINT__
int main() {
  getXYPositions();
  return 0;
}
# endif
