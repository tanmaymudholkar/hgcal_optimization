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
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"

// #include "HGCSSEvent.hh"
// #include "HGCSSInfo.hh"
#include "HGCSSRecoHit.hh"
// #include "HGCSSSimHit.hh"
// #include "HGCSSSamplingSection.hh"
// #include "HGCSSGenParticle.hh"
// #include "HGCSSDetector.hh"
// #include "HGCSSGeometryConversion.hh"

// const TString dataDir = Form("/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");

const TString dataDir = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");
const TString outputDir = Form("/afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias/version_33/PythiaTest/");
const TString outputFileName = Form("combinedTreeMinimal.root");

const unsigned hardcodedMaxDataCounter = 3;

bool testFile(TString inputPath) {
  TFile *testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    return false;
  }
  return true;
}

void readInputFiles(TChain *lRecTree) {
  TString digiPrefix = Form("DigiPu200_IC3_version33_000");
  TString suffix = Form(".root");
  unsigned dataCounter = 0;
  unsigned consecutiveNonexistentFiles = 0;
  while (consecutiveNonexistentFiles < 4 && (hardcodedMaxDataCounter == 0 || dataCounter < hardcodedMaxDataCounter)) {
    // while (consecutiveNonexistentFiles < 4) {
    std::cout << "Reading input files for data counter " << dataCounter << std::endl;
    // TFile *inputFile;
    TString fileName = dataDir + digiPrefix + Form("%03d", dataCounter) + suffix;
    bool dataExists = testFile(fileName);
    if (dataExists) {
      lRecTree->AddFile(fileName);
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

void fillTree(std::vector<HGCSSRecoHit> *rechitvec, TTree *outputTree, Double_t &energy, Double_t &eta, Double_t &xpos, Double_t &ypos, Int_t &layer) {
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
    energy = ((*rechitvec)[iH]).energy();
    eta = ((*rechitvec)[iH]).eta();
    xpos = ((*rechitvec)[iH]).get_x();
    ypos = ((*rechitvec)[iH]).get_y();
    layer = ((*rechitvec)[iH]).layer();
    // if (eta < 10) {
    outputTree->Fill();
      // if (layer > 42) std::cout << "layer number " << layer << std::endl;
    // }
    if (layer >=39 && layer <= 41) std::cout << "xpos = " << xpos << ", ypos = " << ypos << ", layer = " << layer << std::endl;
  }
}

void readDataIntoTree(TTree *outputTree, Double_t &energy, Double_t &eta, Double_t &xpos, Double_t &ypos, Int_t &layer) {
  TChain  *lRecTree = new TChain("RecoTree");
  readInputFiles(lRecTree);
  std::vector<HGCSSRecoHit> *rechitvec = 0;
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  unsigned nEvts_total = lRecTree->GetEntries();
  Double_t fractionCompleted;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  std::time_t initialTimeClocked = std::time(nullptr);
  for (unsigned ievt(0); ievt<nEvts_total; ++ievt){// loop on events
    if (ievt%10==0) {
      if (ievt == 0) {
        initialTimeClocked = std::time(nullptr);
      }
      else {
        fractionCompleted = float(ievt*1./nEvts_total);
        timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
        guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
        printETA(guess_timeRemainingInSec, ievt, nEvts_total);
      }
    }
    lRecTree->GetEntry(ievt);
    fillTree(rechitvec, outputTree, energy, eta, xpos, ypos, layer);
  }
  std::cout << std::endl;
}

void createBranches(TTree *outputTree, Double_t &energy, Double_t &eta, Double_t &xpos, Double_t &ypos, Int_t &layer) {
  TBranch *energyBranch = outputTree->Branch("energy", &energy);
  TBranch *etaBranch = outputTree->Branch("eta", &eta);
  TBranch *xposBranch = outputTree->Branch("xpos", &xpos);
  TBranch *yposBranch = outputTree->Branch("ypos", &ypos);
  TBranch *layerBranch = outputTree->Branch("layer", &layer);

  (void)(energyBranch);
  (void)(etaBranch);
  (void)(xposBranch);
  (void)(yposBranch);
  (void)(layerBranch);
}

void populateTree() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  // Create addresses for variables to put in corresponding branches
  Double_t energy = 0;
  Double_t eta = 0;
  Double_t xpos = 0;
  Double_t ypos = 0;
  Int_t layer = 0;

  TFile *outputFile = new TFile(outputDir+outputFileName, "RECREATE");

  TTree *outputTree = new TTree("baseTree", "baseTree");
  createBranches(outputTree, energy, eta, xpos, ypos, layer);

  readDataIntoTree(outputTree, energy, eta, xpos, ypos, layer);
  
  outputFile->WriteTObject(outputTree);
  outputFile->Close();
  // delete outputFile;
  // delete outputTree;
}

# ifndef __CINT__
int main() {
  populateTree();
  return 0;
}
# endif
