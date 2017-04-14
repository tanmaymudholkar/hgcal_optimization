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

#include "HGCSSRecoHit.hh"
#include "HGCSSHardcodedConstants.hh"

void printAxisBinning(TAxis *inputAxis) {
  unsigned nBins = inputAxis->GetNbins();
  for (unsigned binCounter = 0; binCounter <= 1+nBins; ++binCounter) {
    std::cout << "Bin number " << binCounter << " ranges from " << inputAxis->GetBinLowEdge(binCounter) << " to " << inputAxis->GetBinUpEdge(binCounter) << " and its center is " << inputAxis->GetBinCenter(binCounter) << std::endl;
  }
}

bool testFile(TString TInputPath) {
  TFile *testFile = TFile::Open(TInputPath);
  if ( !testFile ) {
    return false;
  }
  return true;
}

void readInputFiles(TChain *lRecTree, TString TDataDir, TString TInputFileNamePrefix, int hardcodedMaxDataCounter) {
  TString suffix = Form(".root");
  int dataCounter = 1;
  unsigned consecutiveNonexistentFiles = 0;
  bool continueReading = true;
  while (continueReading) {
    std::cout << "Reading input files for data counter " << dataCounter << std::endl;
    TString fileName = TDataDir + TInputFileNamePrefix + Form("%04d", dataCounter) + suffix;
    bool dataExists = testFile(fileName);
    if (dataExists) {
      lRecTree->AddFile(fileName);
      consecutiveNonexistentFiles = 0;
    }
    else {
      ++consecutiveNonexistentFiles;
    }
    ++dataCounter;
    continueReading = (consecutiveNonexistentFiles < 4);
    if (hardcodedMaxDataCounter != 0) continueReading = (continueReading && dataCounter <= hardcodedMaxDataCounter);
  }
}

void addEventToProfile(std::vector<HGCSSRecoHit> *rechitvec, TProfile2D *noiseProfile, double threshold, TAxis *axisWithRadialBinning) {
  // std::map<unsigned, std::map<int, Double_t> > noiseFracMapByEvent;
  // for (unsigned layerCounter = 0; layerCounter <= BHLASTLAYER; ++layerCounter) {
  //   for (int radialBinCounter = 0; radialBinCounter <= (1+nRadialBins); ++radialBinCounter) {
  //     noiseFracMapByEvent[layerCounter][radialBinCounter] = 0;
  //   }
  // }
  for (unsigned hitCounter(0); hitCounter<(*rechitvec).size(); ++hitCounter){ // loop over rechits
    Double_t energy = ((*rechitvec)[hitCounter]).energy();
    if (energy > threshold) {
      Double_t xpos = ((*rechitvec)[hitCounter]).get_x();
      Double_t ypos = ((*rechitvec)[hitCounter]).get_y();
      Double_t r = sqrt(xpos*xpos + ypos*ypos);
      unsigned radialBinForRecHit = axisWithRadialBinning->FindFixBin(r);
      unsigned layerForRecHit = ((*rechitvec)[hitCounter]).layer();
      double noiseFrac = ((*rechitvec)[hitCounter]).noiseFraction();
      // occupiedCellsMapByEvent[layerForRecHit][radialBinForRecHit] += 1;
      noiseProfile->Fill(double(layerForRecHit), double(radialBinForRecHit), noiseFrac);
    }
  }
  // for (unsigned layerCounter = 0; layerCounter <= BHLASTLAYER; ++layerCounter) {
  //   for (int radialBinCounter = 1; radialBinCounter <= (nRadialBins); ++radialBinCounter) { // ONLY RESTRICTED TO NON-OUTLIER BINS
  //     noiseProfile->Fill(double(layerCounter), double(radialBinCounter), double(occupiedCellsMapByEvent[layerCounter][radialBinCounter]));
  //   }
  // }
}

void printETA(Double_t timeRemainingInSec, unsigned ievt, unsigned nEvts_total) {
  Int_t hoursRemaining = static_cast<Int_t>(timeRemainingInSec/3600.);
  Int_t minutesRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining)/60.);
  Int_t secondsRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining - 60*minutesRemaining) + 0.5);
  std::cout << Form("\rApproximately %02d h: %02d m: %02d s remaining. Analyzing %04d/%04d", hoursRemaining, minutesRemaining, secondsRemaining, ievt, nEvts_total) << std::flush;
}

void fillNoiseProfileFromData(TProfile2D *noiseProfile, std::string dataDir, std::string inputFileNamePrefix, int hardcodedMaxDataCounter, double threshold, TAxis *axisWithRadialBinning) {
  std::cout << "Filling map of noise contribution to rechits from data..." << std::endl;
  TChain  *lRecTree = new TChain("RecoTree");
  TString TDataDir = Form(dataDir.c_str());
  TString TInputFileNamePrefix = Form(inputFileNamePrefix.c_str());
  readInputFiles(lRecTree, TDataDir, TInputFileNamePrefix, hardcodedMaxDataCounter);
  std::vector<HGCSSRecoHit> *rechitvec = 0;
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  unsigned nEvts_total = lRecTree->GetEntries();
  std::cout << "Total number of available events = " << nEvts_total << std::endl;
  Double_t fractionCompleted;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  std::time_t initialTimeClocked = std::time(nullptr);
  lRecTree->GetEntry(0);
  addEventToProfile(rechitvec, noiseProfile, threshold, axisWithRadialBinning);
  for (unsigned ievt(1); ievt<nEvts_total; ++ievt){// loop on events
    fractionCompleted = float(ievt*1./nEvts_total);
    timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
    guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
    printETA(guess_timeRemainingInSec, ievt, nEvts_total);
    
    lRecTree->GetEntry(ievt);
    addEventToProfile(rechitvec, noiseProfile, threshold, axisWithRadialBinning);
  }
  std::cout << std::endl;
  delete lRecTree;
}

void writeNoiseProfileToFiles(TProfile2D* noiseProfile, std::string outputDir, int nRadialBins, double threshold) {
  std::cout << "Writing noise profile to files..." << std::endl;
  // TString TThresholdStr = Form("/threshold_%.1f", threshold);
  std::string thresholdStr = std::string(Form("/threshold_%.1f", threshold));
  if (system(("mkdir -p " + outputDir + thresholdStr + " && rm -f " + outputDir + thresholdStr + "/*").c_str())) {
    std::cout << "Error, unable to create a new directory: " << outputDir << thresholdStr << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else {
    std::cout << "Output directory now fresh!" << std::endl;
  }
  for (unsigned layerCounter = 0; layerCounter <= BHLASTLAYER; ++layerCounter) {
    ofstream noiseFracOutputFile;
    ofstream noiseFracErrorsOutputFile;
    std::string noiseFracFileName = std::string(Form("/noiseFracInRadialBins_layer%i.dat", layerCounter));
    std::string noiseFracErrorsFileName = std::string(Form("/noiseFracErrorsInRadialBins_layer%i.dat", layerCounter));
    noiseFracOutputFile.open((outputDir + thresholdStr + noiseFracFileName).c_str());
    noiseFracErrorsOutputFile.open((outputDir + thresholdStr + noiseFracErrorsFileName).c_str());

    noiseFracOutputFile << 0 << "    " << 0 << std::endl; // underflow radial bin
    noiseFracErrorsOutputFile << 0 << "    " << 0 << std::endl; // underflow radial bin
    
    for (int radialBinCounter = 1; radialBinCounter <= nRadialBins; ++radialBinCounter) {
      Int_t globalBinNumber = noiseProfile->FindFixBin(double(layerCounter), double(radialBinCounter));
      double noise = noiseProfile->GetBinContent(globalBinNumber);
      double noiseError = noiseProfile->GetBinError(globalBinNumber);
      noiseFracOutputFile << radialBinCounter << "    " << noise << std::endl;
      noiseFracErrorsOutputFile << radialBinCounter << "    " << noiseError << std::endl;
    }

    noiseFracOutputFile << (1+nRadialBins) << "    " << 0 << std::endl; // overflow radial bin
    noiseFracErrorsOutputFile << (1+nRadialBins) << "    " << 0 << std::endl; // overflow radial bin
    
    noiseFracOutputFile.close();
    noiseFracErrorsOutputFile.close();
  }
}

void getNoiseContribToRechitsPassingThreshold(double minR, double maxR, int nRadialBins, std::string dataDir, std::string inputFileNamePrefix, std::string outputDir, double threshold, int hardcodedMaxDataCounter) {
  std::cout << "Beginning to populate files from data..." << std::endl;
  TAxis *axisWithRadialBinning = new TAxis(nRadialBins, minR, maxR);
  printAxisBinning(axisWithRadialBinning);
  TProfile2D *noiseProfile = new TProfile2D("noiseProfile", "Noise Profile;Layer;r(mm)", 1+BHLASTLAYER, -0.5, 0.5+BHLASTLAYER, nRadialBins, 0.5, 0.5+nRadialBins);
  fillNoiseProfileFromData(noiseProfile, dataDir, inputFileNamePrefix, hardcodedMaxDataCounter, threshold, axisWithRadialBinning);
  writeNoiseProfileToFiles(noiseProfile, outputDir, nRadialBins, threshold);
  delete noiseProfile;
}

# ifndef __CINT__
int main(int argc, char ** argv) {
  int nExpectedArguments = 8;
  if (argc != 1+nExpectedArguments) {
    std::cout << "Error: number of arguments expected: " << nExpectedArguments << "; number of arguments supplied: " << argc-1 << std::endl;
    std::exit(EXIT_FAILURE);
  }
  double minR = std::atof(argv[1]);
  double maxR = std::atof(argv[2]);
  int nRadialBins = std::atoi(argv[3]);
  std::string dataDir = std::string(argv[4]);
  std::string inputFileNamePrefix = std::string(argv[5]);
  // std::string totalNumberOfCellsInRadialBinsDir = std::string(argv[6]);
  std::string outputDir = std::string(argv[6]);
  // std::string outputFileNamePrefix = std::string(argv[8]);
  double threshold = std::atof(argv[7]);
  int hardcodedMaxDataCounter = std::atoi(argv[8]);
  std::cout << "Arguments supplied: " << std::endl
            << "minR = " << minR << std::endl
            << "maxR = " << maxR << std::endl
            << "nRadialBins = " << nRadialBins << std::endl
            << "dataDir = " << dataDir << std::endl
            << "inputFileNamePrefix = " << inputFileNamePrefix << std::endl
    // << "totalNumberOfCellsInRadialBinsDir = " << totalNumberOfCellsInRadialBinsDir << std::endl
            << "outputDir = " << outputDir << std::endl
    // << "outputFileNamePrefix = " << outputFileNamePrefix << std::endl
            << "threshold = " << threshold << std::endl
            << "hardcodedMaxDataCounter = " << hardcodedMaxDataCounter << std::endl
            << std::endl;
  getNoiseContribToRechitsPassingThreshold(minR, maxR, nRadialBins, dataDir, inputFileNamePrefix, outputDir, threshold, hardcodedMaxDataCounter);
  return 0;
}
# endif
