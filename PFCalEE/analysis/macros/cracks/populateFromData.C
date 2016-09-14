#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<stdlib.h>
#include<map>
#include<ctime>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TVector.h"
#include "TLine.h"

#include "../../../userlib/include/HGCSSEvent.hh"
#include "../../../userlib/include/HGCSSInfo.hh"
#include "../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../userlib/include/HGCSSSimHit.hh"
#include "../../../userlib/include/HGCSSSamplingSection.hh"
#include "../../../userlib/include/HGCSSGenParticle.hh"

const Int_t version_number = 100;
// const TString datadir = Form("/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/cmst3/group/hgcal/HGCalCracks/gitV06e-04-06/gamma");
const TString datadir = Form("/afs/cern.ch/user/t/tmudholk/private/temp/testcracks");
const TString pathToLayerWeights = Form("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/analysis/macros/cracks/layer_weights_version100_averaged.dat");
const TString et_portion = Form("et%.0f", 60.0);
const TString eta_portion = Form("_eta%.3f", 2.0);
const TString version_name = Form("cracks_study_2mmgap_averaged_weights");
const Double_t targetBinSize_mm = 2.;
const Double_t particleXMax_mm = 290.;
const Double_t particleXMin_mm = -170.;
const Double_t eta = 2.0;
const Double_t sigmas_up = 6.0;
const Double_t sigmas_down = 6.0;
const Int_t numberOfBinsInEnergyHistograms = 50;
Double_t idealNBins;
Int_t nBinsRequired;
Double_t lowerEdgeFirstBin;
Double_t upperEdgeLastBin;
std::vector<Double_t> weights;
// std::vector<Double_t> binCenters;

void setBinning() {
  idealNBins = (particleXMax_mm-particleXMin_mm)/targetBinSize_mm;
  lowerEdgeFirstBin = particleXMin_mm;
  nBinsRequired = 1+static_cast<Int_t>(idealNBins);
  upperEdgeLastBin = particleXMin_mm + targetBinSize_mm*(nBinsRequired);
  if (static_cast<Double_t>(static_cast<Int_t>(idealNBins)) == idealNBins) {
    nBinsRequired = static_cast<Int_t>(idealNBins);
    upperEdgeLastBin = particleXMin_mm + targetBinSize_mm*(nBinsRequired);
  }
  // Double_t runningCenter = lowerEdgeFirstBin + 0.5*targetBinSize_mm;
  // for (Int_t binCounter = 0; binCounter < nBinsRequired; ++binCounter) {
  //   binCenters.push_back(runningCenter);
  //   runningCenter += targetBinSize_mm;
  // }
}

void readWeights() {
  std::cout << "Reading weights..." << std::endl;
  std::vector<Double_t> weights_raw;
  ifstream f_layer_weights_raw;
  f_layer_weights_raw.open(pathToLayerWeights);
  std::string line;
  Double_t weight;
  if (f_layer_weights_raw.is_open()) {
    while (getline(f_layer_weights_raw,line)) {
      weight = std::atof(line.c_str());
      weights_raw.push_back(weight);
    }
  }
  f_layer_weights_raw.close();
  if (weights_raw.empty()) {
    std::cout << "Unable to read weights from file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for (unsigned layer_counter = 0; layer_counter != (weights_raw.size()); layer_counter++) {
    weights.push_back(weights_raw[layer_counter]);
  }
  std::cout << std::endl;
  std::cout << "Weights:" << std::endl;
  for (unsigned layer_counter = 0; layer_counter != weights.size(); layer_counter++) {
    std::cout << "layer " << layer_counter << ": " << weights[layer_counter] << std::endl;
  }
}

void initializeEnergyHistogramsMap(std::map<Int_t, TH1F*> &energyHistogramsMap, std::map<Int_t, std::pair<Double_t, Double_t> > histogramRangesMap) {
  for(Int_t binCounter=0; binCounter<=(1+nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    energyHistogramsMap[binCounter] = new TH1F(Form("h_energyHistogram_binNumber_%i",binCounter), Form("number of hits; energy deposited"), numberOfBinsInEnergyHistograms, histogramRangesMap[binCounter].first, histogramRangesMap[binCounter].second);
  }
}

void writeEnergyHistogramsToFile(std::map<Int_t, TH1F*> &energyHistogramsMap, TFile *histogramsOutputFile) {
  for (Int_t binCounter = 0; binCounter <= (1+nBinsRequired); ++binCounter) {
    histogramsOutputFile->WriteTObject(energyHistogramsMap[binCounter]);
  }
}

void cleanEnergyHistogramsMap(std::map<Int_t, TH1F*> &energyHistogramsMap) {
  for(Int_t binCounter=0; binCounter<=(1+nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    delete energyHistogramsMap[binCounter];
  }
}

bool testInputFile(TString inputPath, TFile* testFile) {
  testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    // std::cout << " -- Error, input file " << inputPath << " cannot be opened. Skipping..." << std::endl;
    return false;
  } // else std::cout << " -- input file " << testFile->GetName() << " successfully opened." << std::endl;
  delete testFile;
  return true;
}

void printHistogramData(TH1F *inputHistogram) {
  std::cout << "Printing data for histogram with name: " << inputHistogram->GetName() << std::endl;
  TAxis *h_Xaxis = inputHistogram->GetXaxis();
  Int_t h_Nbins = h_Xaxis->GetNbins();
  std::cout << "Underflow:    " << inputHistogram->GetBinContent(0) << std::endl;
  for(Int_t binCounter=1; binCounter<=h_Nbins; ++binCounter) {
    std::cout << h_Xaxis->GetBinCenter(binCounter) << ":    " << inputHistogram->GetBinContent(binCounter) << std::endl;
  }
  std::cout << "Overflow:    " << inputHistogram->GetBinContent(h_Nbins+1) << std::endl;
}

void read_input_files(TChain *lSimTree, TChain *lRecTree) {

  TString Digi_common_prefix, HGcal_common_prefix;
  TString common_suffix = Form(".root");

  Digi_common_prefix = datadir + Form("/DigiIC3_Si2__version%i_model4_BOFF_",version_number);
  HGcal_common_prefix = datadir + Form("/HGcal__version%i_model4_BOFF_",version_number);
    
  unsigned run_no = 1;
  Int_t consecutive_nonexistent_files = 0;

  while (consecutive_nonexistent_files < 4) {
    std::cout << "Reading input files for run number " << run_no << std::endl;
      
    TFile *testFile_hgcal(0);
    TFile *testFile_digi(0);

    bool hgcal_data_exists, digi_data_exists;
    
    hgcal_data_exists=testInputFile(HGcal_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix, testFile_hgcal);
    digi_data_exists=testInputFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix, testFile_digi);
    
    
    if (digi_data_exists && hgcal_data_exists) {
      consecutive_nonexistent_files = 0;
      lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
      lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
    }
    else {
      consecutive_nonexistent_files++;
      if (!(hgcal_data_exists)) {
        std::cout << "HGcal data does not exist for run no " << run_no << std::endl;
      }
      if (!(digi_data_exists)) {
        std::cout << "Digi data does not exist for run no " << run_no << std::endl;
      }
    }
    delete testFile_digi;
    delete testFile_hgcal;
    run_no++;
  }// ends loop over runs
}

void analyzeEvent(std::vector<HGCSSRecoHit> *rechitvec, std::vector<Double_t> &energiesPerEventInBin, Double_t &sumEnergyInBin, Double_t &sumEnergySquareInBin, Int_t &numberOfEventsInBin) {
  Double_t eventEnergy = 0;
  for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
    unsigned hitLayer = ((*rechitvec)[iH]).layer();
    Double_t hitEnergy = ((*rechitvec)[iH]).energy();
    Double_t hitEnergy_weighted = hitEnergy*weights[hitLayer]/tanh(eta);
    eventEnergy += hitEnergy_weighted;
  }// end loop over rechits
  energiesPerEventInBin.push_back(eventEnergy);
  sumEnergyInBin += eventEnergy;
  sumEnergySquareInBin += eventEnergy*eventEnergy;
  numberOfEventsInBin += 1;
}

void initializeEventAnalyzerMaps(std::map<Int_t, Double_t> &sumEnergy, std::map<Int_t, Double_t> &sumEnergySquare, std::map<Int_t, Int_t> &numberOfEvents) {
  for(Int_t binCounter=0; binCounter<=(1+nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    // std::vector<Double_t> energiesPerEvent[binCounter];
    sumEnergy[binCounter] = 0;
    sumEnergySquare[binCounter] = 0;
    numberOfEvents[binCounter] = 0;
  }
}

void populateHistogramRangesMap(std::map<Int_t, std::pair<Double_t, Double_t> > &histogramRangesMap, std::map<Int_t, Double_t> sumEnergy, std::map<Int_t, Double_t> sumEnergySquare, std::map<Int_t, Int_t> numberOfEvents) {
  for(Int_t binCounter=0; binCounter<=(1+nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    Double_t sumEnergyInBin = sumEnergy[binCounter];
    Double_t sumEnergySquareInBin = sumEnergySquare[binCounter];
    Int_t numberOfEventsInBin = numberOfEvents[binCounter];

    Double_t statisticalAverageEnergy = sumEnergyInBin/numberOfEventsInBin;
    Double_t statisticalRMSEnergy = sqrt((sumEnergySquareInBin/numberOfEventsInBin)-statisticalAverageEnergy*statisticalAverageEnergy);
    (histogramRangesMap[binCounter]).first = statisticalAverageEnergy - sigmas_down*statisticalRMSEnergy;
    (histogramRangesMap[binCounter]).second = statisticalAverageEnergy + sigmas_up*statisticalRMSEnergy;
  }
}

void fillEnergyHistograms(std::map<Int_t, TH1F*> &energyHistogramsMap, std::map<Int_t, std::vector<Double_t> > energiesPerEvent) {
  for(Int_t binCounter=0; binCounter<=(1+nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    if ((energiesPerEvent[binCounter]).size() > 0) {
      for (std::vector<Double_t>::iterator eventEnergyIterator = (energiesPerEvent[binCounter]).begin(); eventEnergyIterator != (energiesPerEvent[binCounter]).end(); ++eventEnergyIterator) {
        (energyHistogramsMap[binCounter])->Fill(*eventEnergyIterator);
      }
    }
  }
}

void printETA(Double_t timeRemainingInSec) {
  Double_t timeRemainingInSecLeftover = timeRemainingInSec;
  Int_t hoursRemaining = static_cast<Int_t>(timeRemainingInSecLeftover/3600.);
  timeRemainingInSecLeftover = fmod(timeRemainingInSecLeftover, 3600.);
  Int_t minutesRemaining = static_cast<Int_t>(timeRemainingInSecLeftover/60.);
  timeRemainingInSecLeftover = fmod(timeRemainingInSecLeftover, 60.);
  std::cout << Form("\rApproximately %02d h: %02d m: %02.2f s remaining", hoursRemaining, minutesRemaining, timeRemainingInSecLeftover) << std::flush;
}

void populateEnergyHistogramsMap() {
  TChain  *lSimTree = new TChain("HGCSSTree");
  TChain  *lRecTree = new TChain("RecoTree");
  read_input_files (lSimTree, lRecTree);

  std::vector<HGCSSRecoHit> *rechitvec = 0;
  std::vector<HGCSSSimHit> *simhitvec = 0;
  // std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSGenParticle> *genvec = 0; // get the generated particle vector
  HGCSSEvent * event = 0; // no problem: clearly the header files are being read and the class HGCSSEvent can be used

  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  // lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lSimTree->SetBranchAddress("HGCSSEvent",&event); // Problem ??!!
    
  unsigned nEvts_total = lSimTree->GetEntries();
  std::cout << "number of events: " << nEvts_total << std::endl;
  
  TH1F *h_generatedParticlesX = new TH1F(Form("h_generatedParticlesX"), Form("number of events; x"), nBinsRequired, lowerEdgeFirstBin, upperEdgeLastBin);
  TH1F *h_generatedParticlesY = new TH1F(Form("h_generatedParticlesY"), Form("number of events; y"), nBinsRequired, 0, 0);
  // TH1F *h_generatedParticlesZ = new TH1F(Form("h_generatedParticlesZ"), Form("number of events; z"), nBinsRequired, 0, 0);

  TAxis *axisWithCorrectBinning = h_generatedParticlesX->GetXaxis();

  std::map<Int_t, std::vector<Double_t> > energiesPerEvent;
  std::map<Int_t, Double_t> sumEnergy;
  std::map<Int_t, Double_t> sumEnergySquare;
  std::map<Int_t, Int_t> numberOfEvents;

  initializeEventAnalyzerMaps(sumEnergy, sumEnergySquare, numberOfEvents);

  Double_t fractionCompleted;
  // Double_t timeElapsedInSec;
  Int_t timeElapsedInSec;
  Double_t guess_timeRemainingInSec;
  // clock_t initialCyclesClocked;
  std::time_t initialTimeClocked = std::time(nullptr);

  for (unsigned ievt(0); ievt<nEvts_total; ++ievt){// loop on events
    if (ievt%100==0) {
      if (ievt == 0) {
        // initialCyclesClocked = clock();
        initialTimeClocked = std::time(nullptr);
      }
      else {
        fractionCompleted = float(ievt*1./nEvts_total);
        // timeElapsedInSec = float((clock() - initialCyclesClocked)/CLOCKS_PER_SEC);
        timeElapsedInSec = int(std::time(nullptr)-initialTimeClocked);
        guess_timeRemainingInSec = timeElapsedInSec*(1.-fractionCompleted)/fractionCompleted;
        printETA(guess_timeRemainingInSec);
      }
    }
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);
    Double_t posX_gen_ini = ((*genvec)[0]).x();
    Double_t posY_gen_ini = ((*genvec)[0]).y();
    // Double_t posX_gen_ini = event.vtx_x();
    // Double_t posY_gen_ini = event.vtx_y();
    // Double_t posZ_gen_ini = event.vtx_z();
    h_generatedParticlesX->Fill(posX_gen_ini);
    Double_t binNumber = axisWithCorrectBinning->FindBin(posX_gen_ini);
    // std::cout << "Checking... for initial position: " << posX_gen_ini << ", bin number: " << binNumber << ", corresponding bin position: " << axisWithCorrectBinning->GetBinCenter(binNumber) << std::endl;
    analyzeEvent(rechitvec, energiesPerEvent[binNumber], sumEnergy[binNumber], sumEnergySquare[binNumber], numberOfEvents[binNumber]);
    h_generatedParticlesY->Fill(posY_gen_ini);
    // h_generatedParticlesZ->Fill(posZ_gen_ini);
  }
  std::cout << std::endl;
  TFile *histogramsOutputFile = new TFile(Form("root_histograms/histograms_")+version_name+eta_portion+Form("_results.root"),"RECREATE");
  histogramsOutputFile->WriteTObject(h_generatedParticlesX);
  std::cout << "Printing histograms data for histogram of x-positions of generated particles:" << std::endl;
  printHistogramData(h_generatedParticlesX);
  histogramsOutputFile->WriteTObject(h_generatedParticlesY);
  // histogramsOutputFile->WriteTObject(h_generatedParticlesZ);
  delete h_generatedParticlesX;
  delete h_generatedParticlesY;
  // delete h_generatedParticlesZ;

  std::map<Int_t, std::pair<Double_t, Double_t> > histogramRangesMap;
  populateHistogramRangesMap(histogramRangesMap, sumEnergy, sumEnergySquare, numberOfEvents);

  std::map<Int_t, TH1F*> energyHistogramsMap;
  initializeEnergyHistogramsMap(energyHistogramsMap, histogramRangesMap);
  fillEnergyHistograms(energyHistogramsMap, energiesPerEvent);
  writeEnergyHistogramsToFile(energyHistogramsMap, histogramsOutputFile);
  cleanEnergyHistogramsMap(energyHistogramsMap);
  delete histogramsOutputFile;
}

void populateFromData() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");
  setBinning();
  readWeights();
  // std::cout << "Checking... bin centers: " << std::endl;
  // for (std::vector<Double_t>::iterator binCenters_iterator = binCenters.begin(); binCenters_iterator != binCenters.end(); ++binCenters_iterator) {
  //   std::cout << *binCenters_iterator << std::endl;
  // }
  populateEnergyHistogramsMap();
  std::cout << "Finished populating histograms from data!" << std::endl;
}
