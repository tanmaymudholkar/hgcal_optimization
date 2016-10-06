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
#include "TFitResult.h"
#include "TAxis.h"

// #include "../../../userlib/include/HGCSSEvent.hh"
// #include "../../../userlib/include/HGCSSInfo.hh"
// #include "../../../userlib/include/HGCSSRecoHit.hh"
// #include "../../../userlib/include/HGCSSSimHit.hh"
// #include "../../../userlib/include/HGCSSSamplingSection.hh"
// #include "../../../userlib/include/HGCSSGenParticle.hh"

#include "HGCSSEvent.hh"
#include "HGCSSInfo.hh"
#include "HGCSSRecoHit.hh"
#include "HGCSSSimHit.hh"
#include "HGCSSSamplingSection.hh"
#include "HGCSSGenParticle.hh"


const Int_t version_number = 100;
const TString datadir = Form("/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/cmst3/group/hgcal/HGCalCracks/gitV06e-04-06/gamma");
// const TString datadir = Form("/afs/cern.ch/user/t/tmudholk/private/temp/testcracks");
const TString pathToLayerWeights = Form("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/analysis/macros/cracks/layer_weights_version100.dat");
const TString et_portion = Form("et%.0f", 60.0);
const Double_t incomingParticleEnergy = 225.731;
const TString eta_portion = Form("_eta%.3f", 2.0);
const Double_t particleXMax_mm = 290.;
const Double_t particleXMin_mm = -170.;
const Double_t crackWidth_mm = 6.;
const Double_t crack1Position_mm = 224.;
const Double_t crack2Position_mm = 254.;
const Double_t targetBinSizeAwayFromCracks_mm = 9.;
const Double_t targetBinSizeCloseToCracks_mm = 3.;
const Double_t targetBinSizeInsideCracks_mm = 2.;
const Double_t targetBinSizeBetweenCracks_mm = 3.;
const Double_t eta = 2.0;
const Double_t sigmas_up = 6.0;
const Double_t sigmas_down = 6.0;
const Int_t numberOfBinsInEnergyHistograms = 50;

// TString m_versionName;
TString m_versionName = Form("cracks_study_with_gev");
// Double_t m_targetBinSize_mm;

Int_t m_nBinsRequired;
// Double_t lowerEdgeFirstBin = 0.;
// Double_t upperEdgeLastBin = 0.;

Double_t m_baseEnergy = 0;
Double_t m_averageWidth = 0;
Double_t m_mipsToGeVConversion = 0;

Double_t * m_binEdges;

std::vector<Double_t> weights;
// std::vector<Double_t> binCenters;

std::pair<Int_t,Double_t> getRequiredBinningParameters(Double_t length_mm, Double_t targetBinSize_mm) {
  std::pair<Int_t,Double_t> requiredBinningParameters;
  requiredBinningParameters.first = static_cast<Int_t>(length_mm/targetBinSize_mm + 0.5);
  requiredBinningParameters.second = length_mm/requiredBinningParameters.first;
  return requiredBinningParameters;
}

void buildBinning(std::pair<Int_t,Double_t>binningParameters, Int_t &runningIndex, Double_t startingFromPosition) {
  Double_t runningPosition = startingFromPosition;
  Int_t maxRunningIndex = runningIndex + binningParameters.first;
  for (Int_t regionCounter = runningIndex; regionCounter < maxRunningIndex; ++regionCounter) {
    m_binEdges[regionCounter] = runningPosition;
    runningPosition += binningParameters.second;
    ++runningIndex;
  }
}

void setBinning() {
  // Double_t idealNBins = (particleXMax_mm-particleXMin_mm)/m_targetBinSize_mm;
  // lowerEdgeFirstBin = particleXMin_mm;
  // m_nBinsRequired = 1+static_cast<Int_t>(idealNBins);
  // upperEdgeLastBin = particleXMin_mm + m_targetBinSize_mm*(m_nBinsRequired);
  // if (static_cast<Double_t>(static_cast<Int_t>(idealNBins)) == idealNBins) {
  //   m_nBinsRequired = static_cast<Int_t>(idealNBins);
  //   upperEdgeLastBin = particleXMin_mm + m_targetBinSize_mm*(m_nBinsRequired);
  // }
  // // Double_t runningCenter = lowerEdgeFirstBin + 0.5*m_targetBinSize_mm;
  // // for (Int_t binCounter = 0; binCounter < m_nBinsRequired; ++binCounter) {
  // //   binCenters.push_back(runningCenter);
  // //   runningCenter += m_targetBinSize_mm;
  // // }

  // Regions:
  // Away from cracks -- CloseToCracks -- Crack -- Between Cracks -- Crack -- CloseToCracks

  // Double_t binSizeBetweenCracks_mm = (crack2Position_mm-(crack1Position_mm + crackWidth_mm))/nBinsBetweenCracks;
    
  Double_t distanceOfMaxPositionFromCrack2_mm = particleXMax_mm - (crack2Position_mm + crackWidth_mm);
  Double_t positionFirstBinningChange_mm = crack1Position_mm - distanceOfMaxPositionFromCrack2_mm;

  std::pair<Int_t,Double_t> binningParametersAwayFromCracks = getRequiredBinningParameters(positionFirstBinningChange_mm - particleXMin_mm, targetBinSizeAwayFromCracks_mm);

  // Int_t nBinsAwayFromCracks = static_cast<Int_t>((positionFirstBinningChange_mm - particleXMin_mm)/targetBinSizeAwayFromCracks_mm + 0.5);
  // Double_t binSizeAwayFromCracks = (positionFirstBinningChange_mm - particleXMin_mm)/nBinsAwayFromCracks;

  std::pair<Int_t,Double_t> binningParametersCloseToCracks = getRequiredBinningParameters(distanceOfMaxPositionFromCrack2_mm, targetBinSizeCloseToCracks_mm);

  std::pair<Int_t,Double_t> binningParametersInsideCracks = getRequiredBinningParameters(crackWidth_mm, targetBinSizeInsideCracks_mm);

  std::pair<Int_t,Double_t> binningParametersBetweenCracks = getRequiredBinningParameters(crack2Position_mm - (crack1Position_mm + crackWidth_mm), targetBinSizeBetweenCracks_mm);

  m_nBinsRequired = binningParametersAwayFromCracks.first + binningParametersCloseToCracks.first + binningParametersInsideCracks.first + binningParametersBetweenCracks.first + binningParametersInsideCracks.first + binningParametersCloseToCracks.first;

  m_binEdges = new Double_t[1+m_nBinsRequired];

  Int_t runningIndex = 0;

  buildBinning(binningParametersAwayFromCracks, runningIndex, particleXMin_mm);
  buildBinning(binningParametersCloseToCracks, runningIndex, positionFirstBinningChange_mm);
  buildBinning(binningParametersInsideCracks, runningIndex, crack1Position_mm);
  buildBinning(binningParametersBetweenCracks, runningIndex, crack1Position_mm + crackWidth_mm);
  buildBinning(binningParametersInsideCracks, runningIndex, crack2Position_mm);
  buildBinning(binningParametersCloseToCracks, runningIndex, crack2Position_mm + crackWidth_mm);

  m_binEdges[runningIndex] = particleXMax_mm; // <-- Right edge of last bin
  ++runningIndex;

  if (runningIndex != 1+m_nBinsRequired) {
    std::cout << "Something has gone horribly wrong in the binning! At the end of building, runningIndex is " << runningIndex << " while 1+m_nBinsRequired is " << 1+m_nBinsRequired << std::endl;
    std::exit(EXIT_FAILURE);
  }

  for (Int_t edgeCounter = 0; edgeCounter <= m_nBinsRequired; ++edgeCounter) {
    std::cout << "edge[" << edgeCounter << "] = " << Form("%.2f", m_binEdges[edgeCounter]) << std::endl;
  }
}

// void settargetBinSize_mm(Double_t targetBinSize_mm) {
//   m_targetBinSize_mm = targetBinSize_mm;
// }

// void setversionName(TString versionName) {
//   m_versionName = versionName;
// }

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

void calculateBaseEnergyParameters(std::map<Int_t, TH1F*> energyHistogramsMap, TAxis *axisWithCorrectBinning) {
  m_baseEnergy = 0;

  Double_t firstBinX = axisWithCorrectBinning->GetBinCenter(1);
  Int_t minBinForBaseEnergy = axisWithCorrectBinning->FindBin(firstBinX + 0.1*(crack1Position_mm - firstBinX));
  Int_t maxBinForBaseEnergy = axisWithCorrectBinning->FindBin(firstBinX + 0.7*(crack1Position_mm - firstBinX));
    
  std::cout << "minBinForBaseEnergy = " << minBinForBaseEnergy << "; maxBinForBaseEnergy = " << maxBinForBaseEnergy << std::endl;
    
  Int_t binsToCount = 0;
  for (Int_t binCounter = minBinForBaseEnergy; binCounter <= maxBinForBaseEnergy; ++binCounter) {
    if (energyHistogramsMap[binCounter] != 0) {
      TFitResultPtr gaussianFit = (energyHistogramsMap[binCounter])->Fit("gaus", "IMES");
      Double_t fitMean = gaussianFit->Parameter(1);
      Double_t fitWidth = gaussianFit->Parameter(2);
      m_baseEnergy += fitMean;
      m_averageWidth += fitWidth;
      // std::cout << "fitMean = " << fitMean << std::endl;
      ++binsToCount;
    }
  }
  m_baseEnergy = m_baseEnergy/binsToCount;
  m_averageWidth = m_averageWidth/binsToCount;
  m_mipsToGeVConversion = incomingParticleEnergy/m_baseEnergy;
}

void initializeEnergyHistogramsMap(std::map<Int_t, TH1F*> &energyHistogramsMap, std::map<Int_t, std::pair<Double_t, Double_t> > histogramRangesMap) {
  for(Int_t binCounter=0; binCounter<=(1+m_nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    energyHistogramsMap[binCounter] = new TH1F(Form("h_energyHistogram_binNumber_%i",binCounter), Form("number of hits; energy deposited"), numberOfBinsInEnergyHistograms, histogramRangesMap[binCounter].first, histogramRangesMap[binCounter].second);
  }
}

void initializeAllEventsEnergyHistogram(TH1F* &allEventsEnergyHistogram) {
  allEventsEnergyHistogram = new TH1F(Form("allEventsEnergyHistogram"), Form("allEventsEnergyHistogram"), numberOfBinsInEnergyHistograms*5, m_baseEnergy - 7*m_averageWidth, m_baseEnergy + 7*m_averageWidth);
}

void writeEnergyHistogramsToFile(std::map<Int_t, TH1F*> energyHistogramsMap, TH1F *allEventsEnergyHistogram, TFile *histogramsOutputFile) {
  for (Int_t binCounter = 0; binCounter <= (1+m_nBinsRequired); ++binCounter) {
    histogramsOutputFile->WriteTObject(energyHistogramsMap[binCounter]);
  }
  histogramsOutputFile->WriteTObject(allEventsEnergyHistogram);
}

void cleanEnergyHistogramsMap(std::map<Int_t, TH1F*> &energyHistogramsMap) {
  for(Int_t binCounter=0; binCounter<=(1+m_nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    delete energyHistogramsMap[binCounter];
  }
}

void cleanAllEventsEnergyHistogram(TH1F* &allEventsEnergyHistogram) {
  delete allEventsEnergyHistogram;
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
  for(Int_t binCounter=0; binCounter<=(1+m_nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    // std::vector<Double_t> energiesPerEvent[binCounter];
    sumEnergy[binCounter] = 0;
    sumEnergySquare[binCounter] = 0;
    numberOfEvents[binCounter] = 0;
  }
}

void populateHistogramRangesMap(std::map<Int_t, std::pair<Double_t, Double_t> > &histogramRangesMap, std::map<Int_t, Double_t> sumEnergy, std::map<Int_t, Double_t> sumEnergySquare, std::map<Int_t, Int_t> numberOfEvents) {
  for(Int_t binCounter=0; binCounter<=(1+m_nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    Double_t sumEnergyInBin = sumEnergy[binCounter];
    Double_t sumEnergySquareInBin = sumEnergySquare[binCounter];
    Int_t numberOfEventsInBin = numberOfEvents[binCounter];

    Double_t statisticalAverageEnergy = sumEnergyInBin/numberOfEventsInBin;
    Double_t statisticalRMSEnergy = sqrt((sumEnergySquareInBin/numberOfEventsInBin)-statisticalAverageEnergy*statisticalAverageEnergy);
    (histogramRangesMap[binCounter]).first = statisticalAverageEnergy - sigmas_down*statisticalRMSEnergy;
    (histogramRangesMap[binCounter]).second = statisticalAverageEnergy + sigmas_up*statisticalRMSEnergy;
  }
}

void fillEnergyHistogramsMap(std::map<Int_t, TH1F*> &energyHistogramsMap, std::map<Int_t, std::vector<Double_t> > energiesPerEvent) {
  for(Int_t binCounter=0; binCounter<=(1+m_nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    if ((energiesPerEvent[binCounter]).size() > 0) {
      for (std::vector<Double_t>::iterator eventEnergyIterator = (energiesPerEvent[binCounter]).begin(); eventEnergyIterator != (energiesPerEvent[binCounter]).end(); ++eventEnergyIterator) {
        (energyHistogramsMap[binCounter])->Fill(*eventEnergyIterator);
      }
    }
  }
}

void fillAllEventsEnergyHistogram(TH1F *allEventsEnergyHistogram, std::map<Int_t, Double_t> sumEnergy, std::map<Int_t, Int_t> numberOfEvents, std::map<Int_t, std::vector<Double_t> > energiesPerEvent) {
  for(Int_t binCounter=0; binCounter<=(1+m_nBinsRequired); ++binCounter) { // 0: underflow, 1 --> Nbins: regular bins, (1+Nbins): overflow
    Double_t multFactor = (m_baseEnergy*numberOfEvents[binCounter])/sumEnergy[binCounter];
    if ((energiesPerEvent[binCounter]).size() > 0) {
      for (std::vector<Double_t>::iterator eventEnergyIterator = (energiesPerEvent[binCounter]).begin(); eventEnergyIterator != (energiesPerEvent[binCounter]).end(); ++eventEnergyIterator) {
        Double_t scaledEnergy = (*eventEnergyIterator)*multFactor;
        allEventsEnergyHistogram->Fill(scaledEnergy);
      }
    }
  }
}

void printETA(Double_t timeRemainingInSec) {
  Int_t hoursRemaining = static_cast<Int_t>(timeRemainingInSec/3600.);
  Int_t minutesRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining)/60.);
  Int_t secondsRemaining = static_cast<Int_t>((timeRemainingInSec - 3600.*hoursRemaining - 60*minutesRemaining) + 0.5);
  
  std::cout << Form("\rApproximately %02d h: %02d m: %02d s remaining", hoursRemaining, minutesRemaining, secondsRemaining) << std::flush;
}

void populateEnergyHistograms() {
  TChain  *lSimTree = new TChain("HGCSSTree");
  TChain  *lRecTree = new TChain("RecoTree");
  read_input_files (lSimTree, lRecTree);

  std::vector<HGCSSRecoHit> *rechitvec = 0;
  std::vector<HGCSSSimHit> *simhitvec = 0;
  // std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSGenParticle> *genvec = 0; // get the generated particle vector
  HGCSSEvent * event = 0;

  lSimTree->SetBranchAddress("HGCSSEvent",&event);
  lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
  lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  // lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
  lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    
  unsigned nEvts_total = lSimTree->GetEntries();
  std::cout << "number of events: " << nEvts_total << std::endl;
  
  // TH1F *h_generatedParticlesX = new TH1F(Form("h_generatedParticlesX"), Form("number of events; x"), m_nBinsRequired, lowerEdgeFirstBin, upperEdgeLastBin);
  TH1F *h_generatedParticlesX = new TH1F(Form("h_generatedParticlesX"), Form("number of events; x"), m_nBinsRequired, m_binEdges);
  TH1F *h_generatedParticlesY = new TH1F(Form("h_generatedParticlesY"), Form("number of events; y"), m_nBinsRequired, 0, 0);
  TH1F *h_generatedParticlesZ = new TH1F(Form("h_generatedParticlesZ"), Form("number of events; z"), m_nBinsRequired, 0, 0);

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
    // Double_t posX_gen_ini = ((*genvec)[0]).x();
    // Double_t posY_gen_ini = ((*genvec)[0]).y();
    Double_t posX_gen_ini = event->vtx_x();
    Double_t posY_gen_ini = event->vtx_y();
    Double_t posZ_gen_ini = event->vtx_z();
    h_generatedParticlesX->Fill(posX_gen_ini);
    Double_t binNumber = axisWithCorrectBinning->FindBin(posX_gen_ini);
    // std::cout << "Checking... for initial position: " << posX_gen_ini << ", bin number: " << binNumber << ", corresponding bin position: " << axisWithCorrectBinning->GetBinCenter(binNumber) << std::endl;
    analyzeEvent(rechitvec, energiesPerEvent[binNumber], sumEnergy[binNumber], sumEnergySquare[binNumber], numberOfEvents[binNumber]);
    h_generatedParticlesY->Fill(posY_gen_ini);
    h_generatedParticlesZ->Fill(posZ_gen_ini);
  }
  std::cout << std::endl;
  TFile *histogramsOutputFile = new TFile(Form("root_histograms/histograms_")+m_versionName+eta_portion+Form("_results.root"),"RECREATE");
  histogramsOutputFile->WriteTObject(h_generatedParticlesX);
  std::cout << "Printing histograms data for histogram of x-positions of generated particles:" << std::endl;
  printHistogramData(h_generatedParticlesX);
  histogramsOutputFile->WriteTObject(h_generatedParticlesY);
  histogramsOutputFile->WriteTObject(h_generatedParticlesZ);
  // delete h_generatedParticlesX; // not deleted here because its axis is required below
  delete h_generatedParticlesY;
  delete h_generatedParticlesZ;

  std::map<Int_t, std::pair<Double_t, Double_t> > histogramRangesMap;
  populateHistogramRangesMap(histogramRangesMap, sumEnergy, sumEnergySquare, numberOfEvents);

  std::map<Int_t, TH1F*> energyHistogramsMap;
  initializeEnergyHistogramsMap(energyHistogramsMap, histogramRangesMap);
  fillEnergyHistogramsMap(energyHistogramsMap, energiesPerEvent);

  calculateBaseEnergyParameters(energyHistogramsMap, axisWithCorrectBinning);
  delete h_generatedParticlesX; // now safe to delete generatedParticlesX

  // Scale to GeV:
  for (Int_t binCounter = 0; binCounter <= (1+m_nBinsRequired); ++binCounter) {
    sumEnergy[binCounter] *= m_mipsToGeVConversion;
    if ((energiesPerEvent[binCounter]).size() > 0) {
      for (std::vector<Double_t>::iterator eventEnergyIterator = (energiesPerEvent[binCounter]).begin(); eventEnergyIterator != (energiesPerEvent[binCounter]).end(); ++eventEnergyIterator) {
        *eventEnergyIterator *= m_mipsToGeVConversion;
      }
    }
  }

  TH1F *allEventsEnergyHistogram;
  initializeAllEventsEnergyHistogram(allEventsEnergyHistogram);
  fillAllEventsEnergyHistogram(allEventsEnergyHistogram, sumEnergy, numberOfEvents, energiesPerEvent);
  
  writeEnergyHistogramsToFile(energyHistogramsMap, allEventsEnergyHistogram, histogramsOutputFile);

  cleanEnergyHistogramsMap(energyHistogramsMap);
  cleanAllEventsEnergyHistogram(allEventsEnergyHistogram);
  delete histogramsOutputFile;
}

// void populateFromData(Double_t targetBinSize_mm) {
void populateFromData() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");
  // settargetBinSize_mm(targetBinSize_mm);
  // setversionName(Form("cracks_study_%.2fmmgap_averaged_weights", targetBinSize_mm));
  setBinning();
  readWeights();
  // std::cout << "Checking... bin centers: " << std::endl;
  // for (std::vector<Double_t>::iterator binCenters_iterator = binCenters.begin(); binCenters_iterator != binCenters.end(); ++binCenters_iterator) {
  //   std::cout << *binCenters_iterator << std::endl;
  // }
  populateEnergyHistograms();
  std::cout << "Finished populating histograms from data!" << std::endl;
}

# ifndef __CINT__
// int main(int argc, char* argv[]) {
int main() {
  // if (argc != 2) {
  //   std::cout << "populateFromData requires 1 argument; " << (argc-1) << " given" << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }
  // populateFromData(atof(argv[1]));
  populateFromData();
  return 0;
}
# endif
