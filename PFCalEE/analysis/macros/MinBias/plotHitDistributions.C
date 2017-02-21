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
#include "TCanvas.h"
#include "TH1F.h"
#include "THStack.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"

// #include "HGCSSEvent.hh"
// #include "HGCSSInfo.hh"
// #include "HGCSSRecoHit.hh"
// #include "HGCSSSimHit.hh"
// #include "HGCSSSamplingSection.hh"
// #include "HGCSSGenParticle.hh"
// #include "HGCSSDetector.hh"
// #include "HGCSSGeometryConversion.hh"

const bool useRaw = true;
const bool getIndividualHitDistributions = false;
const bool getCombinedNormalizedDistribution = true;

const TString plotsDir = Form("plots/");
const TString inputFilePath = Form("/afs/cern.ch/work/t/tmudholk/public/simulation_results/minbias/version_33/PythiaTest/rawHitDistribution_limitedStats_timing.root");
Int_t lowerLayerBoundaries[] = {0,28,40};
Int_t upperLayerBoundaries[] = {27, 39, 51};
// Double_t rawUpperLimits[] = {0.0548*4, 0.0548*4, 0.63*4};
// Double_t rawUpperLimits[] = {0.0548*4, 0.0548*4, 0.0548*4};
Double_t rawUpperLimits[] = {1., 1., 1.};
Double_t digiUpperLimits[] = {5.0, 5.0, 5.0};
Double_t timingCutsRaw[] = {25., 25., 150.};
const Int_t nBinsRaw = 1000;
const Int_t nBinsDigi = 100;
std::map<unsigned, TString> subdetectorNames = {{0, Form("EE")}, {1, Form("FH")}, {2, Form("BH")}};
std::map<unsigned, Int_t> colorsForSubdets = {{0, kBlack}, {1, kRed}, {2, kBlue}};
const Double_t preScale = 5.;

const TString versionName = Form("raw_limitedStats_withTimingCuts");

bool testFile(TString inputPath) {
  TFile *testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    return false;
  }
  return true;
}

void normalize1DHist(TH1F *inputHistogram) {
  Double_t normalizationParameter = 1./inputHistogram->Integral();
  inputHistogram->Scale(normalizationParameter);
}

void plotNormalizedHitDistributions(std::vector<TH1F*> histogramsVector, TFile *outputFile) {
  gStyle->SetOptStat(0);
  TString canvasName = Form("combinedNormalizedHitDistribution");
  TCanvas *hitCanvas = new TCanvas(Form("c_")+canvasName, Form("c_")+canvasName);
  THStack *histogramsCollection = new THStack("histogramsCollection", "Normalized Rechit Distributions for EE, FH, BH");
  if (useRaw) histogramsCollection->SetTitle("Normalized Hit Distributions for EE, FH, BH");

  for (std::vector<TH1F*>::iterator histogramsVectorIterator = histogramsVector.begin(); histogramsVectorIterator != histogramsVector.end(); ++histogramsVectorIterator) {
    histogramsCollection->Add(*histogramsVectorIterator);
  }

  histogramsCollection->Draw("nostack");
  gPad->SetLogy();
  hitCanvas->SaveAs(plotsDir+versionName+Form("_")+canvasName+Form(".png"));
  outputFile->WriteTObject(hitCanvas);
}

void getAndPlotNormalizedHitDistributions(TTree *inputTree, TFile *outputFile) {
  std::vector<TH1F*> histogramsVector;
  for (unsigned layerBoundaryCounter = 0; layerBoundaryCounter < sizeof(lowerLayerBoundaries)/sizeof(lowerLayerBoundaries[0]); ++layerBoundaryCounter) {
    Int_t lowerLayerBoundary = lowerLayerBoundaries[layerBoundaryCounter];
    Int_t upperLayerBoundary = upperLayerBoundaries[layerBoundaryCounter];
    Double_t histUpperLimit = 0.;
    // TString distributionName = Form("hitDistributionBetweenLayers%iAnd%i", layerBoundary, upperLayerBoundary);
    Double_t timingCut = timingCutsRaw[layerBoundaryCounter];
    TString distributionName;
    if (useRaw) {
      histUpperLimit = rawUpperLimits[layerBoundaryCounter];
      distributionName = Form("normalizedHitDistribution_"+subdetectorNames[layerBoundaryCounter]);
      if (timingCut > 0.) inputTree->Draw(Form("energy>>")+distributionName+Form("(%i,0,%.6f)", nBinsRaw, histUpperLimit), Form("energy>0 && layer>=%i && layer<=%i && time<=%.2f && rndm() < 1.0/%.3f", lowerLayerBoundary, upperLayerBoundary, timingCut, preScale));
      else inputTree->Draw(Form("energy>>")+distributionName+Form("(%i,0,%.6f)", nBinsRaw, histUpperLimit), Form("energy>0 && layer>=%i && layer<=%i && rndm() < 1.0/%.3f", lowerLayerBoundary, upperLayerBoundary, preScale));
    }
    else {
      histUpperLimit = digiUpperLimits[layerBoundaryCounter];
      distributionName = Form("normalizedRecHitDistribution_"+subdetectorNames[layerBoundaryCounter]);
      inputTree->Draw(Form("energy>>")+distributionName+Form("(%i,0,%.2f)", nBinsDigi, histUpperLimit), Form("energy>0 && layer>=%i && layer<=%i  && rndm() < 1.0/%.3f", lowerLayerBoundary, upperLayerBoundary, preScale));
    }
    
    TH1F *hitHist = (TH1F*)gDirectory->Get(distributionName);
    hitHist->Scale(1./hitHist->GetEntries());
    hitHist->SetLineColor(colorsForSubdets[layerBoundaryCounter]);
    // histogramsCollection->Add(hitHist);
    histogramsVector.push_back(hitHist);
    
    // if (layerBoundaryCounter == 0) {
    //   if (useRaw) hitHist->SetTitle("Combined Hit Histogram, "+subdetectorNames[layerBoundaryCounter]);
    //   else hitHist->SetTitle("Combined RecHit Histogram, "+subdetectorNames[layerBoundaryCounter]);  
    //   hitHist->Draw();
    // }
    // else hitHist->Draw("same");
  }
  plotNormalizedHitDistributions(histogramsVector, outputFile);
}

void getAndPlotHitDistributions(TTree *inputTree, TFile *outputFile) {
  for (unsigned layerBoundaryCounter = 0; layerBoundaryCounter < sizeof(lowerLayerBoundaries)/sizeof(lowerLayerBoundaries[0]); ++layerBoundaryCounter) {
    Int_t lowerLayerBoundary = lowerLayerBoundaries[layerBoundaryCounter];
    Int_t upperLayerBoundary = upperLayerBoundaries[layerBoundaryCounter];
    Double_t histUpperLimit = 0.;
    // TString distributionName = Form("hitDistributionBetweenLayers%iAnd%i", layerBoundary, upperLayerBoundary);
    Double_t timingCut = timingCutsRaw[layerBoundaryCounter];
    TString distributionName;
    if (useRaw) {
      histUpperLimit = rawUpperLimits[layerBoundaryCounter];
      distributionName = Form("hitDistribution_"+subdetectorNames[layerBoundaryCounter]);
      if (timingCut > 0.) inputTree->Draw(Form("energy>>")+distributionName+Form("(50,0,%.6f)", histUpperLimit), Form("energy>0 && layer>=%i && layer<=%i && time<=%.2f && rndm() < 1.0/%.3f", lowerLayerBoundary, upperLayerBoundary, timingCut, preScale));
      else inputTree->Draw(Form("energy>>")+distributionName+Form("(50,0,%.6f)", histUpperLimit), Form("energy>0 && layer>=%i && layer<=%i && rndm() < 1.0/%.3f", lowerLayerBoundary, upperLayerBoundary, preScale));
    }
    else {
      histUpperLimit = digiUpperLimits[layerBoundaryCounter];
      distributionName = Form("recHitDistribution_"+subdetectorNames[layerBoundaryCounter]);
      inputTree->Draw(Form("energy>>")+distributionName+Form("(100,0,%.2f)", histUpperLimit), Form("energy>0 && layer>=%i && layer<=%i && rndm() < 1.0/%.3f", lowerLayerBoundary, upperLayerBoundary, preScale));
    }
    TCanvas *hitCanvas = new TCanvas(Form("c_")+distributionName, Form("c_")+distributionName);
    TH1F *hitHist = (TH1F*)gDirectory->Get(distributionName);
    if (useRaw) hitHist->SetTitle("Hit Histogram, "+subdetectorNames[layerBoundaryCounter]);
    else hitHist->SetTitle("RecHit Histogram, "+subdetectorNames[layerBoundaryCounter]);
    hitHist->Draw();
    if (!useRaw) hitHist->Fit("landau", "WEM", "", 0.8, 1.5);
    hitCanvas->SaveAs(plotsDir+versionName+Form("_")+distributionName+".png");
    outputFile->WriteTObject(hitCanvas);
    delete hitCanvas;
  }
}

void plotHitDistributions() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_cracks_study/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  gStyle->SetOptFit(1111);

  // Create addresses for variables to put in corresponding branches
  bool dataExists = testFile(inputFilePath);
  if (!dataExists) {
    std::cout << "Input file not found!" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  TFile *inputFile = TFile::Open(inputFilePath);
  TTree *inputTree = (TTree*)inputFile->Get("baseTree");

  TFile *outputFile = new TFile(plotsDir+"collectionHitDistributions.root", "RECREATE");
  if (getIndividualHitDistributions) getAndPlotHitDistributions(inputTree, outputFile);
  if (getCombinedNormalizedDistribution) getAndPlotNormalizedHitDistributions(inputTree, outputFile);
  
  inputFile->Close();
}

# ifndef __CINT__
int main() {
  plotHitDistributions();
  return 0;
}
# endif
