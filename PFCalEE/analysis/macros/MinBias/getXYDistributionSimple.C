#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cstdlib>
#include<map>
#include<ctime>
// #include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TAxis.h"

// #include "HGCSSEvent.hh"
// #include "HGCSSInfo.hh"
// #include "HGCSSRecoHit.hh"
// // #include "HGCSSSimHit.hh"
// #include "HGCSSSamplingSection.hh"
// #include "HGCSSGenParticle.hh"
// #include "HGCSSDetector.hh"
// #include "HGCSSGeometryConversion.hh"

// const TString dataDir = Form("/afs/cern.ch/user/t/tmudholk/private/mount/eos_mount/cms/store/cmst3/group/hgcal/HGCalMinbias/PythiaTest/");

const TString inputFileName = Form("root://eoscms//eos/cms/store/cmst3/group/hgcal/BHStudies_Standalone/HepMCTest_VBFHtoBB/PythiaTest/DiginoiseOff_IC3_version33_squareGeometry.root");
const TString versionName = Form("VBFHtoBB_modifiedGeometry");

// void fillTree(std::vector<HGCSSSimHit> *simhitvec, TTree *outputTree, Double_t &energy, Double_t &xpos, Double_t &ypos, unsigned &layer, unsigned &cellid, HGCSSGeometryConversion &geomConv) {
//   unsigned prevLayer = 10000;
//   bool isScint = false;

//   for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){ // loop over simhits
//     HGCSSSimHit lHit = (*simhitvec)[iH];
//     energy = lHit.energy();
//     if (energy > 0) {
//       layer = lHit.layer();
//       cellid = lHit.cellid();
//       if (layer != prevLayer) {
//         const HGCSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
//         isScint = subdet.isScint;
//         prevLayer = layer;
//       }
      
//       std::pair<double,double> xy = lHit.get_xy(isScint,geomConv);
//       xpos = xy.first;
//       ypos = xy.second;
          
//       // xpos = 0.;
//       // ypos = 0.;
    
//       // xpos = lHit.get_x();
//       // ypos = lHit.get_y();
    
//       // if (eta < 10) {
//       outputTree->Fill();
//       // if (layer > 42) std::cout << "layer number " << layer << std::endl;
//       // }
//       // if (layer >=39 && layer <= 41) {
//       //   std::cout << "energy = " << energy << ", xpos = " << xpos << ", ypos = " << ypos << ", layer = " << layer << ", cellid = " << cellid <<  std::endl;
//       //   std::cout << "isScint: " << (isScint? "true" : "false") << std::endl;
//       // }
//     }
//   }
// }

void getXYDistributionSimple() {
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_minbias_modified_geometry/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  TFile *inputFile = TFile::Open(inputFileName);
  TTree *inputTree = inputFile->Get("RecoTree");




  inputTree->Draw("HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>xyDist_allLayers");

  TGraph *xyDist_allLayers = gDirectory->Get("xyDist_allLayers");

  TCanvas *c_xyDist_allLayers = new TCanvas("c_xyDist_allLayers", "c_xyDist_allLayers", 2048, 1536);
  c_xyDist_allLayers->cd();

  xyDist_allLayers->Draw();

  c_xyDist_allLayers->SaveAs("plots/" + versionName + "_xyDist_allLayers.png");




  inputTree->Draw("HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>xyDist_BH", "HGCSSRecoHitVec.layer_>=40");

  TGraph *xyDist_BH = gDirectory->Get("xyDist_BH");

  TCanvas *c_xyDist_BH = new TCanvas("c_xyDist_BH", "c_xyDist_BH", 2048, 1536);
  c_xyDist_BH->cd();

  xyDist_BH->Draw();

  c_xyDist_BH->SaveAs("plots/" + versionName + "_xyDist_BH.png");



  inputTree->Draw("HGCSSRecoHitVec.ypos_:HGCSSRecoHitVec.xpos_>>xyDist_BH_ZoomedIn", "HGCSSRecoHitVec.layer_>=40");

  TGraph *xyDist_BH_ZoomedIn = gDirectory->Get("xyDist_BH_ZoomedIn");

  TCanvas *c_xyDist_BH_ZoomedIn = new TCanvas("c_xyDist_BH_ZoomedIn", "c_xyDist_BH_ZoomedIn", 2048, 1536);
  c_xyDist_BH_ZoomedIn->cd();

  xyDist_BH_ZoomedIn->Draw();
  // xyDist_BH_ZoomedIn->GetXaxis()->SetRangeUser(-100.,0.);
  // xyDist_BH_ZoomedIn->GetYaxis()->SetRangeUser(-400.,-300.);

  c_xyDist_BH_ZoomedIn->SaveAs("plots/" + versionName + "_xyDist_BH_ZoomedIn.png");




  inputFile->Close();
}
