#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<stdlib.h>

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

#include "../../userlib/include/HGCSSEvent.hh"
#include "../../userlib/include/HGCSSInfo.hh"
#include "../../userlib/include/HGCSSRecoHit.hh"
#include "../../userlib/include/HGCSSSimHit.hh"
#include "../../userlib/include/HGCSSSamplingSection.hh"

void plotE_resolution_vlayer(Int_t version_number, TString datadir, Double_t et){//main

  //load the shared library for HGCSS* classes
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

  TString datadir_old=Form("/export/cmss2/tmudholk/HGCal/version34");
  Int_t version_number_old = 34;
  
  std::vector<double> weights;
  std::vector<double> weights_old;
  
  ifstream f_layer_weights;
  std::string line;
  double weight;

  f_layer_weights.open(datadir+Form("/layer_weights.dat"));
  if (f_layer_weights.is_open()) {
    while (getline(f_layer_weights,line)) {
      weight = std::atof(line.c_str());
      weights.push_back(weight);
    }
  }
  f_layer_weights.close();

  f_layer_weights.open(datadir_old+Form("/layer_weights.dat"));
  if (f_layer_weights.is_open()) {
    while (getline(f_layer_weights,line)) {
      weight = std::atof(line.c_str());
      weights_old.push_back(weight);
    }
  }
  f_layer_weights.close();
  

  TString HGcal_common_prefix = datadir + Form("/HGcal__version")+Form("%i",version_number)+Form("_model2_BOFF_");
  TString HGcal_common_prefix_old = datadir_old + Form("/HGcal__version")+Form("%i",version_number_old)+Form("_model2_BOFF_");
  TString Digi_common_prefix = datadir + Form("/Digi__version")+Form("%i",version_number)+Form("_model2_BOFF_");
  TString Digi_common_prefix_old = datadir_old + Form("/Digi__version")+Form("%i",version_number_old)+Form("_model2_BOFF_");
  TString common_suffix = Form(".root");
  
  double eta_values_array[] = {1.6,2.1,2.5};

  std::vector<double> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(double));

  for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
    TCanvas *myc = new TCanvas("Shower Profile","Energy",800,600);
    myc->cd();

    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]);
    TString et_portion = Form("et%.0f",et);
    std::cout << "eta" << eta_values[eta_counter] << std::endl;

    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lSimTree_old = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");
    TChain  *lRecTree_old = new TChain("RecoTree");

    lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
    lSimTree_old->AddFile(HGcal_common_prefix_old+et_portion+eta_portion+common_suffix);
    lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
    lRecTree_old->AddFile(Digi_common_prefix_old+et_portion+eta_portion+common_suffix);
 
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    std::vector<HGCSSSamplingSection> * ssvec_old = 0;
    std::vector<HGCSSSimHit> * simhitvec_old = 0;
    std::vector<HGCSSRecoHit> * rechitvec_old = 0;

    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);

    lSimTree_old->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec_old);
    lSimTree_old->SetBranchAddress("HGCSSSimHitVec",&simhitvec_old);
    lRecTree_old->SetBranchAddress("HGCSSRecoHitVec",&rechitvec_old);
  
    const unsigned nEvts = lSimTree->GetEntries();
    Double_t totalE_overall = 0;
    Double_t totalE_in_layer[24];

    const unsigned nEvts_old = lSimTree_old->GetEntries();
    Double_t totalE_overall_old = 0;
    Double_t totalE_in_layer_old[24];

    for (int layer_counter = 0; layer_counter < 24; layer_counter++) {
      totalE_in_layer[layer_counter] = 0;
      totalE_in_layer_old[layer_counter] = 0;
    }

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      
      if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
      lSimTree->GetEntry(ievt);
      lRecTree->GetEntry(ievt);

      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
	const HGCSSRecoHit lHit = (*rechitvec)[iH];
	unsigned layer = lHit.layer();
	double energy = lHit.energy();
	double weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);
	  
	totalE_in_layer[layer] += weighted_energy/nEvts;
	totalE_overall += weighted_energy/nEvts;
      }
    }//loop on hits

    for (unsigned ievt_old(0); ievt_old<nEvts_old; ++ievt_old){//loop on entries
      
      if (ievt_old%100==0) std::cout << " -- Processing event " << ievt_old << std::endl;
      lSimTree_old->GetEntry(ievt_old);
      lRecTree_old->GetEntry(ievt_old);

      for (unsigned iH_old(0); iH_old<(*rechitvec_old).size(); ++iH_old){//loop over rechits
	const HGCSSRecoHit lHit_old = (*rechitvec_old)[iH_old];
	unsigned layer = lHit_old.layer();
	double energy = lHit_old.energy();
	double weighted_energy = energy*weights_old[layer]/tanh(eta_values[eta_counter]);
	  
	totalE_in_layer_old[layer] += weighted_energy/nEvts_old;
	totalE_overall_old += weighted_energy/nEvts_old;
      }
    }//loop on hits

    Double_t maxval=0;
    for (int layer_counter = 0; layer_counter < 24; layer_counter++) {
      //totalE_in_layer[layer_counter] = totalE_in_layer[layer_counter]/totalE_overall;
      //totalE_in_layer_old[layer_counter] = totalE_in_layer_old[layer_counter]/totalE_overall_old;
      if(maxval < totalE_in_layer[layer_counter]) {
	maxval=totalE_in_layer[layer_counter];
      }
      if(maxval < totalE_in_layer_old[layer_counter]) {
	maxval=totalE_in_layer_old[layer_counter];
      }
      //std::cout << "totalE_frac: " << totalE_in_layer[layer_counter] << std::endl;
      //std::cout << "totalE_old_frac: " << totalE_in_layer_old[layer_counter] << std::endl;
      //std::cout << "maxval is now " << maxval << std::endl;
    }

    //std::cout << "maxval is " << maxval << std::endl;

    Double_t ratio_of_norms = totalE_overall/totalE_overall_old;

    TH1F *p_l = new TH1F("Shower Profile_new", "Profile_new", 24, -0.5, 23.5);
    TH1F *p_l_old = new TH1F("Shower Profile_old", "Profile_old", 24, -0.5, 23.5);
    p_l->SetMaximum(1.1*maxval);
    p_l_old->SetMaximum(1.1*maxval);
    for(unsigned layer_counter = 0; layer_counter < 24; layer_counter++) {
      p_l->Fill(layer_counter,totalE_in_layer[layer_counter]);
      p_l_old->Fill(layer_counter,totalE_in_layer_old[layer_counter]);
    }
    p_l->SetTitle(et_portion + eta_portion + Form("Ratio is: %.3f",ratio_of_norms));
    p_l_old->SetTitle(et_portion + eta_portion);
    p_l->Draw();
    myc->Update();
    p_l_old->SetLineColor(2);
    p_l_old->Draw("same");
    myc->Update();
    myc->SaveAs(Form("plot_shower_profile_withdigi_v")+(Form("%i",version_number)+et_portion)+eta_portion+".png");
    delete myc;
    delete p_l;
    delete p_l_old;
  }
}
