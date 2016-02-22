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

void plotE_resolution_vlayer_nevts(Int_t version_number, TString datadir, Double_t et, Double_t threshold){//main

  //load the shared library for HGCSS* classes
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

  if(threshold<0.5) {
    std::cout << "Threshold provided less than minimum threshold analyzable from Digi file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TString version_name;

  switch(version_number) {
  case 34:
    version_name=Form("vfrac1");
    break;
  case 40:
    version_name=Form("vfrac75");
    break;
  case 41:
    version_name=Form("vfrac50");
    break;
  case 42:
    version_name=Form("vfrac25");
    break;
  case 43:
    version_name=Form("vflat");
    break;
  default:
    std::cout << "version mapping gone wrong" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  TString datadir_old=Form("/export/cmss2/tmudholk/HGCal/vcnncal");
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
  

  TString HGcal_common_prefix = datadir + Form("/HGcal__version%i_model2_BOFF_",version_number);
  TString HGcal_common_prefix_old = datadir_old + Form("/HGcal__version%i_model2_BOFF_",version_number_old);
  TString Digi_common_prefix = datadir + Form("/Digi__version%i_model2_BOFF_",version_number);
  TString Digi_common_prefix_old = datadir_old + Form("/Digi__version%i_model2_BOFF_",version_number_old);
  TString common_suffix = Form(".root");
  
  //double eta_values_array[] = {1.6,2.1,2.5};
  double eta_values_array[] = {2.1};

  std::vector<double> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(double));

  for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
    
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
    Double_t total_n_hits = 0;
    Double_t total_n_hits_in_layer[24];
    Double_t total_E_in_layer[24];
    for(int layer_counter = 0; layer_counter < 24; layer_counter++) {
      total_n_hits_in_layer[layer_counter] = 0;
      total_E_in_layer[layer_counter] = 0;
    }
    std::vector<Double_t> weighted_energy_of_hit;
    std::vector<Double_t> layer_of_hit;
    Double_t putative_weighted_energy_of_hit_mean=0;
    Double_t putative_weighted_energy_of_hit_stdev=0;

    const unsigned nEvts_old = lSimTree_old->GetEntries();
    Double_t total_n_hits_old = 0;
    Double_t total_n_hits_in_layer_old[24];
    Double_t total_E_in_layer_old[24];
    for(int layer_counter = 0; layer_counter < 24; layer_counter++) {
      total_n_hits_in_layer_old[layer_counter] = 0;
      total_E_in_layer_old[layer_counter] = 0;
    }
    std::vector<Double_t> weighted_energy_of_hit_old;
    std::vector<Double_t> layer_of_hit_old;
    Double_t putative_weighted_energy_of_hit_mean_old=0;
    Double_t putative_weighted_energy_of_hit_stdev_old=0;

    // for (int layer_counter = 0; layer_counter < 24; layer_counter++) {
    //   totalE_in_layer[layer_counter] = 0;
    //   totalE_in_layer_old[layer_counter] = 0;
    // }
    
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      
      if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
      lSimTree->GetEntry(ievt);
      lRecTree->GetEntry(ievt);

      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
    	const HGCSSRecoHit lHit = (*rechitvec)[iH];
    	unsigned layer = lHit.layer();
    	double energy = lHit.energy();
    	
	if(energy>threshold) {
	
	  double weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);
	
	  weighted_energy_of_hit.push_back(weighted_energy);
	  //std::cout << "energy is " << weighted_energy << std::endl;
	  total_n_hits_in_layer[layer] += 1.0;
	  total_E_in_layer[layer] += weighted_energy/nEvts;
	  layer_of_hit.push_back(1.0*layer);
	  total_n_hits += 1.0;
	  putative_weighted_energy_of_hit_mean += weighted_energy;
	  putative_weighted_energy_of_hit_stdev += weighted_energy*weighted_energy;
	}
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

	if(energy>threshold) {
	
	  double weighted_energy = energy*weights_old[layer]/tanh(eta_values[eta_counter]);
	
	  weighted_energy_of_hit_old.push_back(weighted_energy);
	  //std::cout << "energy is " << weighted_energy << std::endl;
	  total_n_hits_in_layer_old[layer] += 1.0;
	  total_E_in_layer_old[layer] += weighted_energy/nEvts_old;
	  layer_of_hit_old.push_back(1.0*layer);
	  total_n_hits_old += 1.0;
	  putative_weighted_energy_of_hit_mean_old += weighted_energy;
	  putative_weighted_energy_of_hit_stdev_old += weighted_energy*weighted_energy;
	}
      }
    }//loop on hits

    putative_weighted_energy_of_hit_mean = putative_weighted_energy_of_hit_mean/weighted_energy_of_hit.size();
    putative_weighted_energy_of_hit_mean_old = putative_weighted_energy_of_hit_mean_old/weighted_energy_of_hit_old.size();
    putative_weighted_energy_of_hit_stdev = putative_weighted_energy_of_hit_stdev/weighted_energy_of_hit.size();
    putative_weighted_energy_of_hit_stdev_old = putative_weighted_energy_of_hit_stdev_old/weighted_energy_of_hit_old.size();
    putative_weighted_energy_of_hit_stdev = putative_weighted_energy_of_hit_stdev - putative_weighted_energy_of_hit_mean*putative_weighted_energy_of_hit_mean;
    putative_weighted_energy_of_hit_stdev_old = putative_weighted_energy_of_hit_stdev_old - putative_weighted_energy_of_hit_mean_old*putative_weighted_energy_of_hit_mean_old;
    putative_weighted_energy_of_hit_stdev = sqrt(putative_weighted_energy_of_hit_stdev);
    putative_weighted_energy_of_hit_stdev_old = sqrt(putative_weighted_energy_of_hit_stdev_old);

    //std::cout << "mean is " << putative_weighted_energy_of_hit_mean << std::endl;
    //std::cout << "stdev is " << putative_weighted_energy_of_hit_stdev << std::endl;

    //Double_t lower_bound_for_hist = putative_weighted_energy_of_hit_mean - 0.5*putative_weighted_energy_of_hit_stdev;
    //Double_t lower_bound_for_hist = 0;
    //std::cout << "lhist " << lower_bound_for_hist << std::endl;
    Double_t upper_bound_for_hist = putative_weighted_energy_of_hit_mean + 5.0*putative_weighted_energy_of_hit_stdev;
    //std::cout << "uhist " << upper_bound_for_hist << std::endl;
    
    //Double_t lower_bound_for_hist_old = putative_weighted_energy_of_hit_mean_old - 0.5*putative_weighted_energy_of_hit_stdev_old;
    //Double_t lower_bound_for_hist_old = 0;
    //std::cout << "lhisto " << lower_bound_for_hist_old << std::endl;
    Double_t upper_bound_for_hist_old = putative_weighted_energy_of_hit_mean_old + 5.0*putative_weighted_energy_of_hit_stdev_old;
    //std::cout << "uhisto " << upper_bound_for_hist_old << std::endl;

    Double_t lower_bound_for_net = 0;
    Double_t upper_bound_for_net;

    if(upper_bound_for_hist > upper_bound_for_hist_old) {
      upper_bound_for_net = upper_bound_for_hist;
    }
    else {
      upper_bound_for_net = upper_bound_for_hist_old;
    }

    //std::cout << "u: " << upper_bound_for_net << std::endl << "l: " << lower_bound_for_net << std::endl;

    TCanvas *myc1 = new TCanvas("Energy Distribution","Energy",800,600);
    myc1->cd();

    // TH2F *p_energy_dist_weighted = new TH2F("Energy Distribution_weighted", "Energies", 500, lower_bound_for_net,upper_bound_for_net,24,-0.5,23.5);
    // for(unsigned hist_filler_counter = 0; hist_filler_counter < weighted_energy_of_hit.size(); hist_filler_counter++) {
    //   p_energy_dist_weighted->Fill(weighted_energy_of_hit[hist_filler_counter],layer_of_hit[hist_filler_counter]);
    // }

    TH1F *p_energy_dist = new TH1F("Energy Distribution", "Energies", 500, lower_bound_for_net,upper_bound_for_net);
    for(unsigned hist_filler_counter = 0; hist_filler_counter < weighted_energy_of_hit.size(); hist_filler_counter++) {
      p_energy_dist->Fill(weighted_energy_of_hit[hist_filler_counter]);
    }

    TH1F *p_energy_dist_old = new TH1F("Energy Distribution_old", "Energies", 500, lower_bound_for_net,upper_bound_for_net);
    for(unsigned hist_filler_counter = 0; hist_filler_counter < weighted_energy_of_hit_old.size(); hist_filler_counter++) {
      p_energy_dist_old->Fill(weighted_energy_of_hit_old[hist_filler_counter]);
    }

    Double_t maxval = 0;
    
    Double_t maxval_freq = p_energy_dist->GetMaximum();
    // //std::cout << "maxval_freq is " << maxval_freq << std::endl;
    Double_t maxval_freq_old = p_energy_dist_old->GetMaximum();
    // //std::cout << "maxval_freq_old is " << maxval_freq_old << std::endl;
    
    if(maxval_freq > maxval_freq_old) {
      maxval = maxval_freq;
    }
    else {
      maxval = maxval_freq_old;
    }
    // std::cout << "maxval is " << maxval << std::endl;
    p_energy_dist->SetMaximum(1.1*maxval);
    p_energy_dist_old->SetMaximum(1.1*maxval);

    //p_energy_dist_weighted->SetTitle(Form("v")+(Form("%i",version_number)+et_portion)+eta_portion);
    //gPad->SetLogy();
    //gPad->SetPhi(290);
    p_energy_dist->SetTitle(version_name+Form("_threshold_%.1f",threshold)+et_portion+eta_portion);
    //p_energy_dist_old->SetTitle(et_portion + eta_portion);
    //p_energy_dist_weighted->Draw("LEGO");
    gPad->SetLogy();
    p_energy_dist->Draw();
    myc1->Update();
    p_energy_dist_old->SetLineColor(2);
    p_energy_dist_old->Draw("same");
    myc1->Update();
    //myc1->SaveAs(Form("plot_shower_profile_energy_distribution_2d_withdigi_v")+(Form("%i",version_number)+et_portion)+eta_portion+".png");
    myc1->SaveAs(Form("plot_energy_distribution_")+version_name+Form("_threshold_%.1f_",threshold)+et_portion+eta_portion+Form(".pdf"));
    
    //delete p_energy_dist_weighted;
    delete p_energy_dist;
    delete p_energy_dist_old;
    delete myc1;
    
    // Double_t ratio_of_norms = total_n_hits/total_n_hits_old;

    maxval=0;
    for (int layer_counter = 0; layer_counter < 24; layer_counter++) {
      if(maxval < total_E_in_layer[layer_counter]) {
    	maxval=total_E_in_layer[layer_counter];
      }
      if(maxval < total_E_in_layer_old[layer_counter]) {
    	maxval=total_E_in_layer_old[layer_counter];
      }
    }
    TCanvas *myc = new TCanvas("Shower Profile","Energy",800,600);
    myc->cd();
    TH1F *p_l = new TH1F("Shower Profile_new", "Profile_new", 24, -0.5, 23.5);
    TH1F *p_l_old = new TH1F("Shower Profile_old", "Profile_old", 24, -0.5, 23.5);
    p_l->SetMaximum(1.1*maxval);
    p_l_old->SetMaximum(1.1*maxval);
    for(unsigned layer_counter = 0; layer_counter < 24; layer_counter++) {
      p_l->Fill(layer_counter,total_E_in_layer[layer_counter]);
      p_l_old->Fill(layer_counter,total_E_in_layer_old[layer_counter]);
    }
    Double_t total_energy = p_l->Integral();
    p_l->SetTitle(version_name+Form("_threshold_%.1f",threshold) + et_portion + eta_portion + Form("totalE = %.1f",total_energy));
    // p_l_old->SetTitle(et_portion + eta_portion + Form("_version_%i",version_number));
    p_l->Draw();
    myc->Update();
    p_l_old->SetLineColor(2);
    p_l_old->Draw("same");
    // p_l_old->Draw();
    myc->Update();
    myc->SaveAs(Form("plot_energy_vs_layer_")+version_name+Form("_threshold_%.1f_",threshold)+et_portion+eta_portion+".pdf");
    
    delete p_l;
    delete p_l_old;
    delete myc;

    weighted_energy_of_hit.clear();
    weighted_energy_of_hit_old.clear();
  }
}
