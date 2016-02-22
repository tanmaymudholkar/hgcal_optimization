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

Double_t resolutions_fit(Double_t *x,Double_t *par) {
  return (par[1] + par[0]/sqrt(x[0]));
}

Double_t resolutions_fit_standalone(Double_t *x,Double_t *par) {
  //return sqrt(pow(par[1],2) + pow(par[0],2)/x[0] + pow(par[2]/x[0],2));
  return sqrt(pow(par[1],2) + pow(par[0],2)/x[0]);
}

void plotE_resolution_nocaliboffset(Int_t version_number, TString version_name, TString datadir, Double_t threshold, Double_t sigmas_down, Double_t sigmas_up){//main
  
  //load the shared library for HGCSS* classes
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");
  
  if(threshold<0.5) {
    std::cout << "Threshold provided less than minimum threshold analyzable from Digi file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  std::vector<Double_t> weights;
  ifstream f_layer_weights;
  f_layer_weights.open(datadir+Form("/layer_weights.dat"));
  std::string line;
  Double_t weight;
  if (f_layer_weights.is_open()) {
    while (getline(f_layer_weights,line)) {
      weight = std::atof(line.c_str());
      weights.push_back(weight);
    }
  }
  
  Double_t chisqr;
  Double_t ndfr;
  Double_t chisqdf;
  
  Double_t et_values_array[] = {3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  //Double_t et_values_array[] = {3,5,7,10,20,30,40};
  //Double_t et_values_array[] = {3};
  Double_t eta_values_array[] = {2.1};
  
  std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
  std::vector<Double_t> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(Double_t));
  
  //TString str_threshold = Form(".6f",threshold);
  TString HGcal_common_prefix = datadir + Form("/HGcal__version%i_model2_BOFF_",version_number);
  TString Digi_common_prefix = datadir + Form("/Digi__version%i_model2_BOFF_",version_number);
  TString common_suffix = Form(".root");
  
  for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]);
    std::vector<Double_t> energies_incoming_gev;
    std::vector<Double_t> energies_incoming_gev_errors;
    std::vector<Double_t> mean_energies_wmips;
    std::vector<Double_t> mean_energies_wmips_errors;
    std::vector<Double_t> sigmas_wmips;
    std::vector<Double_t> sigmas_wmips_errors;
    std::vector<Double_t> resolutions;
    std::vector<Double_t> resolutions_errors;
    Double_t resolution;
    Double_t resolution_error;
    for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
      Double_t energy_incoming_gev = et_values[et_counter]*cosh(eta_values[eta_counter]);
      TString et_portion = Form("et%.0f",et_values[et_counter]);
      Double_t energy_incoming_gev_error;
      Double_t mean_energy_wmips;
      Double_t mean_energy_wmips_error;
      Double_t sigma_wmips;
      Double_t sigma_wmips_error;
      std::vector<Double_t> total_energies;
      Double_t meanE_anticipated = 0;
      Double_t sigE_anticipated = 0;
      std::cout << "et" << et_values[et_counter] << " eta" << eta_values[eta_counter] << std::endl;
      
      TChain  *lSimTree = new TChain("HGCSSTree");
      TChain  *lRecTree = new TChain("RecoTree");

      lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
      lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
      
      std::vector<HGCSSSamplingSection> * ssvec = 0;
      std::vector<HGCSSSimHit> * simhitvec = 0;
      std::vector<HGCSSRecoHit> * rechitvec = 0;
      
      lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
      lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
      
      const unsigned nEvts = lSimTree->GetEntries(); 
      Double_t totalE(0);
      
      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
	totalE = 0;
	
	if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	lSimTree->GetEntry(ievt);
	lRecTree->GetEntry(ievt);
	
	//for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop over rechits
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
	  //const HGCSSSimHit lHit = (*simhitvec)[iH];
	  const HGCSSRecoHit lHit = (*rechitvec)[iH];
	  
	  // Double_t posx = lHit.get_x();
	  // Double_t posy = lHit.get_y();
	  //Double_t posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  Double_t energy = lHit.energy();
	  //	  std::cout << "energy is " << energy << "   " << "layer is " << layer << std::endl;
	  
	  if(energy>threshold) {
	    Double_t weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);
	    totalE += weighted_energy;
	  }
	  // weighted_x[layer].push_back(weighted_energy*posx);
	  // weighted_y[layer].push_back(weighted_energy*posy);
	  // totalE_in_layer[layer] += weighted_energy;
	  //p_l_unweighted->Fill(layer);
	}
	
	total_energies.push_back(totalE);
	meanE_anticipated += totalE/nEvts;
	sigE_anticipated += totalE*totalE/nEvts;
      }//loop on hits
      
      sigE_anticipated = sigE_anticipated - meanE_anticipated*meanE_anticipated;
      sigE_anticipated = sqrt(sigE_anticipated);

      Double_t lower_bound_for_hist;
      Double_t upper_bound_for_hist;
      
      lower_bound_for_hist = meanE_anticipated - sigmas_down*sigE_anticipated;
      upper_bound_for_hist = meanE_anticipated + sigmas_up*sigE_anticipated;
            
      TCanvas *myc = new TCanvas("Energy Distribution","Energy",800,600);
      myc->cd();
      TH1F *p_l = new TH1F("Energy Distribution", "Energies", 50, lower_bound_for_hist,upper_bound_for_hist);
      for(unsigned hist_filler_counter = 0; hist_filler_counter < total_energies.size(); hist_filler_counter++) {
	p_l->Fill(total_energies[hist_filler_counter]);
      }
      p_l->Fit("gaus","IME");
      TF1 *gaussian_fit = p_l->GetFunction("gaus");
      mean_energy_wmips = gaussian_fit->GetParameter(1);
      mean_energy_wmips_error = gaussian_fit->GetParError(1);
      sigma_wmips = gaussian_fit->GetParameter(2);
      sigma_wmips_error = gaussian_fit->GetParError(2);
      chisqr = gaussian_fit->GetChisquare();
      ndfr = gaussian_fit->GetNDF();
      std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
      chisqdf = chisqr/ndfr;
      //resolution = stddev/mean_energy;
      //resolution_error = resolution*((stddev_error/stddev) + (mean_energy_error/mean_energy));
      //resolution_error = resolution*sqrt(pow(stddev_error/stddev,2) + pow(mean_energy_error/mean_energy,2));
      //energy_incoming_gev_error = energy_incoming_gev*mean_energy_error/mean_energy;
      energy_incoming_gev_error = 0;
      energies_incoming_gev.push_back(energy_incoming_gev);
      energies_incoming_gev_errors.push_back(energy_incoming_gev_error);
      mean_energies_wmips.push_back(mean_energy_wmips);
      mean_energies_wmips_errors.push_back(mean_energy_wmips_error);
      sigmas_wmips.push_back(sigma_wmips);
      sigmas_wmips_errors.push_back(sigma_wmips_error);
      //resolutions.push_back(resolution);
      //resolutions_errors.push_back(resolution_error);
      
      p_l->SetTitle(version_name + Form("_") + et_portion + eta_portion + Form(". chisq/dof = %3.2f",chisqdf) + Form("_threshold=%.1f",threshold));
      p_l->Draw();
      //p_l_unweighted->Draw();
      //p_E->Draw();
      myc->Print(Form("plot_")+version_name+Form("_nocaliboffset_threshold_%.1f_resolutions_",threshold)+et_portion+eta_portion+Form("_down%.1f_up%.1f.png",sigmas_down,sigmas_up));
      //cout << "min layer = " << min_layer << std::endl;
      //cout << "max layer = " << max_layer << std::endl;
      delete gaussian_fit;
      delete p_l;
      delete myc;
    }
    
    TVectorD Tmean_energies_wmips(mean_energies_wmips.size(),&mean_energies_wmips[0]);
    TVectorD Tmean_energies_wmips_errors(mean_energies_wmips_errors.size(),&mean_energies_wmips_errors[0]);
    TVectorD Tenergies_incoming_gev(energies_incoming_gev.size(),&energies_incoming_gev[0]);
    TVectorD Tenergies_incoming_gev_errors(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
        
    TCanvas *myc = new TCanvas("Energy mips calibration", "mips versus gev",800,600);
    myc->cd();
    TGraphErrors *energy_mips_calibration = new TGraphErrors(Tenergies_incoming_gev,Tmean_energies_wmips,Tenergies_incoming_gev_errors,Tmean_energies_wmips_errors);
    energy_mips_calibration->Fit("pol1");
    TF1 *calibration_fit = energy_mips_calibration->GetFunction("pol1");
    chisqr = calibration_fit->GetChisquare();
    ndfr = calibration_fit->GetNDF();
    std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
    chisqdf = chisqr/ndfr;
    //Double_t offset = calibration_fit->GetParameter(0);
    //Double_t offset_error = calibration_fit->GetParError(0);
    Double_t offset = 0;
    Double_t offset_error = 0;
    Double_t conversion_factor = calibration_fit->GetParameter(1);
    Double_t conversion_factor_error = calibration_fit->GetParError(1);
    TString start_title_conversion_factor = Form("E_mips = (");
    energy_mips_calibration->SetTitle(start_title_conversion_factor + Form("%5.3f +/- %5.3f) * E_GeV + (%.1f +/- %.1f). chisq/dof = %.2f",conversion_factor,conversion_factor_error,offset,offset_error,chisqdf));
    energy_mips_calibration->Draw("AP");
    myc->Update();
    //myc->Print((Form("plot_thr_")+str_threshold+Form("_v")+(Form("%i",version_number)+(Form("_calibration_fit") + eta_portion))) + Form(".pdf"));
    myc->Print(Form("plot_")+version_name+Form("_nocaliboffset_threshold_%.1f_calibration_fit",threshold)+ eta_portion + Form("_down%.1f_up%.1f.pdf",sigmas_down,sigmas_up));
    delete calibration_fit;
    delete energy_mips_calibration;
    delete myc;
    for (unsigned int sigmas_wmips_counter = 0; sigmas_wmips_counter != sigmas_wmips.size(); sigmas_wmips_counter++) {
      // resolution = sigmas_wmips[sigmas_wmips_counter]/(mean_energies_wmips[sigmas_wmips_counter]-offset);
      // resolution_error = resolution*sqrt(pow(sigmas_wmips_errors[sigmas_wmips_counter]/sigmas_wmips[sigmas_wmips_counter],2)+(pow(mean_energies_wmips_errors[sigmas_wmips_counter],2)+pow(offset_error,2))/pow(mean_energies_wmips[sigmas_wmips_counter]-offset,2));
      resolution = sigmas_wmips[sigmas_wmips_counter]/(mean_energies_wmips[sigmas_wmips_counter]);
      resolution_error = resolution*sqrt(pow(sigmas_wmips_errors[sigmas_wmips_counter]/sigmas_wmips[sigmas_wmips_counter],2)+pow(mean_energies_wmips_errors[sigmas_wmips_counter]/mean_energies_wmips[sigmas_wmips_counter],2));      
      resolutions.push_back(resolution);
      resolutions_errors.push_back(resolution_error);
    }
    TVectorD Tresolutions(resolutions.size(),&resolutions[0]);
    TVectorD Tresolutions_errors(resolutions_errors.size(),&resolutions_errors[0]);
    ofstream outfile;
    //TString outfile_name_prefix = Form("data_thr")+str_threshold+Form("_eta_")+eta_portion+Form("_v");
    TString outfile_name_prefix = Form("data_eta_")+eta_portion+Form("_");
    outfile.open(outfile_name_prefix+version_name+Form("_nocaliboffset_threshold_%.1f_down%.1f_up%.1f_resolutions",threshold,sigmas_down,sigmas_up));
    for(unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
      outfile << energies_incoming_gev[et_counter] << "   " << energies_incoming_gev_errors[et_counter] << "   " << resolutions[et_counter] << "   " << resolutions_errors[et_counter] << std::endl;
    }
    outfile.close();
    
    // TVectorD Tenergies_incoming_gev_errors_updated(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
    
    myc = new TCanvas("Resolutions versus Energy","Resolution versus Energy",800,600);
    myc->cd();
    TGraphErrors *resolutions_versus_energy = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors,Tresolutions_errors);
    TF1 *function_to_fit = new TF1("f_to_fit",resolutions_fit,0.0001,1300,2);
    function_to_fit->SetParameters(0.25,0.01);
    function_to_fit->SetParNames("stoch_term","const_term");
    resolutions_versus_energy->Fit(function_to_fit,"IME");
    chisqr = function_to_fit->GetChisquare();
    ndfr = function_to_fit->GetNDF();
    std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
    chisqdf = chisqr/ndfr;
    Double_t fit_stoch = function_to_fit->GetParameter(0);
    Double_t fit_stoch_error = function_to_fit->GetParError(0);
    Double_t fit_const = function_to_fit->GetParameter(1);
    Double_t fit_const_error = function_to_fit->GetParError(1);
    TString start_title_resolutions = version_name + eta_portion + Form(": res = (");
    resolutions_versus_energy->SetTitle(start_title_resolutions + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV). chisq/dof = %.1f",fit_const,fit_const_error,fit_stoch,fit_stoch_error,chisqdf));
    resolutions_versus_energy->Draw("AP");
    myc->Update();
    myc->Print((Form("plot_")+(version_name+(Form("_nocaliboffset_resolutions_versus_energy") + eta_portion))) + Form("_down%.1f_up%.1f.pdf",sigmas_down,sigmas_up));
    delete resolutions_versus_energy;
    delete function_to_fit;
    delete myc;
    std::cout << std::endl;
    
    myc = new TCanvas("Resolutions versus Energy_standalone","Resolution versus Energy",800,600);
    myc->cd();
    TGraphErrors *resolutions_versus_energy_standalone = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors,Tresolutions_errors);
    TF1 *function_to_fit_standalone = new TF1("f_to_fit_standalone",resolutions_fit_standalone,0.0001,1300,2);
    function_to_fit_standalone->SetParameters(0.25,0.01);
    function_to_fit_standalone->SetParNames("stoch_term","const_term");
    resolutions_versus_energy_standalone->Fit(function_to_fit_standalone,"IME");
    chisqr = function_to_fit_standalone->GetChisquare();
    ndfr = function_to_fit_standalone->GetNDF();
    std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
    chisqdf = chisqr/ndfr;
    fit_stoch = function_to_fit_standalone->GetParameter(0);
    fit_stoch_error = function_to_fit_standalone->GetParError(0);
    fit_const = function_to_fit_standalone->GetParameter(1);
    fit_const_error = function_to_fit_standalone->GetParError(1);
    start_title_resolutions = version_name + eta_portion + Form(": res = (");
    resolutions_versus_energy_standalone->SetTitle(start_title_resolutions + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV). chisq/dof = %.1f",fit_const,fit_const_error,fit_stoch,fit_stoch_error,chisqdf));
    resolutions_versus_energy_standalone->Draw("AP");
    myc->Update();
    myc->Print((Form("plot_")+(version_name+(Form("_nocaliboffset_resolutions_versus_energy_standalone") + eta_portion))) + Form("_down%.1f_up%.1f.pdf",sigmas_down,sigmas_up));
    delete resolutions_versus_energy_standalone;
    delete function_to_fit_standalone;
    delete myc;
    std::cout << std::endl;
  }
}
