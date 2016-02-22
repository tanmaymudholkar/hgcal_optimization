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

#include "../../../userlib/include/HGCSSEvent.hh"
#include "../../../userlib/include/HGCSSInfo.hh"
#include "../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../userlib/include/HGCSSSimHit.hh"
#include "../../../userlib/include/HGCSSSamplingSection.hh"

Double_t resolutions_fit(Double_t *x,Double_t *par) {
  return (par[0] + par[1]/sqrt(x[0]));
}

void plotE_resolution(){//main  
   const int layers=28; //number of layers, varied from version to version but remaining constant in one version
   const int truncationwidth=2; // truncate summation of energy in one layer

  //load the shared library for HGCSS* classes
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

  std::vector<double> weights; // vector with double components
  ifstream f_layer_weights; // input file stream
  f_layer_weights.open("layer_weights.dat"); //fopen
  std::string line; // a char string 
  double weight;
  if (f_layer_weights.is_open()) {
    while (getline(f_layer_weights,line)) { //=0 if the last line, line++
      weight = std::atof(line.c_str()); // line.c_str(): a string read from the line, atof: convert string to float
      weights.push_back(weight); // fill the last compoent of vector "weights", with double value "weight" and extend the array by one (like a stack)
    }
  }

  double chisqr;
  double ndfr;
  double chisqdf;

  //double et_values_array[] = {3,5,7,10,30,50,70,100};
  double et_values_array[] = {3,5,7,10,30,50,70,100}; 
  double eta_values_array[] = {1.6,2.1,2.5};

  std::vector<double> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(double)); // set memory range from et_values_array
  std::vector<double> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(double));

  TString HGcal_common_prefix = Form("/export/cmss2/paulini/CMS/HGCal/version30/HGcal_version30_model2_BOFF_");// Form: convert a string to Tstring
  TString Digi_common_prefix = Form("/export/cmss2/paulini/CMS/HGCal/version30/Digi_version30_model2_BOFF_");
  TString common_suffix = Form(".root");

  for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
  //for (unsigned int eta_counter = 1; eta_counter<=1; eta_counter++) {
    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]); 

    std::vector<Double_t> energies_incoming_gev; //calculated by file name
    std::vector<Double_t> energies_incoming_gev_errors;
    std::vector<Double_t> energies;// calculated by summation on events
    std::vector<Double_t> energies_errors;
    std::vector<Double_t> resolutions;
    std::vector<Double_t> resolutions_errors;
    for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
      //for (unsigned int et_counter = 2; et_counter <=5 ; et_counter++) {
      Double_t energy_incoming_gev = et_values[et_counter]*cosh(eta_values[eta_counter]);// total energy E=E_T*ch(eta)
      TString et_portion = Form("et%.0f",et_values[et_counter]); // a string "et30" ....
      double lower_bound_for_hist;
      double upper_bound_for_hist;
      if(energy_incoming_gev < 25) {
	lower_bound_for_hist = 0.4*49*energy_incoming_gev;
	upper_bound_for_hist = 1.2*49*energy_incoming_gev;
      }
      else{
	lower_bound_for_hist = 0.4*49*energy_incoming_gev;
	upper_bound_for_hist = 1.2*49*energy_incoming_gev;
      }
      Double_t energy_incoming_gev_error;
      Double_t mean_energy;
      Double_t mean_energy_error;
      Double_t stddev;
      Double_t stddev_error;
      Double_t resolution;
      Double_t resolution_error;
      std::cout << "et" << et_values[et_counter] << " eta" << eta_values[eta_counter] << std::endl; // printf: et, eta

      TChain  *lSimTree = new TChain("HGCSSTree");
      TChain  *lRecTree = new TChain("RecoTree");

      lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix); // simulated events
      lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);  // reconstructed events
 
      std::vector<HGCSSSamplingSection> * ssvec = 0;
      std::vector<HGCSSSimHit> * simhitvec = 0;
      std::vector<HGCSSRecoHit> * rechitvec = 0;

      lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
      lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
      const unsigned nEvts = lSimTree->GetEntries(); 
      double totalE(0);
      TCanvas *myc = new TCanvas("Energy Distribution","Energy",800,600);
      myc->cd();
      TH1F *p_l = new TH1F("Energy Distribution","Energies",50,lower_bound_for_hist,upper_bound_for_hist); // statistics on energy

      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on events
	totalE = 0;

	if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	lSimTree->GetEntry(ievt);
	lRecTree->GetEntry(ievt);
        
        double centerx[layers]={0},centery[layers]={0},sumx[layers]={0},sumy[layers]={0},sumE[layers]={0};

        for(unsigned iHt(0); iHt<(*rechitvec).size(); ++iHt){ // calculate hit center for one layer
          const HGCSSRecoHit lHit = (*rechitvec)[iHt];
    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  //double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();
          sumx[layer]+=posx*energy;
          sumy[layer]+=posy*energy;
          sumE[layer]+=energy;
        }
        
        for(unsigned ilayer(0); ilayer<layers; ilayer++){
          if(sumE[ilayer]>0){
            centerx[ilayer]=sumx[ilayer]/sumE[ilayer];
            centery[ilayer]=sumy[ilayer]/sumE[ilayer];
            }
        }

        //if(ievt%100==0) std::cout << "eventcenter (" << centerx[15] << "," << centery[15] << ")" << std::endl; //test center
        
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over reconstructed hits
	  const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  //double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();

	  double weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);// weights= thickness, th(eta)=cos(theta), weights/th(eta)=effective thickness
          // energy= energy deposited in one atom , energy*effective thickness= overall energy deposit in one layer

	  if(layer>3 && layer<26 && sumE[layer]>0 && abs(int(posx)/10-(int)(centerx[layer])/10)<=truncationwidth && abs((int)posy/10-(int)(centery[layer])/10)<=truncationwidth) // truncate around energy center 
	    totalE += weighted_energy;
	  
	}
	p_l->Fill(totalE);
      }//loop on events

      p_l->Fit("gaus");
      TF1 *fit = p_l->GetFunction("gaus");
      mean_energy = fit->GetParameter(1);
      mean_energy_error = fit->GetParError(1);
      stddev = fit->GetParameter(2);
      stddev_error = fit->GetParError(2);
      chisqr = fit->GetChisquare();
      ndfr = fit->GetNDF();
      std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
      chisqdf = chisqr/ndfr;
      resolution = stddev/mean_energy;
      resolution_error = resolution*((stddev_error/stddev) + (mean_energy_error/mean_energy));
      energy_incoming_gev_error = energy_incoming_gev*mean_energy_error/mean_energy;

      energies_incoming_gev.push_back(energy_incoming_gev);
      energies_incoming_gev_errors.push_back(energy_incoming_gev_error);
      energies.push_back(mean_energy);
      energies_errors.push_back(mean_energy_error);
      resolutions.push_back(resolution);
      resolutions_errors.push_back(resolution_error);

      p_l->SetTitle(et_portion + eta_portion + Form(". chisq/dof = %3.2f",chisqdf));
      p_l->Draw();
      myc->SaveAs(Form("plot_resolutions_")+et_portion+eta_portion+Form(".png"));
    }// end for et

    TVectorD Tenergies(energies.size(),&energies[0]);
    TVectorD Tresolutions(resolutions.size(),&resolutions[0]);
    TVectorD Tenergies_incoming_gev(energies_incoming_gev.size(),&energies_incoming_gev[0]);
    TVectorD Tenergies_incoming_gev_errors(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
    TVectorD Tenergies_errors(energies_errors.size(),&energies_errors[0]);
    TVectorD Tresolutions_errors(resolutions_errors.size(),&resolutions_errors[0]);

    TCanvas *myc = new TCanvas("Energy mips calibration", "mips versus gev",800,600);
    myc->cd();
    TGraphErrors *energy_mips_calibration = new TGraphErrors(Tenergies_incoming_gev,Tenergies,Tenergies_incoming_gev_errors,Tenergies_errors);
    energy_mips_calibration->Fit("pol1");
    TF1 *calibration_fit = energy_mips_calibration->GetFunction("pol1");
    chisqr = calibration_fit->GetChisquare();
    ndfr = calibration_fit->GetNDF();
    std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
    chisqdf = chisqr/ndfr;
    double conversion_const = calibration_fit->GetParameter(0);
    double conversion_const_error = calibration_fit->GetParError(0);
    double conversion_factor = calibration_fit->GetParameter(1);
    double conversion_factor_error = calibration_fit->GetParError(1);
    TString start_title_conversion_factor = Form("Best Fit: E_mips = (");
    energy_mips_calibration->SetTitle(start_title_conversion_factor + Form("%5.3f +/- %5.3f) * E_GeV + (%5.3f +/- %5.3f). chisq/dof = %3.2f",conversion_factor,conversion_factor_error,conversion_const,conversion_const_error,chisqdf));
    energy_mips_calibration->Draw("AP");
    myc->Update();
    myc->Print(Form("plot_calibration_fit") + eta_portion + Form(".pdf"));
    for (unsigned int energies_incoming_gev_errors_counter = 0; energies_incoming_gev_errors_counter != energies_incoming_gev_errors.size(); energies_incoming_gev_errors_counter++) {
      energies_incoming_gev_errors[energies_incoming_gev_errors_counter] += energies_incoming_gev[energies_incoming_gev_errors_counter]*conversion_factor_error/conversion_factor;
    }
    TVectorD Tenergies_incoming_gev_errors_updated(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
    myc = new TCanvas("Resolutions versus Energy","Resolution versus Energy",800,600);
    myc->cd();
    TGraphErrors *resolutions_versus_energy = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors_updated,Tresolutions_errors);
    TF1 *function_to_fit = new TF1("f_to_fit",resolutions_fit,0.0001,1000,2);
    function_to_fit->SetParameters(0.015,0.1);
    function_to_fit->SetParNames("const_offset","coeff_sqrt");
    resolutions_versus_energy->Fit("f_to_fit");
    chisqr = function_to_fit->GetChisquare();
    ndfr = function_to_fit->GetNDF();
    std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
    chisqdf = chisqr/ndfr;
    TString start_title_resolutions = Form("Best Fit: resolution = (");
    double bestfit_constant = function_to_fit->GetParameter(0);
    double bestfit_constant_error = function_to_fit->GetParError(0);
    double bestfit_coeff = function_to_fit->GetParameter(1);
    double bestfit_coeff_error = function_to_fit->GetParError(1);
    resolutions_versus_energy->SetTitle(start_title_resolutions + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV). chisq/dof = %3.2f",bestfit_constant,bestfit_constant_error,bestfit_coeff,bestfit_coeff_error,chisqdf));
    resolutions_versus_energy->Draw("AP");
    myc->Update();
    myc->Print(Form("plot_resolutions_versus_energy") + eta_portion + Form(".pdf"));                      
   }// end for eta
  
}
