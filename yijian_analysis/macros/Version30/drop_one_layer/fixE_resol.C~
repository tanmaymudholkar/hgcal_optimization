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

void plotE_resolution(int droplayer,int etcounter,double &resol,double &resol_err){
   TString droplayer_portion=Form("%d",droplayer);
   std::cout << "layer="  << droplayer << std::endl;
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

 // for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
  for (unsigned int eta_counter = 1; eta_counter<=1; eta_counter++) {
    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]); 

    std::vector<Double_t> energies_incoming_gev; //calculated by file name
    std::vector<Double_t> energies_incoming_gev_errors;
    std::vector<Double_t> energies;// calculated by summation on events
    std::vector<Double_t> energies_errors;
    std::vector<Double_t> resolutions;
    std::vector<Double_t> resolutions_errors;
    //for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
      for (unsigned int et_counter = etcounter; et_counter <= etcounter ; et_counter++) {
      Double_t energy_incoming_gev = et_values[et_counter]*cosh(eta_values[eta_counter]);// total energy E=E_T*ch(eta)
      TString et_portion = Form("et%.0f",et_values[et_counter]); // a string "et30" ....
      double lower_bound_for_hist;
      double upper_bound_for_hist;
      if(energy_incoming_gev < 25) {
	lower_bound_for_hist = 0.5*49*energy_incoming_gev;
	upper_bound_for_hist = 1.4*49*energy_incoming_gev;
      }
      else{
	lower_bound_for_hist = 0.5*49*energy_incoming_gev;
	upper_bound_for_hist = 1.4*49*energy_incoming_gev;
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

	  if(sumE[layer]>0 && layer!=droplayer)// && abs(int(posx)/10-(int)(centerx[layer])/10)<=truncationwidth && abs((int)posy/10-(int)(centery[layer])/10)<=truncationwidth) // truncate around energy center 
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

      resol=resolution;
      resol_err=resolution_error;
    }// end for et

   }// end for eta
  
}
void fixE_resol(){
   const int layers=28;
   TCanvas *myc=new TCanvas("myc","myc",1000,1000);
   myc->Divide(3,1);
   for(int et=2;et<=6;et+=2){
   myc->cd(et/2);
   std::vector<double> resol;
   std::vector<double> resol_err;
   std::vector<double> layer_array;
   std::vector<double> layer_array_err;
   double one_resol,one_err;
   for(int i=0;i<28;i++){
       plotE_resolution(i,et,one_resol,one_err);
       layer_array.push_back((double)i);
       resol.push_back(one_resol);
       resol_err.push_back(one_err);
       layer_array_err.push_back(0);
   }
   TVectorD Tresol(resol.size(),&resol[0]);
   TVectorD Tresol_err(resol_err.size(),&resol_err[0]);
   TVectorD Tlayer_array(layer_array.size(),&layer_array[0]);
   TVectorD Tlayer_array_err(layer_array_err.size(),&layer_array_err[0]);
   TGraphErrors *fixE_resol_layer = new TGraphErrors(Tlayer_array,Tresol,Tlayer_array_err,Tresol_err);
   fixE_resol_layer->SetTitle("et7/30/70,eta2.1,resol_droplayer");
   //TAxis *ax=fixE_resol_layer->GetYaxis();
   //ax->SetLimits(0.01,0.07);
   fixE_resol_layer->Draw("AC*");
   }
}
