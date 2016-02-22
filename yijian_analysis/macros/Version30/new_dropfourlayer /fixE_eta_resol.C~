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

#include "../../../../userlib/include/HGCSSEvent.hh"
#include "../../../../userlib/include/HGCSSInfo.hh"
#include "../../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../../userlib/include/HGCSSSimHit.hh"
#include "../../../../userlib/include/HGCSSSamplingSection.hh"

Double_t resolutions_fit(Double_t *x,Double_t *par) {
  return (par[0] + par[1]/sqrt(x[0]));
}

void plotE_resolution(int droplayer,int etcounter,int etacounter, double &resol, double &resol_err, double &dropenergy){//main  
  const int layers=28;
  const int truncationwidth=2;
  //load the shared library for HGCSS* classes
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

  std::vector<double> weights;
  ifstream f_layer_weights;
  f_layer_weights.open("layer_weights.dat");
  std::string line;
  double weight;
  if (f_layer_weights.is_open()) {
    while (getline(f_layer_weights,line)) {
      weight = std::atof(line.c_str());
      weights.push_back(weight);
    }
  }

  double chisqr;
  double ndfr;
  double chisqdf;

  double et_values_array[] = {3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  double eta_values_array[] = {1.6,2.1,2.5};

  std::vector<double> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(double));
  std::vector<double> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(double));

  TString HGcal_common_prefix = Form("/export/cmss2/tmudholk/HGCal/version30/HGcal__version30_model2_BOFF_");
  TString Digi_common_prefix = Form("/export/cmss2/tmudholk/HGCal/version30/Digi__version30_model2_BOFF_");
  TString common_suffix = Form(".root");

  //for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
   for(int eta_counter=etacounter; eta_counter<=etacounter; eta_counter++){  
    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]);
    std::vector<Double_t> energies_incoming_gev;
    std::vector<Double_t> energies_incoming_gev_errors;
    std::vector<Double_t> energies;
    std::vector<Double_t> energies_errors;
    std::vector<Double_t> resolutions;
    std::vector<Double_t> resolutions_errors;
    //for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
     for(int et_counter=etcounter; et_counter<=etcounter; et_counter++){
      Double_t energy_incoming_gev = et_values[et_counter]*cosh(eta_values[eta_counter]);
      TString et_portion = Form("et%.0f",et_values[et_counter]);
      Double_t energy_incoming_gev_error;
      Double_t mean_energy;
      Double_t mean_energy_error;
      Double_t stddev;
      Double_t stddev_error;
      Double_t resolution;
      Double_t resolution_error;
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
      double totalE(0);

      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
	totalE = 0;
        dropenergy=0;
	if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	lSimTree->GetEntry(ievt);
	lRecTree->GetEntry(ievt);
	
        double sumx[layers]={0},sumy[layers]={0},sumE[layers]={0},centerx[layers]={0},centery[layers]={0};
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits to get energycenter
	  const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  //double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();
          sumx[layer]+=posx*energy;
          sumy[layer]+=posy*energy;
          sumE[layer]+=energy;
	}//loop on hits
        dropenergy+=sumE[droplayer]*weights[droplayer]/tanh(eta_values[eta_counter]);
        for(unsigned ilayer=0; ilayer<layers; ilayer++)
           if(sumE[ilayer]!=0){
               centerx[ilayer]=sumx[ilayer]/sumE[ilayer];
               centery[ilayer]=sumy[ilayer]/sumE[ilayer];
           }

              
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
	  const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  //double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();
          
	  double weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);
          if(layer==droplayer+1) weighted_energy = energy*(weights[droplayer]+weights[droplayer+1])/tanh(eta_values[eta_counter]);

	  if(sumE[layer]>0 && layer!=droplayer && abs(int(posx)/10-(int)(centerx[layer])/10)<=truncationwidth && abs((int)posy/10-(int)(centery[layer])/10)<=truncationwidth) // truncate around energy center 
	    totalE += weighted_energy;
	}//loop on hits

	total_energies.push_back(totalE);
	meanE_anticipated += totalE/nEvts;
	sigE_anticipated += totalE*totalE/nEvts;
      }//loop on evts

      sigE_anticipated = sigE_anticipated - meanE_anticipated*meanE_anticipated;
      sigE_anticipated = sqrt(sigE_anticipated);
      double lower_bound_for_hist;
      double upper_bound_for_hist;
      //      if(energy_incoming_gev < 200) {
      lower_bound_for_hist = meanE_anticipated - 6.5*sigE_anticipated;
      upper_bound_for_hist = meanE_anticipated + 6.5*sigE_anticipated;
	//}
      
      TCanvas *myc = new TCanvas("Energy Distribution","Energy",800,600);
      myc->cd();
      TH1F *p_l = new TH1F("Energy Distribution", "Energies", 50, lower_bound_for_hist,upper_bound_for_hist);
      for(unsigned hist_filler_counter = 0; hist_filler_counter < total_energies.size(); hist_filler_counter++) {
	p_l->Fill(total_energies[hist_filler_counter]);
      }
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
      delete fit;
      delete p_l;
      delete myc;
    }// end for et

  }//end for eta
}

void fixE_eta_resol(){
    const int layers=28;
    const int et[3]={1,6,14};//et=5,40,175
    const int eta[2]={1,2};//eta=2.1
    TCanvas *myc=new TCanvas("myc","myc",1000,1000);
    std::vector<double> resol;
    std::vector<double> resol_err;
    std::vector<double> layer_array;
    std::vector<double> layer_array_err;
    std::vector<double> showerenergy;
    for(int ilayer=1; ilayer<layers; ilayer++){
        std::cout<<"droplayer="<<ilayer<<std::endl;
        double resol_value,resol_err_value,dropenergy;
        plotE_resolution(ilayer,et[2],eta[1],resol_value,resol_err_value,dropenergy);
        layer_array_err.push_back(0);
        layer_array.push_back((double)ilayer);
        resol.push_back(resol_value);
        resol_err.push_back(resol_err_value);
        showerenergy.push_back(dropenergy);
    }
    //////////change heights of showerenergy////////
    double maxresol=0,maxshower=0;
    for(int ilayer=1; ilayer<layers; ilayer++)
       if(resol[ilayer-1]>maxresol) maxresol=resol[ilayer-1];
    for(int ilayer=1; ilayer<layers; ilayer++)      
       if(showerenergy[ilayer-1]>maxshower) maxshower=showerenergy[ilayer-1];
    for(int ilayer=1; ilayer<layers; ilayer++)    
       showerenergy[ilayer-1]/=(maxshower/maxresol)/1.2;
    /////////////////end//////////////////////
    TVectorD Tresol(resol.size(),&resol[0]);
    TVectorD Tresol_err(resol_err.size(),&resol_err[0]);
    TVectorD Tlayer_array(layer_array.size(),&layer_array[0]);
    TVectorD Tlayer_array_err(layer_array_err.size(),&layer_array_err[0]);
    TVectorD Tshowerenergy(showerenergy.size(),&showerenergy[0]);
    TGraph *shower=new TGraph(Tlayer_array,Tshowerenergy);
    TGraphErrors *fixE_resol_droplayer=new TGraphErrors(Tlayer_array,Tresol,Tlayer_array_err,Tresol_err);
    TString name=Form("droplayer_resol_et150eta2.5");
    myc->cd();
    shower->SetLineColor(3);
    shower->SetTitle(name);
    shower->Draw("AC*");
    fixE_resol_droplayer->SetLineColor(4);
    fixE_resol_droplayer->Draw("CP");  
    myc->SaveAs(name+Form(".png"));
}
