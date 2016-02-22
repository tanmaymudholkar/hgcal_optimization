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
#include "../../userlib/include/HGCSSGenParticle.hh"

Double_t cellSize;
const Double_t radlim = 750;
const Double_t threshold = 0.5;
const Double_t sigmas_down = 3.0;
const Double_t sigmas_up = 3.0;
const Int_t signal_region = 1000;

// Double_t resolutions_fit(Double_t *x,Double_t *par) {
//   return (par[1] + par[0]/sqrt(x[0]));
// }

// Double_t resolutions_fit_standalone_withnoise(Double_t *x,Double_t *par) {
//   return sqrt(pow(par[1],2) + pow(par[0],2)/x[0] + pow(par[2]/x[0],2));
//   // return sqrt(pow(par[1],2) + pow(par[0],2)/x[0]);
// }

// Double_t resolutions_fit_standalone(Double_t *x,Double_t *par) {
//   // return sqrt(pow(par[1],2) + pow(par[0],2)/x[0] + pow(par[2]/x[0],2));
//   return sqrt(pow(par[1],2) + pow(par[0],2)/x[0]);
// }

// Double_t resolutions_fit_standalone_const(Double_t *x,Double_t *par) {
//   return sqrt(pow(0.009,2) + pow(par[0],2)/x[0] + pow(par[1]/x[0],2));
//   // return sqrt(pow(0.009,2) + pow(par[0],2)/x[0]);
// }

Double_t cellsizeat(Double_t x, Double_t y) {
  Double_t r = sqrt(pow(x,2)+pow(y,2));
  if (r<radlim) {
    return (3.0*cellSize);
  }
  else {
    return (4.0*cellSize);
  }
}

void plotE_resolution_high_stats_syst_opt(Int_t version_number, TString version_name, TString datadir, Double_t et) { // main
  
  // load the shared library for HGCSS* classes:
  gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

  // if(threshold<0.5) {
  //   std::cout << "Threshold provided less than minimum threshold analyzable from Digi file" << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }
    
  Double_t chisqr;
  Double_t ndfr;
  Double_t chisqdf;

  // Double_t et_values_array[] = {5,50}; //For systematic optimization studies
  // Double_t et_values_array[] = {3,5,7,10,20,30,40,50,60,70,100,125,150}; //Intersection version 30 & 34
  // Double_t et_values_array[] = {3,5,7,10,20};
  // Double_t et_values_array[] = {3,5,10,20,30,40,50,60,80,100,125,175};
  // Double_t et_values_array[] = {3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  // Double_t et_values_array[] = {3,7,20,35,50,70,100,125};
  // Double_t et_values_array[] = {100,125,150,175,200};
  // Double_t et_values_array[] = {3};
  Double_t eta_values_array[] = {2.1};
  
  // std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
  std::vector<Double_t> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(Double_t));
  
  // TString str_threshold = Form(".6f",threshold);
  TString HGcal_common_prefix_firstpart = datadir + Form("/HGcal__version%i_model2_BOFF_",version_number);
  TString Digi_common_prefix_firstpart = datadir + Form("/Digi__version%i_model2_BOFF_",version_number);
  TString common_suffix = Form(".root");

  ofstream outfile;
  TString outfile_name_prefix = Form("data_et%.0f_",et);
  outfile.open(outfile_name_prefix+version_name+Form("_resolutions"));

  for (Double_t rwcuf = 0; rwcuf<=1.0; rwcuf += 0.1) {
    for (Double_t rwcum = 0; rwcum<=1.0; rwcum += 0.1) {
      if (rwcuf < 0.95 || rwcum < 0.95) {

	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "rwcuf =  " << rwcuf << "; rwcum = " << rwcum << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	std::cout << "___________________________________________________________________________" << std::endl;
	
	std::vector<Double_t> weights;
	ifstream f_layer_weights;
	f_layer_weights.open(datadir+Form("/layer_weights_rwcuf%.1f_rwcum%.1f.dat",rwcuf,rwcum));
	std::string line;
	Double_t weight;
	if (f_layer_weights.is_open()) {
	  while (getline(f_layer_weights,line)) {
	    weight = std::atof(line.c_str());
	    weights.push_back(weight);
	  }
	}

	TString HGcal_common_prefix = HGcal_common_prefix_firstpart + Form("rwcuf_%.1f_rwcum_%.1f_",rwcuf,rwcum);
	TString Digi_common_prefix = Digi_common_prefix_firstpart + Form("rwcuf_%.1f_rwcum_%.1f_",rwcuf,rwcum);
	
	// std::cout << HGcal_common_prefix << std::endl;
	// std::cout << Digi_common_prefix << std::endl;
  
	for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
	  TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]);
	  // std::vector<Double_t> energies_incoming_gev;
	  // std::vector<Double_t> energies_incoming_gev_errors;
	  // std::vector<Double_t> mean_energies_wmips;
	  // std::vector<Double_t> mean_energies_wmips_errors;
	  // std::vector<Double_t> sigmas_wmips;
	  // std::vector<Double_t> sigmas_wmips_errors;
	  // std::vector<Double_t> resolutions;
	  // std::vector<Double_t> resolutions_errors;
	  Double_t resolution;
	  Double_t resolution_error;
	  // for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
	  // Double_t energy_incoming_gev = et*cosh(eta_values[eta_counter]);
	  TString et_portion = Form("et%.0f",et);
	  // Double_t energy_incoming_gev_error;
	  Double_t mean_energy_wmips;
	  Double_t mean_energy_wmips_error;
	  Double_t sigma_wmips;
	  Double_t sigma_wmips_error;
	  std::vector<Double_t> total_energies;
	  unsigned nEvts_to_count = 0;
	  Double_t meanE_anticipated = 0;
	  Double_t sigE_anticipated = 0;
	  std::cout << "et" << et << " eta" << eta_values[eta_counter] << std::endl;
	  
	  for (unsigned run_no = 1; run_no <= 10; run_no++) {
	    std::cout << "___________________________________________________________________________" << std::endl;
	    std::cout << "run number " << run_no << std::endl;
	    std::cout << "___________________________________________________________________________" << std::endl;
	    TChain  *lSimTree = new TChain("HGCSSTree");
	    TChain  *lRecTree = new TChain("RecoTree");
	  
	    lSimTree->AddFile(HGcal_common_prefix+Form("runno_%i_",run_no)+et_portion+eta_portion+common_suffix);
	    lRecTree->AddFile(Digi_common_prefix+Form("runno_%i_",run_no)+et_portion+eta_portion+common_suffix);
	  
	    TFile *inputFile = TFile::Open(HGcal_common_prefix+Form("runno_%i_",run_no)+et_portion+eta_portion+common_suffix);
	    HGCSSInfo * info=(HGCSSInfo*)inputFile->Get("Info");
	    cellSize = info->cellSize();
	    const unsigned versionNumber = info->version();
	    const unsigned model = info->model();
	    delete inputFile;
	    
	    std::cout << "cellSize is " << cellSize << std::endl;
	    std::cout << "versionNumber is " << versionNumber << std::endl;
	    std::cout << "model is " << model << std::endl;
	    
	    std::vector<HGCSSSamplingSection> * ssvec = 0;
	    std::vector<HGCSSSimHit> * simhitvec = 0;
	    std::vector<HGCSSRecoHit> * rechitvec = 0;
	    std::vector<HGCSSGenParticle> * genvec = 0; // get the generated particle vector
      
	    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
	    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
	    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
	    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
      
	    const unsigned nEvts = lSimTree->GetEntries();
	    Double_t totalE(0);
      
	    for (unsigned ievt(0); ievt<nEvts; ++ievt){// loop on entries
	      totalE = 0;
	
	      if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	      lSimTree->GetEntry(ievt);
	      lRecTree->GetEntry(ievt);

	      //std::cout << "track id " << (*genvec)[0].trackID() << "  genvec size " << (*genvec).size() << std::endl;

	      // if((*genvec).size() > 1) continue;  // Delete events which has converted photons.
	      if((*genvec)[0].trackID() >= 2) continue;
	      else{
		nEvts_to_count += 1;
		const HGCSSGenParticle gHit = (*genvec)[0];
		//std::cout << "ACCEPTED: " << "track id " << (*genvec)[0].trackID() << "  genvec size " << (*genvec).size() << std::endl;
		Double_t posx_gen_ini = gHit.x();
		Double_t posy_gen_ini = gHit.y();
		Double_t posz_gen_ini = gHit.z();
		Double_t px_gen = gHit.px();
		Double_t py_gen = gHit.py();
		Double_t pz_gen = gHit.pz();
	
		// for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){// loop over simhits
		for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
		  // const HGCSSSimHit lHit = (*simhitvec)[iH];
		  const HGCSSRecoHit lHit = (*rechitvec)[iH];
		  
		  Double_t posx = lHit.get_x();
		  Double_t posy = lHit.get_y();
		  Double_t posz = lHit.get_z();
		  unsigned layer = lHit.layer();
		  Double_t energy = lHit.energy();
		  // std::cout << "energy is " << energy << "   " << "layer is " << layer << std::endl;
	  
		  if(energy>threshold) {
		    Double_t posx_gen = posx_gen_ini + (px_gen/pz_gen)*(posz - posz_gen_ini);
		    Double_t posy_gen = posy_gen_ini + (py_gen/pz_gen)*(posz - posz_gen_ini);
		    Double_t dx = posx - posx_gen;
		    Double_t dy = posy - posy_gen;
		    Double_t halfCell = 0.5*cellsizeat(posx,posy);
		    if ((fabs(dx) <= signal_region*halfCell) && (fabs(dy) <= signal_region*halfCell)){
		      Double_t weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);
		      // totalE_in_layer[layer] += weighted_energy;
		      totalE += weighted_energy;
		    }// end if condition for counting energies in a signal region
		  }//end if condition for counting energies above a threshold
		}// end loop over rechits
		total_energies.push_back(totalE);
		// meanE_anticipated += totalE/nEvts;
		meanE_anticipated += totalE;
		// sigE_anticipated += totalE*totalE/nEvts;
		sigE_anticipated += totalE*totalE;
	      }// end else condition for counting only meaningful genvecs
	    }// loop on events

	    delete lSimTree;
	    delete lRecTree;
	  }// ends loop over runs

	  meanE_anticipated = meanE_anticipated/nEvts_to_count;
	  sigE_anticipated = sigE_anticipated/nEvts_to_count;
      
	  sigE_anticipated = sigE_anticipated - meanE_anticipated*meanE_anticipated;
	  sigE_anticipated = sqrt(sigE_anticipated);

	  // std::cout << "mean_anticipated is " << meanE_anticipated << std::endl;
	  // std::cout << "sigma_anticipated is " << sigE_anticipated << std::endl;

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
	  p_l->SetTitle(version_name + Form("_rwcuf%.1f_rwcum%.1f . chisq/dof = %3.2f",rwcuf,rwcum,chisqdf));
	  p_l->Draw();
	  myc->Print(Form("plot_rwcuf%.1f_rwcum%.1f_",rwcuf,rwcum)+version_name+Form("_distribution_")+et_portion+eta_portion+Form(".png"));
	  myc->Print(Form("plot_rwcuf%.1f_rwcum%.1f_",rwcuf,rwcum)+version_name+Form("_distribution_")+et_portion+eta_portion+Form(".pdf"));
	  delete gaussian_fit;
	  delete p_l;
	  delete myc;
	  
	  resolution = sigma_wmips/mean_energy_wmips;
	  resolution_error = resolution*sqrt(pow(sigma_wmips_error/sigma_wmips,2)+pow(mean_energy_wmips_error/mean_energy_wmips,2));

	  outfile << rwcuf << "    " << rwcum << "    " << resolution << "    " << resolution_error << std::endl;
	} // ends loop over eta
      } // ends if condition taking care of the extreme corner
    } // ends loop over rwcum
  } // ends loop over rwcuf
  outfile.close();
} // ends main
