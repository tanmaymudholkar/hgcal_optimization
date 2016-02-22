#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TVector.h"

#include "../../userlib/include/HGCSSEvent.hh"
#include "../../userlib/include/HGCSSInfo.hh"
#include "../../userlib/include/HGCSSRecoHit.hh"
#include "../../userlib/include/HGCSSSimHit.hh"
#include "../../userlib/include/HGCSSSamplingSection.hh"

// Double_t radius_fit(Double_t *x, Double_t *y, Double_t *par) {
//   return (par[0]*(exp(-(par[1]*(x[0]-par[4])*(x[0]-par[4])+par[2]*(y[0]-par[5])*(y[0]-par[5])+2*par[3]*(x[0]-par[4])*(y[0]-par[5])))));
// }

void plotE_radius_fit(){//main  

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

  //double chisqr;
  //double ndfr;
  //double chisqdf;

  //double et_values_array[] = {3,5,7,10,30,50,70,100};
  //double et_values_array[] = {3,5,7,10,50,70,100,200};
  //double eta_values_array[] = {1.6,2.1,2.5};

  //std::vector<double> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(double));
  //std::vector<double> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(double));

  TString HGcal_common_prefix = Form("/export/cmss2/paulini/CMS/HGCal/version30/HGcal_version30_model2_BOFF_");
  TString Digi_common_prefix = Form("/export/cmss2/paulini/CMS/HGCal/version30/Digi_version30_model2_BOFF_");
  TString common_suffix = Form(".root");

  //for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
  //TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]);
  double eta = 2.1;
  double et = 30;
  TString eta_portion = Form("_eta%.3f",eta);
  //std::vector<Double_t> energies_incoming_gev;
  //std::vector<Double_t> energies_incoming_gev_errors;
  //std::vector<Double_t> energies;
  //std::vector<Double_t> energies_errors;
  //std::vector<Double_t> resolutions;
  //std::vector<Double_t> resolutions_errors;
  //for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
  //Double_t energy_incoming_gev = et*cosh(eta);
  TString et_portion = Form("et%.0f",et);
  //double lower_bound_for_hist;
  //double upper_bound_for_hist;
  // if(energy_incoming_gev < 25) {
  //   lower_bound_for_hist = 0.6*49.5*energy_incoming_gev;
  //   upper_bound_for_hist = 1.4*49.5*energy_incoming_gev;
  // }
  // else{
  //   lower_bound_for_hist = 0.8*49.5*energy_incoming_gev;
  //   upper_bound_for_hist = 1.2*49.5*energy_incoming_gev;
  // }
  //Double_t energy_incoming_gev_error;
  // Double_t mean_energy;
  // Double_t mean_energy_error;
  // Double_t stddev;
  // Double_t stddev_error;
  // Double_t resolution;
  // Double_t resolution_error;
  //std::cout << "et" << et_values[et_counter] << " eta" << eta_values[eta_counter] << std::endl;

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
  //double totalE(0);
  // TCanvas *myc = new TCanvas("Energy Distribution","Energy",800,600);
  // myc->cd();
  // TH1F *p_l = new TH1F("Energy Distribution","Energies",50,lower_bound_for_hist,upper_bound_for_hist);
  //TH1F *p_l_unweighted = new TH1F("Layers unw","layer",28,-0.5,27.5);
  //TH2F *p_E_l = new TH2F("total E and layer number in ECAL","layer number",500,0,1500,28,-0.5,27.5);

  //std::vector<double> xrms;
  //std::vector<double> yrms;
  //Double_t xrms[28];
  //Double_t yrms[28];
  Double_t radius[28];
  std::vector<Double_t> radii[28];
  Double_t radius_error[28];
  Double_t layer_number[28];
  Double_t layer_number_error[28];
  //int n_of_measurements_of_radii_in_layer[28];

  // for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
  //   xrms[layer_counter] = 0;
  //   yrms[layer_counter] = 0;
  //   how_many_to_count[layer_counter] = 0;
  // }

  for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
    radius[layer_counter] = 0;
    layer_number[layer_counter] = layer_counter + 1.0;
    layer_number_error[layer_counter] = 0;
    //n_of_measurements_of_radii_in_layer[layer_counter] = 0;
  }


  //unsigned ievt=1;
  
  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
  //totalE = 0;
    
    if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);
    
    std::vector<Double_t> xpos[28];
    std::vector<Double_t> ypos[28];
    std::vector<Double_t> E_of_hit[28];
    Double_t totE_in_layer_for_one_evt[28];
    for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
      totE_in_layer_for_one_evt[layer_counter] = 0;
    }
  
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
      const HGCSSRecoHit lHit = (*rechitvec)[iH];
      
      double posx = lHit.get_x();
      double posy = lHit.get_y();
      //double posz = lHit.get_z();
      unsigned layer = lHit.layer();
      double energy = lHit.energy();
      
      double weighted_energy = energy*weights[layer]/tanh(eta);
      
      //totalE += weighted_energy;
      xpos[layer].push_back(posx);
      ypos[layer].push_back(posy);
      E_of_hit[layer].push_back(weighted_energy);
      totE_in_layer_for_one_evt[layer] += weighted_energy;
      //p_l_unweighted->Fill(layer);
    }
    
    for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
      if(totE_in_layer_for_one_evt[layer_counter] > 0) {
	//n_of_measurements_of_radii_in_layer[layer_counter]++;
	std::vector<Double_t> weighted_xpos[28];
	std::vector<Double_t> weighted_ypos[28];
	for(unsigned hit_counter = 0; hit_counter != xpos[layer_counter].size(); hit_counter++) {
	  weighted_xpos[layer_counter].push_back(xpos[layer_counter][hit_counter]*E_of_hit[layer_counter][hit_counter]/totE_in_layer_for_one_evt[layer_counter]);
	  weighted_ypos[layer_counter].push_back(ypos[layer_counter][hit_counter]*E_of_hit[layer_counter][hit_counter]/totE_in_layer_for_one_evt[layer_counter]);
	}
	//TVectorD Tweighted_xpos(weighted_xpos[layer_counter].size(),&weighted_xpos[layer_counter][0]);
	//TVectorD Tweighted_ypos(weighted_ypos[layer_counter].size(),&weighted_ypos[layer_counter][0]);
	
	//TGraph *weighted_positions = new TGraph(Tweighted_xpos,Tweighted_ypos);
	//double meanx = weighted_positions->GetMean(1);
	//double meany = weighted_positions->GetMean(2);
	double meanx = 0;
	double meany = 0;
	for(unsigned mean_counter = 0; mean_counter < weighted_xpos[layer_counter].size(); mean_counter++) {
	  meanx += weighted_xpos[layer_counter][mean_counter];
	  meany += weighted_ypos[layer_counter][mean_counter];
	}
	meanx = meanx/weighted_xpos[layer_counter].size();
	meany = meany/weighted_ypos[layer_counter].size();	
	double stdevx = 0;
	double stdevy = 0;
	for(unsigned stdev_counter = 0; stdev_counter < weighted_xpos[layer_counter].size(); stdev_counter++) {
	  stdevx += (weighted_xpos[layer_counter][stdev_counter]-meanx)*(weighted_xpos[layer_counter][stdev_counter]-meanx);
	  stdevy += (weighted_ypos[layer_counter][stdev_counter]-meany)*(weighted_ypos[layer_counter][stdev_counter]-meany);
	}
	stdevx = sqrt(stdevx/weighted_xpos[layer_counter].size());
	stdevy = sqrt(stdevy/weighted_ypos[layer_counter].size());
	radii[layer_counter].push_back(0.5*(stdevx+stdevy));
      }
    }
    
    
    // for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
    //   if(n_evts_in_layer[layer_counter]>8) {
    // 	TGraph2D *energies = new TGraph2D(xpos[layer_counter],ypos[layer_counter],E_of_hit[layer_counter]);
    // 	TH2D *energies_histogram = energies->GetHistogram();
    // 	TF2 *gaussian_to_fit = new TF2("g_to_fit",radius_fit,-1500,1500,-1500,1500,6);
    // 	double expected_xcenter = energies_histogram->GetMean(1);
    // 	double expected_ycenter = energies_histogram->GetMean(2);
    // 	double expected_sigx = energies_histogram->GetRMS(1);
    // 	double expected_sigy = energies_histogram->GetRMS(2);
    // 	gaussian_to_fit->SetParameters(energies_histogram->GetZmax()/(expected_sigx*expected_sigy),1/sqrt(expected_sigx),1/sqrt(expected_sigy),0,expected_xcenter,expected_ycenter);
    // 	energies_histogram->Fit("g_to_fit");
    // 	double a = gaussian_to_fit->GetParameter(1);
    // 	double b = gaussian_to_fit->GetParameter(2);
    // 	double c = gaussian_to_fit->GetParameter(3);
    // 	if(a*b>c*c) {
    // 	  n_of_fits[layer_counter]++;
    // 	  radius[layer_counter] += (0.25/sqrt(a*b-c*c))*(sqrt(a+b+sqrt((a-b)*(a-b)+4*c*c))+sqrt(a+b-sqrt((a-b)*(a-b)+4*c*c)));
    // 	}
    //   }
    // }
    //p_l->Fill(totalE);
  }//loop on hits
  
  // for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
  //   if (n_of_fits[layer_counter]>0) {
  //     radius[layer_counter] = radius[layer_counter]/n_of_fits[layer_counter];
  //     layer_number[layer_counter] = layer_counter + 1.0;
  //   }
  // }
  
  // for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
  //   if(n_of_measurements_of_radii_in_layer[layer_counter]>0) {
  //     radius[layer_counter] = radius[layer_counter]/n_of_measurements_of_radii_in_layer[layer_counter];
  //   }
  // }
  
  for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
    std::cout << "layer = " << 1+layer_counter << std::endl;
    //TVectorD Tradii(radii[layer_counter].size(),&radii[layer_counter][0]);
    //TGraph *gradii = new TGraph(Tradii,Tradii);
    //radius[layer_counter] = gradii->GetMean(1);
    for(unsigned mean_counter = 0; mean_counter < radii[layer_counter].size(); mean_counter++) {
      radius[layer_counter] += radii[layer_counter][mean_counter];
      std::cout << "radius is " << radii[layer_counter][mean_counter] << std::endl;
    }
    radius[layer_counter]=radius[layer_counter]/radii[layer_counter].size();
    double stdevradii = 0;
    for(unsigned stdev_counter = 0; stdev_counter < radii[layer_counter].size(); stdev_counter++) {
      stdevradii += (radii[layer_counter][stdev_counter]-radius[layer_counter])*(radii[layer_counter][stdev_counter]-radius[layer_counter]);
    }
    stdevradii = sqrt(stdevradii/radii[layer_counter].size());
    radius_error[layer_counter] = stdevradii;
  }
  
  for(unsigned layer_counter = 0; layer_counter != 28; layer_counter++) {
    //std::cout << layer_number[layer_counter] << "   " << radius[layer_counter] << "   " << n_of_measurements_of_radii_in_layer[layer_counter] << "   " << totE_in_layer_for_one_evt[layer_counter] << std::endl;
    std::cout << layer_number[layer_counter] << "   " << radius[layer_counter] << "   +/-   " << radius_error[layer_counter]  << std::endl;
  }
  
  
  TCanvas *myc = new TCanvas("radii","radius",800,600);
  myc->cd();
  
  //TGraphErrors *gradii_final = new TGraphErrors(28,layer_number,layer_number_error,radius,radius_error);
  TGraph *gradii_final = new TGraph(28,layer_number,radius);
  gradii_final->Draw();
  myc->Update();
  myc->SaveAs("test_radii.pdf");
  
//   //p_l->Fit("gaus");
//   TF1 *fit = p_l->GetFunction("gaus");
//   mean_energy = fit->GetParameter(1);
//   mean_energy_error = fit->GetParError(1);
//   stddev = fit->GetParameter(2);
//   stddev_error = fit->GetParError(2);
//   chisqr = fit->GetChisquare();
//   ndfr = fit->GetNDF();
//   std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
//   chisqdf = chisqr/ndfr;
//   resolution = stddev/mean_energy;
//   resolution_error = resolution*((stddev_error/stddev) + (mean_energy_error/mean_energy));
//   energy_incoming_gev_error = energy_incoming_gev*mean_energy_error/mean_energy;
  
//   energies_incoming_gev.push_back(energy_incoming_gev);
//   energies_incoming_gev_errors.push_back(energy_incoming_gev_error);
//   energies.push_back(mean_energy);
//   energies_errors.push_back(mean_energy_error);
//   resolutions.push_back(resolution);
//   resolutions_errors.push_back(resolution_error);
  
//   //Double_t tavg_energy = avg_energy;

//   //TString s_avg_energy;
  
//   //s_avg_energy.Form("%f\n",tavg_energy);
  
//   //std::ostringstream ss;
//   //ss << avg_energy;
//   //std::string title_to_set = "avg_energy = " + ss.str();
  
//   //p_l->SetMaximum(6000000);
//   p_l->SetTitle(et_portion + eta_portion + Form(". chisq/dof = %3.2f",chisqdf));
//   p_l->Draw();
//   //p_l_unweighted->Draw();
//   //p_E->Draw();
//   myc->SaveAs(Form("plot_resolutions_")+et_portion+eta_portion+Form(".png"));
//   //cout << "min layer = " << min_layer << std::endl;
//   //cout << "max layer = " << max_layer << std::endl;
      
// }
  
// TVectorD Tenergies(energies.size(),&energies[0]);
// TVectorD Tresolutions(resolutions.size(),&resolutions[0]);
// TVectorD Tenergies_incoming_gev(energies_incoming_gev.size(),&energies_incoming_gev[0]);
// TVectorD Tenergies_incoming_gev_errors(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
// TVectorD Tenergies_errors(energies_errors.size(),&energies_errors[0]);
// TVectorD Tresolutions_errors(resolutions_errors.size(),&resolutions_errors[0]);
// TCanvas *myc = new TCanvas("Energy mips calibration", "mips versus gev",800,600);
// myc->cd();
// TGraphErrors *energy_mips_calibration = new TGraphErrors(Tenergies_incoming_gev,Tenergies,Tenergies_incoming_gev_errors,Tenergies_errors);
// energy_mips_calibration->Fit("pol1");
// TF1 *calibration_fit = energy_mips_calibration->GetFunction("pol1");
// chisqr = calibration_fit->GetChisquare();
// ndfr = calibration_fit->GetNDF();
// std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
// chisqdf = chisqr/ndfr;
// double conversion_factor = calibration_fit->GetParameter(1);
// double conversion_factor_error = calibration_fit->GetParError(1);
// TString start_title_conversion_factor = Form("Best Fit: E_mips = (");
// energy_mips_calibration->SetTitle(start_title_conversion_factor + Form("%5.3f +/- %5.3f) * E_GeV. chisq/dof = %3.2f",conversion_factor,conversion_factor_error,chisqdf));
// energy_mips_calibration->Draw("AP");
// myc->Update();
// myc->Print(Form("plot_calibration_fit") + eta_portion + Form(".pdf"));
// for (unsigned int energies_incoming_gev_errors_counter = 0; energies_incoming_gev_errors_counter != energies_incoming_gev_errors.size(); energies_incoming_gev_errors_counter++) {
//   energies_incoming_gev_errors[energies_incoming_gev_errors_counter] += energies_incoming_gev[energies_incoming_gev_errors_counter]*conversion_factor_error/conversion_factor;
//  }
// TVectorD Tenergies_incoming_gev_errors_updated(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
// myc = new TCanvas("Resolutions versus Energy","Resolution versus Energy",800,600);
// myc->cd();
// TGraphErrors *resolutions_versus_energy = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors_updated,Tresolutions_errors);
// TF1 *function_to_fit = new TF1("f_to_fit",resolutions_fit,0.0001,1000,2);
// function_to_fit->SetParameters(0.015,0.1);
// function_to_fit->SetParNames("const_offset","coeff_sqrt");
// resolutions_versus_energy->Fit("f_to_fit");
// chisqr = function_to_fit->GetChisquare();
// ndfr = function_to_fit->GetNDF();
// std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
// chisqdf = chisqr/ndfr;
// TString start_title_resolutions = Form("Best Fit: resolution = (");
// double bestfit_constant = function_to_fit->GetParameter(0);
// double bestfit_constant_error = function_to_fit->GetParError(0);
// double bestfit_coeff = function_to_fit->GetParameter(1);
// double bestfit_coeff_error = function_to_fit->GetParError(1);
// resolutions_versus_energy->SetTitle(start_title_resolutions + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV). chisq/dof = %3.2f",bestfit_constant,bestfit_constant_error,bestfit_coeff,bestfit_coeff_error,chisqdf));
// resolutions_versus_energy->Draw("AP");
// myc->Update();
// myc->Print(Form("plot_resolutions_versus_energy") + eta_portion + Form(".pdf"));
// // std::cout << std::endl;
    
// // std::cout << "____________________________________________________________________________________" << std::endl;
// // std::cout << "eta = " << eta_values[eta_counter] << std::endl;
// // std::cout << std::endl;
  
// // for(unsigned int resolutions_counter = 0; resolutions_counter != resolutions.size(); resolutions_counter++) {
// //   std::cout << energies[resolutions_counter] << "    " << resolutions[resolutions_counter] << std::endl;
// // }
// }
}
