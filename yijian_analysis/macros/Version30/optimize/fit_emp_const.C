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

#include "/export/home/yijianz/research/HGCstandalone/userlib/include/HGCSSEvent.hh"
#include "/export/home/yijianz/research/HGCstandalone/userlib/include/HGCSSInfo.hh"
#include "/export/home/yijianz/research/HGCstandalone/userlib/include/HGCSSRecoHit.hh"
#include "/export/home/yijianz/research/HGCstandalone/userlib/include/HGCSSSimHit.hh"
#include "/export/home/yijianz/research/HGCstandalone/userlib/include/HGCSSSamplingSection.hh"
#include "/export/home/yijianz/research/HGCstandalone/userlib/include/HGCSSGenParticle.hh"

void fit_emp_const(){ //read data

    const int layers=24;
    const int et_counter=1; //et=5
    const int eta_counter=1; //eta=2.1
    gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");

    std::vector<double> weights;
    ifstream f_layer_weights;
    f_layer_weights.open("weights_v34.dat");
    std::string line;
    double weight;
    if (f_layer_weights.is_open()) {
    while (getline(f_layer_weights,line)) {
      weight = std::atof(line.c_str());
      weights.push_back(weight);
      }
    }

  double et_values_array[] = {3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200};
  double eta_values_array[] = {1.6,2.1,2.5};

  std::vector<double> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(double));
  std::vector<double> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(double));

  TString HGcal_common_prefix = Form("/export/cmss2/tmudholk/HGCal/version34/HGcal__version34_model2_BOFF_");// Form: convert a string to Tstring
  TString Digi_common_prefix = Form("/export/cmss2/tmudholk/HGCal/version34/Digi__version34_model2_BOFF_");
  TString common_suffix = Form(".root");

      TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]);
      Double_t energy_incoming_gev = et_values[et_counter]*cosh(eta_values[eta_counter]);
      TString et_portion = Form("et%.0f",et_values[et_counter]);
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
      
      double totalE[layers]={0};

      for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

	if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	lSimTree->GetEntry(ievt);
	lRecTree->GetEntry(ievt);  
        const unsigned nEvts = lSimTree->GetEntries(); 
	for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
	  const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
	  // double posx = lHit.get_x();
	  // double posy = lHit.get_y();
	  //double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();
	  //	  std::cout << "energy is " << energy << "   " << "layer is " << layer << std::endl;

	  double weighted_energy = energy*weights[layer]/tanh(eta_values[eta_counter]);
	  
	    totalE[layer] += weighted_energy;
        }//end for hits
     }//end for evts
     
     for (unsigned ilayer(0); ilayer<layers; ilayer++) // average mips
         totalE[ilayer]/=nEvts;
     double avermip=0;
     for (unsigned ilayer(0); ilayer<layers; ilayer++)
         avermip+=totalE[ilayer];
     double calibration_constant=energy_incoming_gev/avermip;
     for (unsigned ilayer(0); ilayer<layers; ilayer++)
         totalE[ilayer]*=calibration_constant;
     TH1F *myt=new TH1F("el","el",24,0,24);
     for (unsigned ilayer(0); ilayer<layers; ilayer++)
         myt->Fill(ilayer,totalE[ilayer]);
     myt->Draw("COLZ");
}//end for readdata     
