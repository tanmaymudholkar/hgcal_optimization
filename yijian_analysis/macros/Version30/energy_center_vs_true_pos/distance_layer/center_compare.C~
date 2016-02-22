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
#include "TGraph.h"
#include "TGraph2D.h"

#include "../../../../userlib/include/HGCSSEvent.hh"
#include "../../../../userlib/include/HGCSSInfo.hh"
#include "../../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../../userlib/include/HGCSSSimHit.hh"
#include "../../../../userlib/include/HGCSSSamplingSection.hh"
#include "../../../../userlib/include/HGCSSGenParticle.hh"

void center_compare(){//main  
   const int layers=28; //number of layers, varied from version to version but remaining constant in one version
      TH1F *distance_evts[layers];
      TString THname[layers];
      for(int l=0;l<layers;l++)
          THname[l]=strcat(Form("h"),Form("%d",l));
      for(int l=0;l<layers;l++)
          distance_evts[l]=new TH1F(THname[l],THname[l],20,-0.5,39.5);
   double layer_array[layers];
   for(unsigned layernum=0; layernum<layers; layernum++){
      layer_array[layernum] = (double)layernum;
   }
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

  //for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
  for (unsigned int eta_counter = 2; eta_counter<=2; eta_counter++) {

    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]); 

    //for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
      for (unsigned int et_counter = 5; et_counter <=5 ; et_counter++) {

      TString et_portion = Form("et%.0f",et_values[et_counter]); // a string "et30" ....

      std::cout << "et" << et_values[et_counter] << " eta" << eta_values[eta_counter] << std::endl; // printf: et, eta

      TChain  *lSimTree = new TChain("HGCSSTree");
      TChain  *lRecTree = new TChain("RecoTree");

      lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix); // simulated events
      lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);  // reconstructed events
 
      std::vector<HGCSSSamplingSection> * ssvec = 0;
      std::vector<HGCSSSimHit> * simhitvec = 0;
      std::vector<HGCSSRecoHit> * rechitvec = 0;
      std::vector<HGCSSGenParticle> * genvec = 0; // get the generated particle vector 

      lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
      lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
      lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
      lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec); 

      const unsigned nEvts = lSimTree->GetEntries(); 

      for (unsigned ievt(0); ievt<nEvts; ievt+=1){//loop on events

	if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	lSimTree->GetEntry(ievt);
	lRecTree->GetEntry(ievt);
        
        double centerx[layers]={0},centery[layers]={0},centerz[layers]={0},sumx[layers]={0},sumy[layers]={0},sumz[layers]={0},sumE[layers]={0};

        for(unsigned iHt(0); iHt<(*rechitvec).size(); ++iHt){ 
          const HGCSSRecoHit lHit = (*rechitvec)[iHt];
    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();
          sumx[layer]+=posx*energy;
          sumy[layer]+=posy*energy;
          sumz[layer]+=posz*energy;
          sumE[layer]+=energy;
        }//end for iHt

        for(unsigned ilayer(0); ilayer<layers; ilayer++){
          if(sumE[ilayer]>0){
            centerx[ilayer]=sumx[ilayer]/sumE[ilayer];
            centery[ilayer]=sumy[ilayer]/sumE[ilayer];
            centerz[ilayer]=sumz[ilayer]/sumE[ilayer];
            }
        }//end for ilayer

        double gen_x[layers],gen_y[layers],gen_posx0,gen_posy0,gen_posz0,gen_px,gen_py,gen_pz;

        if(1||!((*genvec).size() > 1)){  // Delete events which has converted photons.
           const HGCSSGenParticle gHit = (*genvec)[0];
    
	   gen_posx0 = gHit.x(); //initial position 
	   gen_posy0 = gHit.y();
           gen_posz0 = gHit.z();
           gen_px=gHit.px();
           gen_py=gHit.py();
           gen_pz=gHit.pz();
        }
        
        //calculate true position

        for(unsigned ilayer=0; ilayer<layers; ilayer++){
          if(sumE[ilayer]>0){
            gen_x[ilayer]=gen_posx0+gen_px/gen_pz*(centerz[ilayer]-gen_posz0);
            gen_y[ilayer]=gen_posy0+gen_py/gen_pz*(centerz[ilayer]-gen_posz0);
            double d=sqrt((centerx[ilayer]-gen_x[ilayer])*(centerx[ilayer]-gen_x[ilayer])+(centery[ilayer]-gen_y[ilayer])*(centery[ilayer]-gen_y[ilayer]));
            distance_evts[ilayer]->Fill(d);
          }
        }
      }//loop on events

    }// end for et
                    
   }// end for eta
   for(int ilayer=0;ilayer<layers;ilayer++){
       TCanvas *myc=new TCanvas("myc","myc",800,600);
       myc->cd();
       distance_evts[ilayer]->Draw("COLZ");
       myc->SaveAs((THname[ilayer]+Form(".png")));
   }
}
  
