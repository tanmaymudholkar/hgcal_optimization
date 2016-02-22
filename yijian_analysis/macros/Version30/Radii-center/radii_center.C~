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

#include "../../../userlib/include/HGCSSEvent.hh"
#include "../../../userlib/include/HGCSSInfo.hh"
#include "../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../userlib/include/HGCSSSimHit.hh"
#include "../../../userlib/include/HGCSSSamplingSection.hh"
#include "../../../userlib/include/HGCSSGenParticle.hh"


typedef struct{
   int layernum;
   double x;
   double y;
   double E;
   double E_suminside; //sum over energy of hits inside the circle passing by the hit.
   double distance; //between hit pos and E center
}onehit;

typedef std::vector<onehit> layerhits;

void sortlayer(layerhits &simlayerhits){ //bubbling sort by distance
   int length=simlayerhits.size();
   if(length==0||length==1) return;
   onehit temp;// used to exchange two hits
   for(unsigned bubbletime=1;bubbletime<length;bubbletime++)
     for(unsigned count=0;count<length-bubbletime;count++)
       if(simlayerhits[count].distance>simlayerhits[count+1].distance){ //exchange
           temp=simlayerhits[count];
           simlayerhits[count]=simlayerhits[count+1];
           simlayerhits[count+1]=temp;
       }
}
       
void radii_center(){//main  
   const int layers=28; //number of layers, varied from version to version but remaining constant in one version

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

  TString HGcal_common_prefix = Form("/export/cmss2/tmudholk/HGCal/version30/HGcal__version30_model2_BOFF_");// Form: convert a string to Tstring
  TString Digi_common_prefix = Form("/export/cmss2/tmudholk/HGCal/version30/Digi__version30_model2_BOFF_");
  TString common_suffix = Form(".root");

  for (unsigned int eta_counter = 0; eta_counter != eta_values.size(); eta_counter++) {
  //for (unsigned int eta_counter = 1; eta_counter<=1; eta_counter++) {

    TString eta_portion = Form("_eta%.3f",eta_values[eta_counter]); 

      for (unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
      //for (unsigned int et_counter = 5; et_counter <=5 ; et_counter++) {

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

      double radius68[layers]={0},radius90[layers]={0};
      int nohit[layers]={0}; //number of events that contain no hit in one layer

      for (unsigned ievt(0); ievt<nEvts; ievt+=1){//loop on events


	if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	lSimTree->GetEntry(ievt);
	lRecTree->GetEntry(ievt);
        
        double centerx[layers]={0},centery[layers]={0},sumx[layers]={0},sumy[layers]={0},sumE[layers]={0};

        for(unsigned iHt(0); iHt<(*rechitvec).size(); ++iHt){ 
          const HGCSSRecoHit lHit = (*rechitvec)[iHt];
    
	  double posx = lHit.get_x();
	  double posy = lHit.get_y();
	  double posz = lHit.get_z();
	  unsigned layer = lHit.layer();
	  double energy = lHit.energy();
          sumx[layer]+=posx*energy;
          sumy[layer]+=posy*energy;
          sumE[layer]+=energy;
        }//end for iHt

        for(unsigned ilayer(0); ilayer<layers; ilayer++){
          if(sumE[ilayer]>0){
            centerx[ilayer]=sumx[ilayer]/sumE[ilayer];
            centery[ilayer]=sumy[ilayer]/sumE[ilayer];
            }
        }//end for ilayer

        layerhits hitdata[layers];
        for(unsigned iHt(0); iHt<(*rechitvec).size(); ++iHt){ 
          const HGCSSRecoHit lHit = (*rechitvec)[iHt];
          onehit simhit;
          simhit.x = lHit.get_x();
	  simhit.y = lHit.get_y();
	  simhit.layernum = lHit.layer();
	  simhit.E = lHit.energy();
          simhit.distance=(simhit.x-centerx[simhit.layernum])*(simhit.x-centerx[simhit.layernum])+(simhit.y-centery[simhit.layernum])*(simhit.y-centery[simhit.layernum]);
          hitdata[simhit.layernum].push_back(simhit);
        }

        for(unsigned layer_num=0;layer_num<layers;layer_num++)
          sortlayer(hitdata[layer_num]);

        for(unsigned layer_num=0;layer_num<layers;layer_num++)
          if(hitdata[layer_num].size()>0){
             (hitdata[layer_num])[0].E_suminside=(hitdata[layer_num])[0].E;
             unsigned hitnum=1;
             while(hitnum<hitdata[layer_num].size()){
                (hitdata[layer_num])[hitnum].E_suminside=(hitdata[layer_num])[hitnum-1].E_suminside+(hitdata[layer_num])[hitnum].E;
                hitnum++;
             }//end while
          }//endif end for layernum

        for(unsigned layer_num=0;layer_num<layers;layer_num++){
          if(hitdata[layer_num].size()==0)
            nohit[layer_num]++;
          else{
            double totalE=(hitdata[layer_num])[hitdata[layer_num].size()-1].E_suminside;
            double E68=totalE*0.68;
            double E90=totalE*0.90;
            int hitnum;
            for(hitnum=hitdata[layer_num].size()-1;hitnum>=0;hitnum--)
               if((hitdata[layer_num])[hitnum].E_suminside<E90) break;
            radius90[layer_num]+=sqrt((hitdata[layer_num])[hitnum+1].distance);
            for(;hitnum>=0;hitnum--)
               if((hitdata[layer_num])[hitnum].E_suminside<E68) break;
            radius68[layer_num]+=sqrt((hitdata[layer_num])[hitnum+1].distance);         
           }//end else
        }//end for layer_num
      }//loop on events
      
      for(unsigned layer_num=0;layer_num<layers;layer_num++){
         radius90[layer_num]/=(nEvts-nohit[layer_num]);
         radius68[layer_num]/=(nEvts-nohit[layer_num]);
      }
      
      TCanvas *myc=new TCanvas("radii","radii",800,600); 
      TGraph *r90=new TGraph(layers,layer_array,radius90);
      TGraph *r68=new TGraph(layers,layer_array,radius68);
      r68->SetLineColor(4);
      r68->SetMarkerStyle(21);
      r68->SetTitle("68 percent containment");
      r90->SetLineColor(2);
      r90->SetTitle("90 percent containment");
      myc->Divide(2,1);
      myc->cd(1);
      r90->Draw("AC*");
      myc->cd(2);
      r68->Draw("ACP");      
      myc->SaveAs(Form("plot_radii_")+et_portion+eta_portion+Form("_tanmay.png"));

    }// end for et
                    
   }// end for eta

}
  
