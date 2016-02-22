#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"

#include "../../../../userlib/include/HGCSSEvent.hh"
#include "../../../../userlib/include/HGCSSInfo.hh"
#include "../../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../../userlib/include/HGCSSSimHit.hh"
#include "../../../../userlib/include/HGCSSSamplingSection.hh"
#include "../../../../userlib/include/HGCSSGenParticle.hh"

typedef struct{
    unsigned int hitweight;//total energy of hits in the layer
    double hit_x;
    double hit_y; //hit position
    double hit_z;
}layerstat;// statistics on one layer

typedef struct{
    double init_x; //position of hit in the first layer
    double init_y;
    double init_z;
    double x_vs_layer; // slope of reconstructed line of energy center
    double y_vs_layer;
    double z_vs_layer;
    double tan_theta; //
    double tan_phi;
    double reco_eta;
}Reco_init_status;

double weight(int layer) {
  return 1.0;
}

void convert_num_to_string(unsigned int x,char numchar[]){
  if(x>=1000) {numchar[0]='\0';}
  numchar[0]='0'+x/100;
  numchar[1]='0'+(x/10)%10;
  numchar[2]='0'+x%10;
  numchar[3]='\0';
}

void createfilename(int eventnum,char filename[],int mode){
  // for a given event and given layernum, set the file name. mode=x,y,z
  char prestring[100]="layer_vs_enengy_center_et50_eta2pt1_event";
  char eventnumchar[4];
  char modestring[7]="_x_fit";
  const char filefix[5]=".png";
  switch(mode){
     case 2: modestring[1]='y';break;
     case 3 :modestring[1]='z';break;
     default: modestring[1]='x';
  }
  convert_num_to_string(eventnum,eventnumchar);
  strcat(prestring,eventnumchar);
  strcat(prestring,modestring);
  strcat(prestring,filefix);
  strcpy(filename,prestring);
}

void layer_vs_energy_center_sig_event(int eventnum,double aver_x[],double aver_y[],double aver_z[],double &rangemin,double &rangemax,double &truex,double &truey,double &truez,double &truephi,double &trueeta){
    //load the shared library for HGCSS* classes
    gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");  

    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");

    lSimTree->AddFile("/export/cmss2/paulini/CMS/HGCal/version30/HGcal_version30_model2_BOFF_et50_eta2.100.root");
    lRecTree->AddFile("/export/cmss2/paulini/CMS/HGCal/version30/Digi_version30_model2_BOFF_et50_eta2.100.root");
 
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;
    std::vector<HGCSSGenParticle> * genvec = 0; // get the generated particle vector 

    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);   

    const unsigned nEvts = lSimTree->GetEntries(); 
    //double totalE(0);
    //unsigned layerno(0);
    //unsigned min_layer(0);
    //unsigned max_layer(0);
    //min_layer=10000;
    //max_layer=0;
    
    //TH1F *num_hit_layer =new TH1F("hit number versus layer number","statistics",28,-0.5,27.5); 
    unsigned ievt(0);
    //choode a proper event(>0,<=1000) as the data to make statistics on energy/position and etc.
    //note that if we want to sum over all 1000 events, we need a 'for' cycle.
    ievt=eventnum;
    layerstat singleevt[28];
    unsigned layernum;
    for(layernum=0;layernum<=27;layernum++){
         singleevt[layernum].hitweight=0;
         singleevt[layernum].hit_x=0;
         singleevt[layernum].hit_y=0;
         singleevt[layernum].hit_z=0;
    }//initializing 
    //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    //totalE = 0;
    
    if (ievt%50==0) std::cout << " -- Processing event " << ievt << std::endl;
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);
    
    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
      const HGCSSRecoHit lHit = (*rechitvec)[iH];
      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();
      unsigned layer = lHit.layer();
      double energy = weight(layer)*lHit.energy();
      singleevt[layer].hitweight+=energy;
      singleevt[layer].hit_x+=posx*energy;
      singleevt[layer].hit_y+=posy*energy;
      singleevt[layer].hit_z+=posz*energy;
      
      //if (layer==layernum)//choose a proper layer
      //x_y_E->Fill(posx,posy,energy);
      //endif
      //num_hit_layer->Fill(layer);
    }
    //p_E_l->Fill(totalE,layerno);
    //p_E->Fill(totalE);
    //p_E->Fill(layerno);
    //}//loop on hits
    
    //x_y_E->Draw("LEGO");//draw the diagram of TH2F x_y_E // Note: if you do not need this histogram, just delete it~ 
    //num_hit_layer->Draw("LEGO");
    //char filename[100];
    //createfilename(eventnum,layernum,filename);
    //myc->SaveAs(filename);
    //cout << "min layer = " << min_layer << std::endl;
    //cout << "max layer = " << max_layer << std::endl;
    //obtain the average hit position for each layer
    for(layernum=0;layernum<=27;layernum++){
        if(singleevt[layernum].hitweight!=0){
            aver_x[layernum]=singleevt[layernum].hit_x/singleevt[layernum].hitweight;
            aver_y[layernum]=singleevt[layernum].hit_y/singleevt[layernum].hitweight;
            aver_z[layernum]=singleevt[layernum].hit_z/singleevt[layernum].hitweight;
        }
        else aver_z[layernum]=aver_x[layernum]=aver_y[layernum]=0;
    }
    for(layernum=0;layernum<27;layernum++)
          if(singleevt[layernum].hitweight>500) break;
    rangemin=(double)layernum;
    for(layernum=27;layernum>=0;layernum--)
          if(singleevt[layernum].hitweight>1000) break;
    rangemax=(double)layernum;

        if(1||!((*genvec).size() > 1)){  // Delete events which has converted photons.
           const HGCSSGenParticle gHit = (*genvec)[0];
    
	   truex = gHit.x(); //initial position 
	   truey = gHit.y();
           truez = gHit.z();
           trueeta = gHit.eta();
           truephi = gHit.phi();
        }    
}

void allevts_layer_vs_energy_center(){//main
  int eventnum=500;
  int layernum;
  double layers[28],aver_x[28],aver_y[28],aver_z[28],rangemin,rangemax;
  double init_x,x_vs_layer,init_y,y_vs_layer,init_z,z_vs_layer,tan_phi,tan_theta,reco_eta;
  double truex,truey,truez,truephi,trueeta;
  TCanvas *myc;
  TGraph *tt;
  TF1 *fit;
  char filename[100];
  FILE *fp=fopen("fit.txt","w");
  fprintf(fp,"evtnum     x         truex         y          truey         z          truez           tanphi         truetanphi     eta        trueeta\n");
  for(int k=0;k<=27;k++) 
    layers[k]=k;
  for(eventnum=0;eventnum<=999;eventnum+=100){
     layer_vs_energy_center_sig_event(eventnum,aver_x,aver_y,aver_z,rangemin,rangemax,truex,truey,truez,truephi,trueeta);
     tt =new TGraph(28,layers,aver_x);
     tt->Fit("pol1","O","",rangemin,rangemax);
     fit=tt->GetFunction("pol1");
     init_x=fit->GetParameter(0);
     x_vs_layer=fit->GetParameter(1);

     tt=new TGraph(28,layers,aver_y);
     tt->Fit("pol1","O","",rangemin,rangemax);
     fit=tt->GetFunction("pol1");
     init_y=fit->GetParameter(0);
     y_vs_layer=fit->GetParameter(1);

     tt=new TGraph(28,layers,aver_z);
     tt->Fit("pol1","O","",rangemin,rangemax);
     fit=tt->GetFunction("pol1");
     init_z=fit->GetParameter(0);
     z_vs_layer=fit->GetParameter(1);
     
     tan_phi=y_vs_layer/x_vs_layer;
     tan_theta=sqrt((y_vs_layer*y_vs_layer+x_vs_layer*x_vs_layer)/(z_vs_layer*z_vs_layer));
     reco_eta=-log(tan(atan(tan_theta)/2));
     fprintf(fp,"%d   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf\n",eventnum,init_x,truex,init_y,truey,init_z,truez,tan_phi,tan(truephi),reco_eta,trueeta);
  }
  fclose(fp);
}
