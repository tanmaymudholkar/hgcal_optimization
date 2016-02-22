#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"

#include "../../../userlib/include/HGCSSEvent.hh"
#include "../../../userlib/include/HGCSSInfo.hh"
#include "../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../userlib/include/HGCSSSimHit.hh"
#include "../../../userlib/include/HGCSSSamplingSection.hh"

double weight(int layer) {
  return 1.0;
}

void plotE_singleevt_copy(){//main  

    //load the shared library for HGCSS* classes
    gSystem->Load("/export/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");  

    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");

    lSimTree->AddFile("/export/cmss2/paulini/CMS/HGCal/version30/HGcal_version30_model2_BOFF_et50_eta2.100.root");
    lRecTree->AddFile("/export/cmss2/paulini/CMS/HGCal/version30/Digi_version30_model2_BOFF_et50_eta2.100.root");
 
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
    const unsigned nEvts = lSimTree->GetEntries(); 
    //double totalE(0);
    //unsigned layerno(0);
    //unsigned min_layer(0);
    //unsigned max_layer(0);
    //min_layer=10000;
    //max_layer=0;
    
    //TCanvas *myc = new TCanvas("Total E","Total E",800,600);
    TCanvas *myc = new TCanvas("Layers","Layers",800,800);
    myc->cd();
    //TH1F *p_E = new TH1F("Layers","layer",29,0,29);
    TH2F *x_y_E = new TH2F("","hit postion versus energy",100,670,730,100,-460,-390);
    
    //TH1F *num_hit_layer =new TH1F("hit number versus layer number","statistics",28,-0.5,27.5); 

    unsigned ievt(0);
    //choode a proper event(>0,<=1000) as the data to make statistics on energy/position and etc.
    //note that if we want to sum over all 1000 events, we need a 'for' cycle.
    ievt=500;

    //for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
    //totalE = 0;
    
    if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
    lSimTree->GetEntry(ievt);
    lRecTree->GetEntry(ievt);

    for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
      const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
      double posx = lHit.get_x();
      double posy = lHit.get_y();
      double posz = lHit.get_z();
      unsigned layer = lHit.layer();
      double energy = weight(layer)*lHit.energy();
	  
      //totalE += energy;
      //layerno += layer;
      //if(layer < min_layer) min_layer = layer;
      //if(layer > max_layer) max_layer = layer;
      //cout << "layer is " << layer << std::endl;
      if (layer==27)//choose a proper layer
      x_y_E->Fill(posx,posy,energy);
      //endif
      //num_hit_layer->Fill(layer);
    }
    //p_E_l->Fill(totalE,layerno);
    //p_E->Fill(totalE);
    //p_E->Fill(layerno);
    //}//loop on hits
    
    x_y_E->Draw("LEGO");//draw the diagram of TH2F x_y_E // Note: if you do not need information from this histogram, just delete it~ 
    //num_hit_layer->Draw("LEGO");

    myc->SaveAs("position_et50_eta2pt1_layer4.png");
    //cout << "min layer = " << min_layer << std::endl;
    //cout << "max layer = " << max_layer << std::endl;

}
