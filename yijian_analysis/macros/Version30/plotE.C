#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TChain.h"
#include "TSystem.h"

#include "../../userlib/include/HGCSSEvent.hh"
#include "../../userlib/include/HGCSSInfo.hh"
#include "../../userlib/include/HGCSSRecoHit.hh"
#include "../../userlib/include/HGCSSSimHit.hh"
#include "../../userlib/include/HGCSSSamplingSection.hh"

void plotE(){//main  

    //load the shared library for HGCSS* classes
    gSystem->Load("/afs/cern.ch/work/m/msun/TEST/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");  

    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");

    lSimTree->AddFile("HGcal_version12_model2_BOFF_et20_alpha0.110.root");
    lRecTree->AddFile("Digi_version12_model2_BOFF_et20_alpha0.110.root");
 
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
    const unsigned nEvts = lSimTree->GetEntries(); 
    double totalE(0);
  
    TCanvas *myc = new TCanvas("Total E","Total E",800,600);
    myc->cd();
    TH1F *p_E = new TH1F("total E in ECAL","total_E",1000,15000,25000);

    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries
      totalE = 0;

      if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
      lSimTree->GetEntry(ievt);
      lRecTree->GetEntry(ievt);

      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){//loop over rechits
        const HGCSSRecoHit lHit = (*rechitvec)[iH];
    
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	double posz = lHit.get_z();
	unsigned layer = lHit.layer();
	double energy = lHit.energy();
	  
	totalE += energy;
      }
      p_E->Fill(totalE);      
    }//loop on hits
    
    p_E->Draw();
    myc->SaveAs("totalE.png");

}



