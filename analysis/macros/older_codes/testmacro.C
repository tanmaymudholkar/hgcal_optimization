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

void testmacro(){//main  

    //load the shared library for HGCSS* classes
    gSystem->Load("/home/tmudholk/research/HGCstandalone/userlib/lib/libPFCalEEuserlib.so");  

    TChain  *lSimTree = new TChain("HGCSSTree");
    TChain  *lRecTree = new TChain("RecoTree");

    lSimTree->AddFile("/export/cmss2/tmudholk/HGCal/version34/HGcal__version34_model2_BOFF_et50_eta2.100.root");
    lRecTree->AddFile("/export/cmss2/tmudholk/HGCal/version34/Digi__version34_model2_BOFF_et50_eta2.100.root");
 
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0;
    std::vector<HGCSSRecoHit> * rechitvec = 0;

    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  
    std::vector<Double_t> total_energies;
    Double_t expected_mean_energy = 0;
    Double_t expected_stdev_energy = 0;

    const unsigned nEvts = lSimTree->GetEntries(); 
    double totalE(0);

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
      total_energies.push_back(totalE);
      expected_mean_energy += totalE/nEvts;
      expected_stdev_energy += totalE*totalE/nEvts;
    }//loop on evts

    expected_stdev_energy = expected_stdev_energy - expected_mean_energy*expected_mean_energy;
    expected_stdev_energy = sqrt(expected_stdev_energy);
    
    TCanvas *myc = new TCanvas("Total E","Total E",800,600);
    myc->cd();
    TH1F *p_E = new TH1F("total E in ECAL","total_E",50,expected_mean_energy-5*expected_stdev_energy,expected_mean_energy+5*expected_stdev_energy);
    for (unsigned histogram_counter = 0; histogram_counter < total_energies.size(); histogram_counter++) {
      p_E->Fill(total_energies[histogram_counter]);
    }
    p_E->Draw();
    myc->SaveAs("energy_histogram.png");
}
