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
#include "../../../userlib/include/HGCSSGenParticle.hh"


void plotXY(){//main  

    //load the shared library for HGCSS* classes
    //change the path to your own directory
    gSystem->Load("../../../userlib/lib/libPFCalEEuserlib.so");  

    TChain  *lSimTree = new TChain("HGCSSTree");

    lSimTree->AddFile("../export_dir/version34/HGcal__version34_model2_BOFF_et10_eta2.100.root");
 
    std::vector<HGCSSSamplingSection> * ssvec = 0;
    std::vector<HGCSSSimHit> * simhitvec = 0; // get the simulated hit vector
    std::vector<HGCSSGenParticle> * genvec = 0; // get the generated particle vector 

    lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec); 

    const unsigned nEvts = lSimTree->GetEntries(); 
  
    TCanvas *genc = new TCanvas("GenParticle position","GenParticle position",1000,800);
    genc->Divide(2,2);
    TCanvas *simc1 = new TCanvas("Sim position","Sim position",1000,800);
    simc1->Divide(6,5);
    TCanvas *simc2 = new TCanvas("unweighted Sim position","unweighted Sim position",1000,800);
    simc2->Divide(6,5);

    TH2F *p_gen_xy = new TH2F("p_gen_xy","Gen particle position",340,-1700,1700,340,-1700,1700);
    TH1F *p_gen_eta = new TH1F("p_gen_eta","Gen particle eta",160,1.4,3.0);
    TH1F *p_gen_phi = new TH1F("p_gen_phi","Gen particle phi",200,0,2*3.14);

    std::ostringstream simlayer;
    TH2F *p_sim[28];  // one 2D histogram for each layer. p_sim is for weighted hit position, and p_unweight is for unweight hit, i.e. All hits.
    TH2F *p_unweight[28];
    for(unsigned iL(0); iL <28; iL ++){
      simlayer.str("");
      simlayer << "layer_" << iL ;
      p_sim[iL] =  new TH2F(simlayer.str().c_str(), simlayer.str().c_str(),340,-1700,1700,340,-1700,1700);
      simlayer.str("");
      simlayer << "unweight_layer_" << iL ;
      p_unweight[iL] = new TH2F(simlayer.str().c_str(), simlayer.str().c_str(),340,-1700,1700,340,-1700,1700);
   }
    
    for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

      if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
      lSimTree->GetEntry(ievt);

      if((*genvec).size() > 1) continue;  // Delete events which has converted photons.
      else{   
        const HGCSSGenParticle gHit = (*genvec)[0];
    
	double posx = gHit.x();
	double posy = gHit.y();
        double eta = gHit.eta();
        double phi = gHit.phi();

        p_gen_xy->Fill(posx, posy);
        p_gen_eta->Fill(eta);
        p_gen_phi->Fill(phi);
      }

      double WeightCenter[28][2];
      double layerE[28];
      for(unsigned i(0); i<28; i++){
         WeightCenter[i][0] = 0;
         WeightCenter[i][1] = 0;
         layerE[i] = 0;
      }
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop over rechits
        const HGCSSSimHit lHit = (*simhitvec)[iH];
    
	double posx = lHit.get_x();
	double posy = lHit.get_y();
	unsigned layer = lHit.layer();
	double energy = lHit.energy();
   
        p_unweight[layer]->Fill(posx, posy);  // unweighted position
      
        layerE[layer] += energy;
	WeightCenter[layer][0] += energy*posx;  //use energy as weight for the position
        WeightCenter[layer][1] += energy*posy;
      }
      
      for(unsigned iL(0); iL < 28; iL++){
         if(layerE[iL]!=0) p_sim[iL]->Fill(WeightCenter[iL][0]/layerE[iL], WeightCenter[iL][1]/layerE[iL]);   
         else p_sim[iL]->Fill(WeightCenter[iL][0], WeightCenter[iL][1]);   
      }
    }//loop on hits
    
    for(unsigned iL(0); iL < 28; iL++){
       simc1->cd(iL+1);
       p_sim[iL]->Draw();
       simc2->cd(iL+1);
       p_unweight[iL]->Draw("colz");
    }
    genc->cd(1);
    p_gen_xy->Draw();
    genc->cd(2);
    p_gen_eta->Draw();
    genc->cd(3);
    p_gen_phi->Draw();
    genc->SaveAs("GenPosition.png");
    simc1->SaveAs("SimPosition.png");
    simc2->SaveAs("unweighted_SimPosition.png");
}



