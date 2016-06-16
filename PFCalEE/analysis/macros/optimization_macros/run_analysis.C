// Refactored version of the code, supposed to be all-encompassing. Assumes "hexagonal" geometry, and that the digi files are generated with a threshold of 0.5 mips.

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
#include "TLine.h"

#include "../../../userlib/include/HGCSSEvent.hh"
#include "../../../userlib/include/HGCSSInfo.hh"
#include "../../../userlib/include/HGCSSRecoHit.hh"
#include "../../../userlib/include/HGCSSSimHit.hh"
#include "../../../userlib/include/HGCSSSamplingSection.hh"
#include "../../../userlib/include/HGCSSGenParticle.hh"

const Double_t sigmas_down = 3.0;
const Double_t sigmas_up = 3.0;
const Int_t maximum_nevts_hits_distribution = 1000;
// const Double_t abovecut = 21000;
// const Double_t belowcut = 20500;
// const Int_t signal_region = 1000;

Double_t chisqr;
Double_t ndfr;
Double_t chisqdf;  

// Double_t et_values_array[] = {5,50,100,150}; //For systematic optimization studies
// Double_t et_values_array[] = {3,5,7,10,20,30,40,50,60,70,100,125,150}; //Intersection version 30 & 34
// Double_t et_values_array[] = {3,5,7,10,20};
// Double_t et_values_array[] = {3,5,10,20,30,40,50,60,80,100,125,175};
// Double_t et_values_array[] = {3,5,7,10,20,30,40,50,60,70,80,90,100,125,150,175,200}; // version 34 only
// Double_t et_values_array[] = {3,7,20,35,50,70,100,125};
// Double_t et_values_array[] = {100,125,150,175,200};
// Double_t et_values_array[] = {3};
// Double_t et_values_array[] = {3,5,20,50,100,150}; // high stats vflat and version 34 intersection
// Double_t et_values_array[] = {3,5,20,50,100}; // high stats version 34 only
Double_t et_values_array[] = {3,5,10,30,50,70,100,150}; // version 33 and 30, hexagonal geometry
// Double_t et_values_array[] = {3,5,10,100}; // version 33 and 30, hexagonal geometry
// Double_t et_values_array[] = {3,5,7,10,20,30,40,50,60,70,100,125,150}; // old version 30
// Double_t et_values_array[] = {3,7,10,40,100}; // old version 30 reduced
// Double_t et_values_array[] = {3,5,10,100}; // version 33 and 30, hexagonal geometry, reduced
Double_t eta_values_array[] = {2.1};
  
std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
std::vector<Double_t> eta_values(eta_values_array,eta_values_array+sizeof(eta_values_array)/sizeof(Double_t));

std::vector<Double_t> weights;

Int_t version_number;
TString version_name;
TString datadir;
TString outputdir;
Double_t threshold_adc;
Int_t digi_or_raw_switch;
Bool_t use_modified_weighting_scheme;
Bool_t multiple_runs_present;
Bool_t get_hits_distribution;
Bool_t get_shower_profile;
Bool_t get_xy_positions_gen_particle;
Bool_t get_xy_positions_by_layer;
unsigned max_events;

Double_t threshold_mips;
TString eta_portion;
TString et_portion;

// Double_t resolutions_fit(Double_t *x,Double_t *par) {
//   return (par[1] + par[0]/sqrt(x[0]));
// }

Double_t resolutions_fit_standalone_withnoise(Double_t *x,Double_t *par) {
  return sqrt(pow(par[1],2) + pow(par[0],2)/x[0] + pow(par[2]/x[0],2));
  // return sqrt(pow(par[1],2) + pow(par[0],2)/x[0]);
}

Double_t resolutions_fit_standalone(Double_t *x,Double_t *par) {
  // return sqrt(pow(par[1],2) + pow(par[0],2)/x[0] + pow(par[2]/x[0],2));
  return sqrt(pow(par[1],2) + pow(par[0],2)/x[0]);
}

// Double_t resolutions_fit_standalone_const(Double_t *x,Double_t *par) {
//   return sqrt(pow(0.009,2) + pow(par[0],2)/x[0] + pow(par[1]/x[0],2));
//   // return sqrt(pow(0.009,2) + pow(par[0],2)/x[0]);
// }

// Double_t cellsizeat(Double_t x, Double_t y) {
//   Double_t r = sqrt(pow(x,2)+pow(y,2));
//   if (r<radlim) {
//     return (3.0*cellSize);
//   }
//   else {
//     return (4.0*cellSize);
//   }
// }

bool testInputFile(TString inputPath, TFile* testFile) {
  testFile = TFile::Open(inputPath);
  if ( !testFile ) {
    // std::cout << " -- Error, input file " << inputPath << " cannot be opened. Skipping..." << std::endl;
    return false;
  } // else std::cout << " -- input file " << testFile->GetName() << " successfully opened." << std::endl;
  delete testFile;
  return true;
}

void read_weights(Int_t version_number, TString datadir, Bool_t use_modified_weighting_scheme) {
  std::vector<Double_t> weights_raw;
  ifstream f_layer_weights_raw;
  f_layer_weights_raw.open(datadir+Form("/layer_weights_version%i.dat",version_number));
  std::string line;
  Double_t weight;
  if (f_layer_weights_raw.is_open()) {
    while (getline(f_layer_weights_raw,line)) {
      weight = std::atof(line.c_str());
      weights_raw.push_back(weight);
    }
  }
  f_layer_weights_raw.close();
  if (weights_raw.empty()) {
    std::cout << "Unable to read weights from file" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // std::cout << std::endl;
  // std::cout << std::endl;
  // std::cout << std::endl;
  // std::cout << std::endl;
  // std::cout << "Old weights:" << std::endl;

  // for (unsigned layer_counter = 0; layer_counter != weights_raw.size(); layer_counter++) {
  //   std::cout << "layer " << layer_counter << ": " << weights_raw[layer_counter] << std::endl;
  // }

  if (use_modified_weighting_scheme) {
    weights.push_back(weights_raw[0]+0.5*weights_raw[1]);
    for (unsigned layer_counter = 1; layer_counter != (weights_raw.size()-1); layer_counter++) {
      weights.push_back(0.5*(weights_raw[layer_counter]+weights_raw[layer_counter+1]));
    }
    weights.push_back(0.5*weights_raw[(weights_raw.size()-1)]);
  }
  else {
    for (unsigned layer_counter = 0; layer_counter != (weights_raw.size()); layer_counter++) {
      weights.push_back(weights_raw[layer_counter]);
    }
  }
  
  std::cout << std::endl;
  std::cout << "Weights:" << std::endl;

  for (unsigned layer_counter = 0; layer_counter != weights.size(); layer_counter++) {
    std::cout << "layer " << layer_counter << ": " << weights[layer_counter] << std::endl;
  }
}

void read_input_files(TChain *lSimTree, TChain *lRecTree) {

  TString Digi_common_prefix, HGcal_common_prefix;

  TString common_suffix = Form(".root");
      
  if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) {
    // WHEN YOU CHANGE THIS REMEMBER ALSO TO CHANGE CONDITION REGARDING THRESHOLD IN TOTAL ENERGY CALCULATION
    Digi_common_prefix = datadir + Form("/DigiIC3_thr5.0__version%i_model2_BOFF_",version_number);
    // Digi_common_prefix = datadir + Form("/DigiIC3__version%i_model2_BOFF_",version_number);
  }
  if (digi_or_raw_switch == 2 || digi_or_raw_switch == 3) {
    HGcal_common_prefix = datadir + Form("/HGcal__version%i_model2_BOFF_",version_number);
  }

  if (multiple_runs_present) {
    unsigned run_no = 1;
    Int_t consecutive_nonexistent_files = 0;

    while (consecutive_nonexistent_files < 4) {
      std::cout << "Reading input files for run number " << run_no << std::endl;
      
      TFile *testFile_hgcal(0);
      TFile *testFile_digi(0);

      // if (digi_or_raw_switch == 2) {
      //   TFile *testFile_hgcal(0);
      // }
      // if (digi_or_raw_switch == 1) {
      //   TFile *testFile_digi(0);
      // }

      bool hgcal_data_exists, digi_data_exists;
      if (digi_or_raw_switch == 2 || digi_or_raw_switch == 3) {
        hgcal_data_exists=testInputFile(HGcal_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix, testFile_hgcal);
      }
      if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) {
        digi_data_exists=testInputFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix, testFile_digi);
      }

      if (digi_or_raw_switch == 1) {
        if (digi_data_exists) {
          // if (hgcal_data_exists) {
          consecutive_nonexistent_files = 0;
          // std::cout << "Exists!" << std::endl;
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
          // lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
          lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);

          // TFile *inputFile = TFile::Open(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
          // HGCSSInfo * info=(HGCSSInfo*)inputFile->Get("Info");
          // cellSize = info->cellSize();
          // const unsigned versionNumber = info->version();
          // const unsigned model = info->model();
          // delete inputFile;
        }
        else {
          consecutive_nonexistent_files++;
          std::cout << "Digi data does not exist for run no " << run_no << std::endl;
          //     // std::cout << "HGcal data does not exist" << std::endl; 
        }
      }

      if (digi_or_raw_switch == 2) {
        if (hgcal_data_exists) {
          consecutive_nonexistent_files = 0;
          // std::cout << "Exists!" << std::endl;
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
          // lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
          // lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);

          // TFile *inputFile = TFile::Open(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
          // HGCSSInfo * info=(HGCSSInfo*)inputFile->Get("Info");
          // cellSize = info->cellSize();
          // const unsigned versionNumber = info->version();
          // const unsigned model = info->model();
          // delete inputFile;
        }
        else {
          consecutive_nonexistent_files++;
          std::cout << "HGcal data does not exist for run no " << run_no << std::endl;
          //     // std::cout << "HGcal data does not exist" << std::endl; 
        }
      }

      if (digi_or_raw_switch == 3) {
        if (digi_data_exists && hgcal_data_exists) {
          consecutive_nonexistent_files = 0;
          // std::cout << "Exists!" << std::endl;
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
          // lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
          // lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+Form("_run%i",run_no)+common_suffix);
          // lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);

          // TFile *inputFile = TFile::Open(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
          // HGCSSInfo * info=(HGCSSInfo*)inputFile->Get("Info");
          // cellSize = info->cellSize();
          // const unsigned versionNumber = info->version();
          // const unsigned model = info->model();
          // delete inputFile;
        }
        else {
          consecutive_nonexistent_files++;
          if (!(hgcal_data_exists)) {
            std::cout << "HGcal data does not exist for run no " << run_no << std::endl;
          }
          if (!(digi_data_exists)) {
            std::cout << "Digi data does not exist for run no " << run_no << std::endl;
          }
        }
      }

      // if (digi_or_raw_switch == 1) {
      //   delete testFile_digi;
      // }
      // if (digi_or_raw_switch == 2) {
      //   delete testFile_hgcal;
      // }
      delete testFile_digi;
      delete testFile_hgcal;
      run_no++;
    }// ends loop over runs
  } // ends if condition checking whether multiple runs are present
  else {
    std::cout << "Reading input files..." << std::endl;
      
    TFile *testFile_hgcal(0);
    TFile *testFile_digi(0);

    bool hgcal_data_exists, digi_data_exists;
    if (digi_or_raw_switch == 2 || digi_or_raw_switch == 3) {
      hgcal_data_exists=testInputFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix, testFile_hgcal);
    }
    if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) {
      digi_data_exists=testInputFile(Digi_common_prefix+et_portion+eta_portion+common_suffix, testFile_digi);
    }

    if (digi_or_raw_switch == 1) {
      if (digi_data_exists) {
        lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
      }
      else {
        std::cout << "Digi data does not exist" << std::endl;
      }
    }

    if (digi_or_raw_switch == 2) {
      if (hgcal_data_exists) {
        lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
      }
      else {
        std::cout << "HGcal data does not exist" << std::endl; 
      }
    }

    if (digi_or_raw_switch == 3) {
      if (digi_data_exists && hgcal_data_exists) {
        lSimTree->AddFile(HGcal_common_prefix+et_portion+eta_portion+common_suffix);
        lRecTree->AddFile(Digi_common_prefix+et_portion+eta_portion+common_suffix);
      }
      else {
        if (!(hgcal_data_exists)) {
          std::cout << "HGcal data does not exist" << std::endl;
        }
        if (!(digi_data_exists)) {
          std::cout << "Digi data does not exist" << std::endl;
        }
      }
    }
    delete testFile_digi;
    delete testFile_hgcal;
  } // end loop entered into if multiple runs are not present
}

void set_et_portion(Double_t et) {
  et_portion = Form("et%.0f",et);
}

void analyze_for_given_et(Double_t et,
                          Double_t eta,
                          TFile *histograms_output_file,
                          std::vector<Double_t> &energies_incoming_gev,
                          std::vector<Double_t> &energies_incoming_gev_errors,
                          std::vector<Double_t> &mean_energies_wmips,
                          std::vector<Double_t> &mean_energies_wmips_errors,
                          std::vector<Double_t> &sigmas_wmips,
                          std::vector<Double_t> &sigmas_wmips_errors,
                          std::vector<Int_t> &events_measured
                          ) {

  set_et_portion(et);

  //TString Digi_common_prefix = Digi_common_prefix_firstpart + Form("rwcuf_%.1f_rwcum_%.1f_",rwcuf,rwcum);

  Double_t energy_incoming_gev = et*cosh(eta);
  Double_t energy_incoming_gev_error = 0;

  Double_t mean_energy_wmips;
  Double_t mean_energy_wmips_error;
  Double_t sigma_wmips;
  Double_t sigma_wmips_error;
  std::vector<Double_t> total_energies;
  Double_t total_energies_by_layer[weights.size()];
  std::vector<std::pair<Double_t, Double_t>> xy_positions_by_layer[weights.size()];
  std::vector<std::pair<Double_t, Double_t>> xy_positions_gen_particle;
  Double_t meanE_statistical = 0;
  Double_t sigE_statistical = 0;
  unsigned nEvts_to_count = 0;
  // std::cout << "et" << et << " eta" << eta << std::endl;

  TChain  *lSimTree = new TChain("HGCSSTree");
  TChain  *lRecTree = new TChain("RecoTree");

  read_input_files (lSimTree, lRecTree);
      
  // std::cout << "cellSize is " << cellSize << std::endl;
  // std::cout << "versionNumber is " << versionNumber << std::endl;
  // std::cout << "model is " << model << std::endl;

  std::vector<HGCSSRecoHit> *rechitvec = 0;
  std::vector<HGCSSSimHit> *simhitvec = 0;
  // std::vector<HGCSSSamplingSection> * ssvec = 0;
  std::vector<HGCSSGenParticle> *genvec = 0; // get the generated particle vector
  
  if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) {
    lRecTree->SetBranchAddress("HGCSSRecoHitVec",&rechitvec);
  }
  if (digi_or_raw_switch == 2 || digi_or_raw_switch == 3) {
    lSimTree->SetBranchAddress("HGCSSSimHitVec",&simhitvec);
    lSimTree->SetBranchAddress("HGCSSGenParticleVec",&genvec);
  }
  // lSimTree->SetBranchAddress("HGCSSSamplingSectionVec",&ssvec);
      
  unsigned nEvts_total;
  if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) {
    nEvts_total = lRecTree->GetEntries();
  }
  if (digi_or_raw_switch == 2) {
    nEvts_total = lSimTree->GetEntries();
  }

  unsigned nEvts_over_which_to_loop = (max_events == 0? nEvts_total : (nEvts_total < max_events? nEvts_total : max_events));
  
  Double_t totalE(0);

  Double_t totalE_by_layer[weights.size()];
      
  // std::cout << "no of events over which to loop = " << nEvts_over_which_to_loop << std::endl;

  TCanvas *c_individual_hits_energy_distribution = new TCanvas(Form("c_individual_hits_energy_distribution_")+et_portion,Form("individual_hits_energy_distribution_"),800,600);
  c_individual_hits_energy_distribution->cd();
  TH1F *h_individual_hits_energy_distribution = new TH1F(Form("h_individual_hits_energy_distribution_")+et_portion, Form("individual_hits_energy_distribution_")+et_portion+":hit_raw_energy", 200, 0, 0);
  // if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) p_l_temp = new TH1F("Hits Energy Distribution", "Energies", 200, 0, 30);
  // if (digi_or_raw_switch == 2) p_l_temp = new TH1F("Hits Energy Distribution", "Energies", 200, 0, 5);
  // if (get_hits_distribution) h_individual_hits_energy_distribution

  ofstream o_total_energies;
  TString o_total_energies_name = outputdir+Form("/data/data_total_energies_")+version_name+et_portion+eta_portion+Form("_thr%.1f",threshold_mips);
  o_total_energies.open(o_total_energies_name);
  for (unsigned ievt(0); ievt<nEvts_over_which_to_loop; ++ievt){// loop on events
    // for (unsigned ievt(0); ievt<1000; ++ievt){// loop on events
    totalE = 0;
    // Double_t totalE_with_5_mips_threshold = 0;
    std::vector<std::pair<Double_t, Double_t>> xy_positions_by_layer_for_this_event[weights.size()];

    if (get_shower_profile) {
      for (unsigned int layer_counter = 0; layer_counter != weights.size(); layer_counter++) {
        totalE_by_layer[layer_counter] = 0;
      }
    }
	
    if (ievt%100==0) std::cout << " -- Processing event " << ievt << std::endl;
	
    if (digi_or_raw_switch == 1 || digi_or_raw_switch == 3) {
      lRecTree->GetEntry(ievt);
    }
    if (digi_or_raw_switch == 2 || digi_or_raw_switch == 3) {
      lSimTree->GetEntry(ievt);
    }

    //std::cout << "track id " << (*genvec)[0].trackID() << "  genvec size " << (*genvec).size() << std::endl;

    // if((*genvec).size() > 1) continue;  // Delete events which has converted photons.
    // if((*genvec)[0].trackID() >= 2) continue;
    // else{
    // nEvts_to_count += 1;
    // const HGCSSGenParticle gHit = (*genvec)[0];
    // std::cout << "ACCEPTED: " << "track id " << (*genvec)[0].trackID() << "  genvec size " << (*genvec).size() << std::endl;
    // Double_t posx_gen_ini = ((*genvec)[0]).x();
    // Double_t posy_gen_ini = ((*genvec)[0]).y();
    // Double_t posz_gen_ini = ((*genvec)[0]).z();
    // Double_t px_gen = ((*genvec)[0]).px();
    // Double_t py_gen = ((*genvec)[0]).py();
    // Double_t pz_gen = ((*genvec)[0]).pz();

    if (digi_or_raw_switch == 1) {
      nEvts_to_count += 1;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
        Double_t posx = 0;
        Double_t posy = 0;
        if (get_xy_positions_by_layer) {
          posx = ((*rechitvec)[iH]).get_x();
          posy = ((*rechitvec)[iH]).get_y();
        }
        // Double_t posz = ((*rechitvec)[iH]).get_z();
        unsigned layer = ((*rechitvec)[iH]).layer();
        Double_t energy = ((*rechitvec)[iH]).energy();
            
        if (get_hits_distribution && ievt<=maximum_nevts_hits_distribution && energy>0) {
          h_individual_hits_energy_distribution->Fill(energy);
        }
        
        if(energy>threshold_mips) {
          std::pair<Double_t, Double_t> xypair (posx, posy);
          if (get_xy_positions_by_layer) xy_positions_by_layer_for_this_event[layer].push_back(xypair);
          // std::cout << "layer number: " << layer << std::endl;
          // Double_t posx_gen = posx_gen_ini + (px_gen/pz_gen)*(posz - posz_gen_ini);
          //   Double_t posy_gen = posy_gen_ini + (py_gen/pz_gen)*(posz - posz_gen_ini);
          //   Double_t dx = posx - posx_gen;
          //   Double_t dy = posy - posy_gen;
          //   Double_t halfCell = 0.5*cellsizeat(posx,posy);
          //   if ((fabs(dx) <= signal_region*halfCell) && (fabs(dy) <= signal_region*halfCell)){
          Double_t weighted_energy = energy*weights[layer]/tanh(eta);
          // totalE_in_layer[layer] += weighted_energy;
          totalE += weighted_energy;
          if (get_shower_profile) totalE_by_layer[layer] += weighted_energy;
          //  } // end if condition for counting energies in a signal region
        } //end if condition for counting energies above a threshold
      } // end loop over rechits
    } // End loop entered into if only digi hits to be analyzed
    if (digi_or_raw_switch == 3) {
      if((*genvec)[0].trackID() >= 2) continue;
      Double_t posx_gen_ini = 0;
      Double_t posy_gen_ini = 0;
      if (get_xy_positions_gen_particle) {
        posx_gen_ini = ((*genvec)[0]).x();
        posy_gen_ini = ((*genvec)[0]).y();
      }
      if (posx_gen_ini > 200) continue;
      std::pair<Double_t, Double_t> xy_gen_pair (posx_gen_ini, posy_gen_ini);
      if (get_xy_positions_gen_particle) xy_positions_gen_particle.push_back(xy_gen_pair);
      nEvts_to_count += 1;
      for (unsigned iH(0); iH<(*rechitvec).size(); ++iH){ // loop over rechits
        Double_t posx = 0;
        Double_t posy = 0;
        if (get_xy_positions_by_layer) {
          posx = ((*rechitvec)[iH]).get_x();
          posy = ((*rechitvec)[iH]).get_y();
        }
        // Double_t posz = ((*rechitvec)[iH]).get_z();
        unsigned layer = ((*rechitvec)[iH]).layer();
        Double_t energy = ((*rechitvec)[iH]).energy();
            
        if (get_hits_distribution && ievt<=maximum_nevts_hits_distribution && energy>0) {
          h_individual_hits_energy_distribution->Fill(energy);
        }
        
        if(energy>threshold_mips) {
          std::pair<Double_t, Double_t> xypair (posx, posy);
          if (get_xy_positions_by_layer) xy_positions_by_layer_for_this_event[layer].push_back(xypair);
          // std::cout << "layer number: " << layer << std::endl;
          // Double_t posx_gen = posx_gen_ini + (px_gen/pz_gen)*(posz - posz_gen_ini);
          //   Double_t posy_gen = posy_gen_ini + (py_gen/pz_gen)*(posz - posz_gen_ini);
          //   Double_t dx = posx - posx_gen;
          //   Double_t dy = posy - posy_gen;
          //   Double_t halfCell = 0.5*cellsizeat(posx,posy);
          //   if ((fabs(dx) <= signal_region*halfCell) && (fabs(dy) <= signal_region*halfCell)){
          Double_t weighted_energy = energy*weights[layer]/tanh(eta);
          // totalE_in_layer[layer] += weighted_energy;
          totalE += weighted_energy;
          if (get_shower_profile) totalE_by_layer[layer] += weighted_energy;
          //  }// end if condition for counting energies in a signal region
        }//end if condition for counting energies above a threshold
      }// end loop over rechits
    } // End loop entered into if both digi and raw hits are to be analyzed
    if (digi_or_raw_switch == 2) {
      // if((*genvec).size() > 1) continue;  // Delete events which has converted photons.
      if((*genvec)[0].trackID() >= 2) continue;
      Double_t posx_gen_ini = 0;
      Double_t posy_gen_ini = 0;
      if (get_xy_positions_gen_particle) {
        Double_t posx_gen_ini = ((*genvec)[0]).x();
        Double_t posy_gen_ini = ((*genvec)[0]).y();
      }
      if (posx_gen_ini > 200) continue;
      std::pair<Double_t, Double_t> xy_gen_pair (posx_gen_ini, posy_gen_ini);
      if (get_xy_positions_gen_particle) xy_positions_gen_particle.push_back(xy_gen_pair);
      nEvts_to_count += 1;
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){// loop over simhits
        // Double_t posx = ((*rechitvec)[iH]).get_x();
        // Double_t posy = ((*rechitvec)[iH]).get_y();
        // Double_t posz = ((*rechitvec)[iH]).get_z();
        unsigned layer = ((*simhitvec)[iH]).layer();
        Double_t energy = ((*simhitvec)[iH]).energy();
        // unsigned layer = ((*simhitvec)[iH]).layer();
        // Double_t energy = ((*simhitvec)[iH]).energy();

        if (get_hits_distribution && ievt<=maximum_nevts_hits_distribution && energy>0) {
          h_individual_hits_energy_distribution->Fill(energy);
        }
            
        if(energy>threshold_mips) {
          // std::pair<Double_t, Double_t> xypair (posx, posy);
          // xy_positions_by_layer[layer].push_back(xypair);
          // std::cout << "layer number: " << layer << std::endl;
          // Double_t posx_gen = posx_gen_ini + (px_gen/pz_gen)*(posz - posz_gen_ini);
          //   Double_t posy_gen = posy_gen_ini + (py_gen/pz_gen)*(posz - posz_gen_ini);
          //   Double_t dx = posx - posx_gen;
          //   Double_t dy = posy - posy_gen;
          //   Double_t halfCell = 0.5*cellsizeat(posx,posy);
          //   if ((fabs(dx) <= signal_region*halfCell) && (fabs(dy) <= signal_region*halfCell)){
          Double_t weighted_energy = energy*weights[layer]/tanh(eta);
          // totalE_in_layer[layer] += weighted_energy;
          totalE += weighted_energy;
          if (get_shower_profile) totalE_by_layer[layer] += weighted_energy;
          //  }// end if condition for counting energies in a signal region
        }//end if condition for counting energies above a threshold
      }// end loop over simhits
    }
	
    total_energies.push_back(totalE);
    if (get_shower_profile || get_xy_positions_by_layer) {
      for (unsigned int layer_counter = 0; layer_counter != weights.size(); layer_counter++) {
        if (get_shower_profile) total_energies_by_layer[layer_counter] += totalE_by_layer[layer_counter];
        if (get_xy_positions_by_layer) {
          for (std::vector<std::pair<Double_t, Double_t>>::iterator pair_iterator = xy_positions_by_layer_for_this_event[layer_counter].begin(); pair_iterator != xy_positions_by_layer_for_this_event[layer_counter].end(); ++pair_iterator) {
            xy_positions_by_layer[layer_counter].push_back(*pair_iterator);
          }
        }
      }
    }
    
    meanE_statistical += totalE;
    sigE_statistical += totalE*totalE;
    
    o_total_energies << ievt << "    " << totalE << std::endl;
  }// loop on events
  o_total_energies.close();
  // if (digi_or_raw_switch == 2) {
  //   delete lSimTree;
  // }
  // if (digi_or_raw_switch == 1) {
  //   delete lRecTree;
  // }

  delete lSimTree;
  delete lRecTree;

  if (get_hits_distribution) {
    gPad->SetLogy();
    h_individual_hits_energy_distribution->SetTitle(version_name);
    h_individual_hits_energy_distribution->Draw();
  }
        
  TLine *lthr5 = new TLine(0.5,1,0.5,100000);
  TLine *lthr20 = new TLine(2,1,2,100000);
  TLine *lthr50 = new TLine(5,1,5,100000);
  if (get_hits_distribution) {
    lthr5->Draw();
    lthr20->Draw();
    lthr50->Draw();
    histograms_output_file->WriteTObject(c_individual_hits_energy_distribution);
  }
  // myc_temp->Print(outputdir+Form("/plots/plot_")+version_name+Form("_energy_distribution_")+et_portion+eta_portion+Form("_thr%.1f",threshold_mips)+Form(".png"));
  // myc_temp->Print(outputdir+Form("/plots/plot_")+version_name+Form("_energy_distribution_")+et_portion+eta_portion+Form("_thr%.1f",threshold_mips)+Form(".pdf"));
        
  delete lthr50;
  delete lthr20;
  delete lthr5;
      
  delete h_individual_hits_energy_distribution;
  delete c_individual_hits_energy_distribution;

  meanE_statistical = meanE_statistical/nEvts_to_count;
  sigE_statistical = sigE_statistical/nEvts_to_count;
      
  sigE_statistical = sigE_statistical - meanE_statistical*meanE_statistical;
  sigE_statistical = sqrt(sigE_statistical);

  // std::cout << "mean_statistical is " << meanE_statistical << std::endl;
  // std::cout << "sigma_statistical is " << sigE_statistical << std::endl;

  Double_t lower_bound_for_fit;
  Double_t upper_bound_for_fit;
      
  lower_bound_for_fit = meanE_statistical - sigmas_down*sigE_statistical;
  upper_bound_for_fit = meanE_statistical + sigmas_up*sigE_statistical;

  TCanvas *c_total_energy_distribution = new TCanvas(Form("c_total_energy_distribution")+et_portion,Form("total_energy_distribution")+et_portion,800,600);
  c_total_energy_distribution->cd();
  TH1F *h_total_energy_distribution = new TH1F(Form("h_total_energy_distribution")+et_portion, Form("total_energy_distribution")+et_portion+":total deposited energy(weighted)", 50, lower_bound_for_fit,upper_bound_for_fit);
  for(unsigned hist_filler_counter = 0; hist_filler_counter < total_energies.size(); hist_filler_counter++) {
    h_total_energy_distribution->Fill(total_energies[hist_filler_counter]);
  }
  h_total_energy_distribution->Fit("gaus", "IME", "", lower_bound_for_fit, upper_bound_for_fit);
  TF1 *gaussian_fit = h_total_energy_distribution->GetFunction("gaus");
  mean_energy_wmips = gaussian_fit->GetParameter(1);
  mean_energy_wmips_error = gaussian_fit->GetParError(1);
  sigma_wmips = gaussian_fit->GetParameter(2);
  sigma_wmips_error = gaussian_fit->GetParError(2);
  chisqr = gaussian_fit->GetChisquare();
  ndfr = gaussian_fit->GetNDF();
  std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
  chisqdf = chisqr/ndfr;
  // resolution = stddev/mean_energy;
  // resolution_error = resolution*((stddev_error/stddev) + (mean_energy_error/mean_energy));
  // resolution_error = resolution*sqrt(pow(stddev_error/stddev,2) + pow(mean_energy_error/mean_energy,2));
  // energy_incoming_gev_error = energy_incoming_gev*mean_energy_error/mean_energy;
  // energy_incoming_gev_error = 0;
  energies_incoming_gev.push_back(energy_incoming_gev);
  energies_incoming_gev_errors.push_back(energy_incoming_gev_error);
  mean_energies_wmips.push_back(mean_energy_wmips);
  mean_energies_wmips_errors.push_back(mean_energy_wmips_error);
  sigmas_wmips.push_back(sigma_wmips);
  sigmas_wmips_errors.push_back(sigma_wmips_error);
  events_measured.push_back(nEvts_to_count);
  // resolutions.push_back(resolution);
  // resolutions_errors.push_back(resolution_error);
      
  h_total_energy_distribution->SetTitle(version_name + Form(" . chisq/dof = %.2f",chisqdf));
  h_total_energy_distribution->Draw();
  // p_l_unweighted->Draw();
  // p_E->Draw();
  // c_total_energy_distribution->Print(outputdir+Form("/plots/plot_")+version_name+Form("_distribution_")+et_portion+eta_portion+Form("_thr%.1f",threshold_mips)+Form(".png"));
  // c_total_energy_distribution->Print(outputdir+Form("/plots/plot_")+version_name+Form("_distribution_")+et_portion+eta_portion+Form("_thr%.1f",threshold_mips)+Form(".pdf"));
  histograms_output_file->WriteTObject(c_total_energy_distribution);
  // cout << "min layer = " << min_layer << std::endl;
  // cout << "max layer = " << max_layer << std::endl;
  delete gaussian_fit;
  delete h_total_energy_distribution;
  delete c_total_energy_distribution;
  
  if (get_shower_profile) {
    TCanvas *c_shower_profile = new TCanvas(Form("c_shower_profile")+et_portion,Form("shower_profile")+et_portion,800,600);
    TH1F *h_shower_profile = new TH1F(Form("h_shower_profile")+et_portion, Form("_shower_profile")+et_portion+":layer:total deposited energy(weighted)", weights.size(), -0.5, weights.size() - 0.5);
  
    for (unsigned int layer_counter = 0; layer_counter != weights.size(); layer_counter++) {
      h_shower_profile->Fill(layer_counter, total_energies_by_layer[layer_counter]);
    }
    h_shower_profile->Draw();
    histograms_output_file->WriteTObject(c_shower_profile);
    delete h_shower_profile;
    delete c_shower_profile;
  }
  
  if (get_xy_positions_gen_particle) {
    TCanvas *c_xy_positions_gen_particle = new TCanvas(Form("c_xy_positions_gen_particle_")+et_portion,Form("xy_positions_gen_particle_")+et_portion,800,600);
    TH2F *h_xy_positions_gen_particle = new TH2F(Form("h_xy_positions_gen_particle_")+et_portion, Form("xy_positions_gen_particle_")+et_portion+":gen x:gen y", 500, -1000, 1000, 500, -1000, 1000);
    for (std::vector<std::pair<Double_t, Double_t>>::iterator pair_iterator = xy_positions_gen_particle.begin(); pair_iterator != xy_positions_gen_particle.end(); ++pair_iterator) {
      h_xy_positions_gen_particle->Fill((*pair_iterator).first, (*pair_iterator).second);
    }
    h_xy_positions_gen_particle->Draw();
    histograms_output_file->WriteTObject(c_xy_positions_gen_particle);
    delete h_xy_positions_gen_particle;
    delete c_xy_positions_gen_particle;
  }

  if (get_xy_positions_by_layer) {
    for (unsigned int layer_counter = 0; layer_counter != weights.size(); layer_counter++) {
      TCanvas *c_xy_positions = new TCanvas(Form("c_xy_positions_")+et_portion+Form("_layer_%i",layer_counter), Form("xy_positions_")+et_portion+Form("_layer_%i",layer_counter), 800, 600);
      TH2F *h_xy_positions = new TH2F(Form("h_xy_positions_")+et_portion+Form("_layer_%i",layer_counter), Form("xy_positions_")+et_portion+Form("_layer_%i",layer_counter)+":rechit x:rechit y", 500, -1000, 1000, 500, -1000, 1000);
      for (std::vector<std::pair<Double_t, Double_t>>::iterator pair_iterator = xy_positions_by_layer[layer_counter].begin(); pair_iterator != xy_positions_by_layer[layer_counter].end(); ++pair_iterator) {
        h_xy_positions->Fill((*pair_iterator).first, (*pair_iterator).second);
      }
      h_xy_positions->Draw();
      histograms_output_file->WriteTObject(c_xy_positions);
      delete h_xy_positions;
      delete c_xy_positions;
    }
  }
}

void set_eta_portion(Double_t eta) {
  eta_portion = Form("_eta%.3f",eta);
}

void analyze_for_given_eta(Double_t eta) {
  set_eta_portion(eta);
  std::vector<Double_t> energies_incoming_gev;
  std::vector<Double_t> energies_incoming_gev_errors;
  std::vector<Double_t> mean_energies_wmips;
  std::vector<Double_t> mean_energies_wmips_errors;
  std::vector<Double_t> sigmas_wmips;
  std::vector<Double_t> sigmas_wmips_errors;
  std::vector<Int_t> events_measured;
  
  TFile *histograms_output_file = new TFile(outputdir+Form("/root_histograms/histograms_threshold_%.1f_",threshold_mips)+version_name+eta_portion+Form("_results.root"),"RECREATE");

  for (std::vector<Double_t>::iterator et_iterator = et_values.begin(); et_iterator != et_values.end(); ++et_iterator) {
    std::cout << "Analyzing for et = " << *et_iterator << std::endl;
    analyze_for_given_et(*et_iterator, eta, histograms_output_file, energies_incoming_gev, energies_incoming_gev_errors, mean_energies_wmips, mean_energies_wmips_errors, sigmas_wmips, sigmas_wmips_errors, events_measured);
  }
    
  TVectorD Tmean_energies_wmips(mean_energies_wmips.size(),&mean_energies_wmips[0]);
  TVectorD Tmean_energies_wmips_errors(mean_energies_wmips_errors.size(),&mean_energies_wmips_errors[0]);
  // TVectorD Tsigmas_wmips(sigmas_wmips.size(),&sigmas_wmips[0]);
  // TVectorD Tsigmas_wmips_errors(sigmas_wmips_errors.size(),&sigmas_wmips_errors[0]);
  TVectorD Tenergies_incoming_gev(energies_incoming_gev.size(),&energies_incoming_gev[0]);
  TVectorD Tenergies_incoming_gev_errors(energies_incoming_gev_errors.size(),&energies_incoming_gev_errors[0]);
        
  TCanvas *c_mips_calibration = new TCanvas(Form("c_mips_calibration"), Form("mips_calibration"),800,600);
  c_mips_calibration->cd();
  TGraphErrors *energy_mips_calibration = new TGraphErrors(Tenergies_incoming_gev,Tmean_energies_wmips,Tenergies_incoming_gev_errors,Tmean_energies_wmips_errors);
  energy_mips_calibration->Fit("pol1");
  TF1 *calibration_fit = energy_mips_calibration->GetFunction("pol1");
  chisqr = calibration_fit->GetChisquare();
  ndfr = calibration_fit->GetNDF();
  std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
  chisqdf = chisqr/ndfr;
  Double_t offset = calibration_fit->GetParameter(0);
  Double_t offset_error = calibration_fit->GetParError(0);
  Double_t conversion_factor = calibration_fit->GetParameter(1);
  Double_t conversion_factor_error = calibration_fit->GetParError(1);
  TString start_title_conversion_factor = Form("E_mips = (");
  energy_mips_calibration->SetTitle(start_title_conversion_factor + Form("%5.3f +/- %5.3f) * E_GeV + (%.1f +/- %.1f). chisq/dof = %.2f",conversion_factor,conversion_factor_error,offset,offset_error,chisqdf));
  energy_mips_calibration->Draw("AP");
  c_mips_calibration->Update();
  // c_mips_calibration->Print((Form("plot_thr_")+str_threshold+Form("_v")+(Form("%i",version_number)+(Form("_calibration_fit") + eta_portion))) + Form(".pdf"));
  // c_mips_calibration->Print(outputdir+Form("/plots/plot_")+version_name+Form("_calibration_fit")+ eta_portion + Form("_thr%.1f",threshold_mips) + Form(".pdf"));
  histograms_output_file->WriteTObject(c_mips_calibration);
  delete calibration_fit;
  delete energy_mips_calibration;
  delete c_mips_calibration;

  ofstream o_calibration;
  TString o_calibration_name = outputdir+Form("/data/data_calibrations_")+version_name+eta_portion+Form("_thr%.1f",threshold_mips);
  o_calibration.open(o_calibration_name);
  o_calibration << conversion_factor << "    " << conversion_factor_error << "    " << offset << "    " << offset_error << std::endl;
  o_calibration.close();

  std::vector<Double_t> resolutions;
  std::vector<Double_t> resolutions_errors;
  Double_t resolution;
  Double_t resolution_error;
    
  for (unsigned int sigmas_wmips_counter = 0; sigmas_wmips_counter != sigmas_wmips.size(); sigmas_wmips_counter++) {
    resolution = sigmas_wmips[sigmas_wmips_counter]/(mean_energies_wmips[sigmas_wmips_counter]-offset);
    resolution_error = resolution*sqrt(pow(sigmas_wmips_errors[sigmas_wmips_counter]/sigmas_wmips[sigmas_wmips_counter],2)+(pow(mean_energies_wmips_errors[sigmas_wmips_counter],2)+pow(offset_error,2))/pow(mean_energies_wmips[sigmas_wmips_counter]-offset,2));
    resolutions.push_back(resolution);
    resolutions_errors.push_back(resolution_error);
  }
  TVectorD Tresolutions(resolutions.size(),&resolutions[0]);
  TVectorD Tresolutions_errors(resolutions_errors.size(),&resolutions_errors[0]);
  ofstream o_resolutions;
  TString o_resolutions_name = outputdir+Form("/data/data_resolutions_")+version_name+eta_portion+Form("_thr%.1f",threshold_mips);
  o_resolutions.open(o_resolutions_name);
  for(unsigned int et_counter = 0; et_counter != et_values.size(); et_counter++) {
    o_resolutions << et_values[et_counter] << "   " << resolutions[et_counter] << "   " << resolutions_errors[et_counter] << "    " << events_measured[et_counter] << std::endl;
  }
  o_resolutions.close();
    
  TCanvas *c_resolution_versus_energy_standalone_withnoise = new TCanvas(Form("resolution_versus_energy_standalone_withnoise"),"Resolution versus Energy",800,600);
  c_resolution_versus_energy_standalone_withnoise->cd();
  TGraphErrors *resolutions_versus_energy_standalone_withnoise = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors,Tresolutions_errors);
  TF1 *function_to_fit_standalone_withnoise = new TF1("f_to_fit_standalone_withnoise",resolutions_fit_standalone_withnoise,0.0001,1300,3);
  function_to_fit_standalone_withnoise->SetParameters(0.25,0.01,0.00);
  function_to_fit_standalone_withnoise->SetParLimits(2,0.0,0.2);
  function_to_fit_standalone_withnoise->SetParNames("stoch_term","const_term","noise_term");
  resolutions_versus_energy_standalone_withnoise->Fit(function_to_fit_standalone_withnoise,"IME");
  chisqr = function_to_fit_standalone_withnoise->GetChisquare();
  ndfr = function_to_fit_standalone_withnoise->GetNDF();
  std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
  chisqdf = chisqr/ndfr;
  Double_t fit_stoch_withnoise = function_to_fit_standalone_withnoise->GetParameter(0);
  Double_t fit_stoch_error_withnoise = function_to_fit_standalone_withnoise->GetParError(0);
  Double_t fit_const_withnoise = function_to_fit_standalone_withnoise->GetParameter(1);
  Double_t fit_const_error_withnoise = function_to_fit_standalone_withnoise->GetParError(1);
  Double_t fit_noise_withnoise = function_to_fit_standalone_withnoise->GetParameter(2);
  Double_t fit_noise_error_withnoise = function_to_fit_standalone_withnoise->GetParError(2);
  // fit_stoch = function_to_fit_standalone->GetParameter(0);
  // fit_stoch_error = function_to_fit_standalone->GetParError(0);
  // fit_const = function_to_fit_standalone->GetParameter(1);
  // fit_const_error = function_to_fit_standalone->GetParError(1);
  TString start_title_resolutions_withnoise = version_name + eta_portion + Form(": res = (");
  // start_title_resolutions = version_name + eta_portion + Form(": res = (");
  resolutions_versus_energy_standalone_withnoise->SetTitle(start_title_resolutions_withnoise + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV) + (%4.3f +/- %4.3f)/(E/GeV). chisq/dof = %.1f",fit_const_withnoise,fit_const_error_withnoise,fit_stoch_withnoise,fit_stoch_error_withnoise,fit_noise_withnoise,fit_noise_error_withnoise,chisqdf));
  resolutions_versus_energy_standalone_withnoise->Draw("AP");
  c_resolution_versus_energy_standalone_withnoise->Update();
  // c_resolution_versus_energy_standalone_withnoise->Print(outputdir+(Form("/plots/plot_")+(version_name+(Form("_resolutions_versus_energy_standalone_withnoise_") + eta_portion))) + Form("_thr%.1f",threshold_mips) + Form(".pdf"));
  histograms_output_file->WriteTObject(c_resolution_versus_energy_standalone_withnoise);
  ofstream o_fit_withnoise;
  TString o_fit_withnoise_name = outputdir+Form("data/data_fits_withnoise_")+version_name+eta_portion+Form("_thr%.1f",threshold_mips);
  o_fit_withnoise.open(o_fit_withnoise_name);
  o_fit_withnoise << fit_const_withnoise << "    " << fit_const_error_withnoise << "    " << fit_stoch_withnoise << "    " << fit_stoch_error_withnoise << "    " << fit_noise_withnoise << "    " << fit_noise_error_withnoise << "    " << chisqdf << std::endl;
  o_fit_withnoise.close();
  delete resolutions_versus_energy_standalone_withnoise;
  delete function_to_fit_standalone_withnoise;
  delete c_resolution_versus_energy_standalone_withnoise;
  std::cout << std::endl;

  TCanvas *resolution_versus_energy_standalone = new TCanvas(Form("resolution_versus_energy_standalone"),"Resolution versus Energy",800,600);
  resolution_versus_energy_standalone->cd();
  TGraphErrors *resolutions_versus_energy_standalone = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors,Tresolutions_errors);
  TF1 *function_to_fit_standalone = new TF1("f_to_fit_standalone",resolutions_fit_standalone,0.0001,1300,2);
  function_to_fit_standalone->SetParameters(0.25,0.01);
  // function_to_fit_standalone->SetParLimits(2,0.0,0.2);
  function_to_fit_standalone->SetParNames("stoch_term","const_term");
  resolutions_versus_energy_standalone->Fit(function_to_fit_standalone,"IME");
  chisqr = function_to_fit_standalone->GetChisquare();
  ndfr = function_to_fit_standalone->GetNDF();
  std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
  chisqdf = chisqr/ndfr;
  Double_t fit_stoch = function_to_fit_standalone->GetParameter(0);
  Double_t fit_stoch_error = function_to_fit_standalone->GetParError(0);
  Double_t fit_const = function_to_fit_standalone->GetParameter(1);
  Double_t fit_const_error = function_to_fit_standalone->GetParError(1);
  // Double_t fit_noise = function_to_fit_standalone->GetParameter(2);
  // Double_t fit_noise_error = function_to_fit_standalone->GetParError(2);
  // fit_stoch = function_to_fit_standalone->GetParameter(0);
  // fit_stoch_error = function_to_fit_standalone->GetParError(0);
  // fit_const = function_to_fit_standalone->GetParameter(1);
  // fit_const_error = function_to_fit_standalone->GetParError(1);
  TString start_title_resolutions = version_name + eta_portion + Form(": res = (");
  // start_title_resolutions = version_name + eta_portion + Form(": res = (");
  resolutions_versus_energy_standalone->SetTitle(start_title_resolutions + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV). chisq/dof = %.1f",fit_const,fit_const_error,fit_stoch,fit_stoch_error,chisqdf));
  resolutions_versus_energy_standalone->Draw("AP");
  resolution_versus_energy_standalone->Update();
  // resolution_versus_energy_standalone->Print(outputdir+(Form("/plots/plot_")+(version_name+(Form("_resolutions_versus_energy_standalone_") + eta_portion))) + Form("_thr%.1f",threshold_mips) + Form(".pdf"));
  histograms_output_file->WriteTObject(resolution_versus_energy_standalone);
  ofstream o_fit_withoutnoise;
  TString o_fit_withoutnoise_name = outputdir+Form("/data/data_fits_")+version_name+eta_portion+Form("_thr%.1f",threshold_mips);
  o_fit_withoutnoise.open(o_fit_withoutnoise_name);
  o_fit_withoutnoise << fit_const << "    " << fit_const_error << "    " << fit_stoch << "    " << fit_stoch_error << "    " << chisqdf << std::endl;
  o_fit_withoutnoise.close();
  delete resolutions_versus_energy_standalone;
  delete function_to_fit_standalone;
  delete resolution_versus_energy_standalone;
  std::cout << std::endl;

  // myc = new TCanvas("Resolutions versus Energy_standalone_const","Resolution versus Energy",800,600);
  // myc->cd();
  // TGraphErrors *resolutions_versus_energy_standalone_const = new TGraphErrors(Tenergies_incoming_gev,Tresolutions,Tenergies_incoming_gev_errors,Tresolutions_errors);
  // TF1 *function_to_fit_standalone_const = new TF1("f_to_fit_standalone_const",resolutions_fit_standalone_const,0.0001,1300,2);
  // // const Double_t inistoch = 0.25;
  // // const Double_t * pointer_to_inistoch = &inistoch;
  // function_to_fit_standalone_const->SetParameters(0.25,0.0);
  // function_to_fit_standalone_const->SetParLimits(1,0.0,0.2);
  // function_to_fit_standalone_const->SetParNames("stoch_term","noise");
  // resolutions_versus_energy_standalone_const->Fit(function_to_fit_standalone_const,"IME");
  // chisqr = function_to_fit_standalone_const->GetChisquare();
  // ndfr = function_to_fit_standalone_const->GetNDF();
  // std::cout << "checking ... number of ds of f is " << ndfr << std::endl;
  // chisqdf = chisqr/ndfr;
  // fit_stoch = function_to_fit_standalone_const->GetParameter(0);
  // fit_stoch_error = function_to_fit_standalone_const->GetParError(0);
  // fit_const = 0.009;
  // fit_const_error = 0.0;
  // fit_noise = function_to_fit_standalone_const->GetParameter(1);
  // fit_noise_error = function_to_fit_standalone_const->GetParError(1);
  // // fit_stoch = function_to_fit_standalone->GetParameter(0);
  // // fit_stoch_error = function_to_fit_standalone->GetParError(0);
  // // fit_const = function_to_fit_standalone->GetParameter(1);
  // // fit_const_error = function_to_fit_standalone->GetParError(1);
  // start_title_resolutions = version_name + eta_portion + Form(": res = (");
  // // start_title_resolutions = version_name + eta_portion + Form(": res = (");
  // resolutions_versus_energy_standalone_const->SetTitle(start_title_resolutions + Form("%6.5f +/- %6.5f) + (%4.3f +/- %4.3f)/sqrt(E/GeV) + (%4.3f +/- %4.3f)/(E/GeV). chisq/dof = %.1f",fit_const,fit_const_error,fit_stoch,fit_stoch_error,fit_noise,fit_noise_error,chisqdf));
  // resolutions_versus_energy_standalone_const->Draw("AP");
  // myc->Update();
  // myc->Print((Form("plot_")+(version_name+(Form("_resolutions_versus_energy_standalone_const_withnoise_thr%.1f",threshold_mips) + eta_portion))) + Form("_down%.1f_up%.1f_sigreg%i.pdf",sigmas_down,sigmas_up,signal_region));
  // delete resolutions_versus_energy_standalone_const;
  // delete function_to_fit_standalone_const;
  // delete myc;
  // std::cout << std::endl;

  delete histograms_output_file;
}

void check_parameters(Int_t arg_version_number,
                      TString arg_version_name,
                      TString arg_datadir,
                      TString arg_outputdir,
                      Double_t arg_threshold_adc,
                      Int_t arg_digi_or_raw_switch,
                      Bool_t arg_use_modified_weighting_scheme,
                      Bool_t arg_multiple_runs_present,
                      Bool_t arg_get_hits_distribution,
                      Bool_t arg_get_shower_profile,
                      Bool_t arg_get_xy_positions_gen_particle,
                      Bool_t arg_get_xy_positions_by_layer,
                      unsigned arg_max_events
                      ) {
  std::cout << "Running analysis for version = " << arg_version_number << ";" << std::endl
            << "Data directory = " << arg_datadir << ";" << std::endl
            << "Output directory = " << arg_outputdir << ";" << std::endl
            << "Threshold = " << arg_threshold_adc << " ADCs;" << std::endl
            << "Switch = " << (arg_digi_or_raw_switch == 1? "Digi;" : (arg_digi_or_raw_switch == 2? "Raw;" : (arg_digi_or_raw_switch == 3? "Digi+Raw;" : "None;"))) << std::endl
            << (arg_use_modified_weighting_scheme? "Using modified weighting scheme;" : "Using normal weighting scheme;") << std::endl
            << (arg_multiple_runs_present? "Running analysis for multiple runs;" : "Running analysis for single runs;") << std::endl
            << (arg_get_hits_distribution? "Getting distribution of individual hits;" : "Not getting distribution of individual hits;") << std::endl
            << (arg_get_shower_profile? "Getting shower profile;" : "Not getting shower profile;") << std::endl
            << (arg_get_xy_positions_gen_particle? "Getting x,y distribution of generated particles" : "Not getting x,y distribution of generated particles") << std::endl
            << (arg_get_xy_positions_by_layer? "Getting x,y distribution of rechits" : "Not getting x,y distribution of rechits") << std::endl;

  if (arg_max_events == 0) {
    std::cout << "Analyzing all events" << std::endl;
  }
  else {
    std::cout << "Analyzing at most " << arg_max_events << " events" << std::endl;
  }
  
  if (arg_threshold_adc < 0.5 && arg_digi_or_raw_switch != 2) {
    std::cout << "Threshold provided is less than minimum threshold analyzable from Digi file" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else if (arg_digi_or_raw_switch != 1 && arg_digi_or_raw_switch != 2 && arg_digi_or_raw_switch != 3) {
    std::cout << "Digi hits or raw hits?" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else if (arg_get_xy_positions_gen_particle && (arg_digi_or_raw_switch == 1)) {
    std::cout << "Generated particle position cannot be found from Digis alone" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else if (arg_get_xy_positions_by_layer && (arg_digi_or_raw_switch == 2)) {
    std::cout << "Cannot yet find xy positions for simhits" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  else if (arg_threshold_adc > 0 && arg_digi_or_raw_switch == 2) {
    std::cout << "If you are sure you want to proceed with a nonzero threshold for HGCal hits, please comment this section of the code line" << std::endl;
    std::exit(EXIT_FAILURE);
  }
}

void set_parameters(Int_t arg_version_number,
                    TString arg_version_name,
                    TString arg_datadir,
                    TString arg_outputdir,
                    Double_t arg_threshold_adc,
                    Int_t arg_digi_or_raw_switch,
                    Bool_t arg_use_modified_weighting_scheme,
                    Bool_t arg_multiple_runs_present,
                    Bool_t arg_get_hits_distribution,
                    Bool_t arg_get_shower_profile,
                    Bool_t arg_get_xy_positions_gen_particle,
                    Bool_t arg_get_xy_positions_by_layer,
                    unsigned arg_max_events
                    ) {
  version_number = arg_version_number;
  version_name = arg_version_name;
  datadir = arg_datadir;
  outputdir = arg_outputdir;
  threshold_adc = arg_threshold_adc;
  digi_or_raw_switch = arg_digi_or_raw_switch;
  use_modified_weighting_scheme = arg_use_modified_weighting_scheme;
  multiple_runs_present = arg_multiple_runs_present;
  get_hits_distribution = arg_get_hits_distribution;
  get_shower_profile = arg_get_shower_profile;
  get_xy_positions_gen_particle = arg_get_xy_positions_gen_particle;
  get_xy_positions_by_layer = arg_get_xy_positions_by_layer;
  max_events = arg_max_events;

  threshold_mips = threshold_adc/10.0;
}

void run_analysis(Int_t arg_version_number,
                  TString arg_version_name,
                  TString arg_datadir,
                  TString arg_outputdir,
                  Double_t arg_threshold_adc,
                  Int_t arg_digi_or_raw_switch,
                  Bool_t arg_use_modified_weighting_scheme,
                  Bool_t arg_multiple_runs_present,
                  Bool_t arg_get_hits_distribution,
                  Bool_t arg_get_shower_profile,
                  Bool_t arg_get_xy_positions_gen_particle,
                  Bool_t arg_get_xy_positions_by_layer,
                  unsigned arg_max_events
                  ) { // main
  
  // load the shared library for HGCSS* classes:
  gSystem->Load("/afs/cern.ch/user/t/tmudholk/public/research/hgcal_optimization_latest/PFCal/PFCalEE/userlib/lib/libPFCalEEuserlib.so");

  check_parameters(arg_version_number, arg_version_name, arg_datadir, arg_outputdir, arg_threshold_adc, arg_digi_or_raw_switch, arg_use_modified_weighting_scheme, arg_multiple_runs_present, arg_get_hits_distribution, arg_get_shower_profile, arg_get_xy_positions_gen_particle, arg_get_xy_positions_by_layer, arg_max_events);

  set_parameters(arg_version_number, arg_version_name, arg_datadir, arg_outputdir, arg_threshold_adc, arg_digi_or_raw_switch, arg_use_modified_weighting_scheme, arg_multiple_runs_present, arg_get_hits_distribution, arg_get_shower_profile, arg_get_xy_positions_gen_particle, arg_get_xy_positions_by_layer, arg_max_events);

  // TString str_threshold = Form("%.2f",threshold);
  // TString HGcal_common_prefix_firstpart = datadir + Form("/HGcal__version%i_model2_BOFF_",version_number);

  read_weights(version_number, datadir, use_modified_weighting_scheme);
    
  for (std::vector<Double_t>::iterator eta_iterator = eta_values.begin(); eta_iterator != eta_values.end(); ++eta_iterator) {
    std::cout << "Analyzing for eta = " << *eta_iterator << std::endl;
    analyze_for_given_eta(*eta_iterator);
  }
  // return 0;
} // end main
