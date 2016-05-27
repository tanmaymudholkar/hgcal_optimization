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
// #include "Riostream.h"
// #include "TString.h"

void analyze_data() {
  //  Read data from an ascii file and create a root file with an histogram and an ntuple.
  //   see a variant of this macro in basic2.C
  //Author: Rene Brun


  // read file $ROOTSYS/tutorials/tree/basic.dat
  // this file has 3 columns of float data
  Double_t et_values_array[] = {3,5,10,30,50,70,100,150}; // version 33 and 30, hexagonal geometry
  Double_t base_threshold = 0.5;
  Double_t threshold_to_compare = 5.0;
  Int_t version_number = 30;
  // Double_t thresholds_to_compare_array[] = {2.0, 5.0};
  std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
  // std::vector<Double_t> thresholds_to_compare(thresholds_to_compare_array, thresholds_to_compare_array+sizeof(thresholds_to_compare_array)/sizeof(Double_t));

  std::vector<Double_t> means;
  std::vector<Double_t> rmses;
  TFile *f = new TFile("percent_difference_histograms.root","RECREATE");

  for (std::vector<Double_t>::iterator et = et_values.begin(); et != et_values.end(); ++et) {
    TString et_portion = Form("et%.0f",*et);
    std::cout << "Et = " << *et << std::endl;
    ifstream in_bare;
    // TString filename_in_bare = Form("data_total_energies_version%i_half%s_eta2.100_thr%.1f", version_number, et_portion, base_threshold);
    // std::cout << "filename in bare: " << filename_in_bare << std::endl;
    in_bare.open(Form("data_total_energies_version%i_half", version_number) + et_portion + Form("_eta2.100_thr%.1f", base_threshold));
    // in_bare.open("data_total_energies_version30_halfet10_eta2.100_thr0.5");
    ifstream in_to_compare;
    in_to_compare.open(Form("data_total_energies_version%i_half", version_number) + et_portion + Form("_eta2.100_thr%.1f", threshold_to_compare));
    // in_to_compare.open(Form("data_total_energies_version%i_half%s_eta2.100_thr%.1f", version_number, et_portion, threshold_to_compare));
    // in_to_compare.open("data_total_energies_version30_halfet10_eta2.100_thr5.0");
    Int_t evt_number;
    Int_t evt_number_to_compare;
    Double_t energy_deposited_bare;
    Double_t energy_deposited_to_compare;
    std::vector<Double_t> percent_difference;
    

    while (1) {
      in_bare >> evt_number >> energy_deposited_bare;
      if (!in_bare.good()) break;
      in_to_compare >> evt_number_to_compare >> energy_deposited_to_compare;
      if (!in_to_compare.good()) break;
      if (evt_number < 50 or evt_number > 10150) {
        printf("evt_no_1 = %i, evt_no_2 = %i, energy_deposited_bare=%.2f, energy_deposited_to_compare=%.2f\n",evt_number,evt_number_to_compare,energy_deposited_bare, energy_deposited_to_compare);
      }
      percent_difference.push_back((energy_deposited_bare - energy_deposited_to_compare)/energy_deposited_bare);
    }
    in_bare.close();
    in_to_compare.close();
    TH1F *h1 = new TH1F(Form("percent difference distribution ")+et_portion, "percent difference distribution",50,0,0);
    for (std::vector<Double_t>::iterator percent_difference_iterator = percent_difference.begin(); percent_difference_iterator != percent_difference.end(); ++percent_difference_iterator) {
      h1->Fill(*percent_difference_iterator);
    }
    means.push_back(h1->GetMean());
    rmses.push_back(h1->GetRMS());
    f->WriteTObject(h1);
    delete h1;
  }

  for (unsigned int counter = 0; counter < means.size(); ++counter) {
    std::cout << et_values[counter] << "    " << means[counter] << "    " << rmses[counter] << std::endl;
  }
  
  // TString dir = gSystem->UnixPathName(__FILE__);
  // dir.ReplaceAll("analyze_data.C","");
  // dir.ReplaceAll("/./","/");
  

  // Float_t x,y,z;
  // Int_t nlines = 0;
  // TFile *f = new TFile("basic.root","RECREATE");
  // TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
  // TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y:z");

  // while (1) {
  //   in >> x >> y >> z;
  //   if (!in.good()) break;
  //   if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
  //   h1->Fill(x);
  //   ntuple->Fill(x,y,z);
  //   nlines++;
  // }
  // printf(" found %d points\n",nlines);

  // in.close();

  // f->Write();
}
