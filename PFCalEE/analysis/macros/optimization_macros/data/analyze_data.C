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

  Double_t et_values_array[] = {3,5,10,30,50,70,100,150}; // version 33 and 30, hexagonal geometry
  // Double_t ranges_array[] = {0.25, 0.195, 0.15, 0.09, 0.07, 0.06, 0.05, 0.041};
  Int_t version_number = 30;
  TString version_name = "suppression_study";
  Double_t thresholds_to_compare_array[] = {1.0, 2.0, 3.5, 5.0, 7.5, 10.0};
  Double_t base_threshold = 0.5;
  // Double_t max_threshold = 10.0;
  Double_t number_of_sigmas_for_max_range = 3.6;
  std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
  // std::vector<Double_t> ranges(ranges_array, ranges_array+sizeof(ranges_array)/sizeof(Double_t));
  std::vector<Double_t> range_maxvals;
  std::vector<Double_t> thresholds_to_compare(thresholds_to_compare_array, thresholds_to_compare_array+sizeof(thresholds_to_compare_array)/sizeof(Double_t));
  Double_t max_threshold = thresholds_to_compare[thresholds_to_compare.size()-1];

  // First, calculate means and RMSes
  for (std::vector<Double_t>::iterator thresholds_to_compare_iterator = thresholds_to_compare.begin(); thresholds_to_compare_iterator != thresholds_to_compare.end(); ++thresholds_to_compare_iterator) {
    std::vector<Double_t> means;
    std::vector<Double_t> means_errors;
    std::vector<Double_t> rmses;
    std::vector<Double_t> rmses_errors;
    // for (std::vector<Double_t>::iterator et_iterator = et_values.begin(); et_iterator != et_values.end(); ++et_iterator) {
    for (unsigned int et_counter = 0; et_counter != et_values.size(); ++et_counter) {
      Double_t et = et_values[et_counter];
      TString et_portion = Form("et%.0f",et);
      // std::cout << "Et = " << et << std::endl;
      TString in_base_filename = Form("data_total_energies_version%i_suppression_study", version_number) + et_portion + Form("_eta2.100_thr%.1f", base_threshold);
      // TCanvas c_fractional_difference_histograms = new TCanvas(Form("c_fractional_difference_histograms_")+et_portion, Form("c_fractional_difference_histograms_")+et_portion, 800, 600);
      ifstream in_base;
      // std::cout << "in_base_filename: " << in_base_filename << std::endl;
      TString in_to_compare_filename = Form("data_total_energies_version%i_suppression_study", version_number) + et_portion + Form("_eta2.100_thr%.1f", *thresholds_to_compare_iterator);
      in_base.open(in_base_filename);
      // in_base.open("data_total_energies_version30_halfet10_eta2.100_thr0.5");
      ifstream in_to_compare;
      // in_to_compare.open(Form("data_total_energies_version%i_half", version_number) + et_portion + Form("_eta2.100_thr%.1f", threshold_to_compare));
      in_to_compare.open(in_to_compare_filename);
      // in_to_compare.open(Form("data_total_energies_version%i_half%s_eta2.100_thr%.1f", version_number, et_portion, threshold_to_compare));
      // in_to_compare.open("data_total_energies_version30_halfet10_eta2.100_thr5.0");
      Int_t evt_number;
      Int_t evt_number_to_compare;
      Double_t energy_deposited_base;
      Double_t energy_deposited_to_compare;
      std::vector<Double_t> fractional_difference;

      while (true) {
        in_base >> evt_number >> energy_deposited_base;
        if (!in_base.good()) break;
        in_to_compare >> evt_number_to_compare >> energy_deposited_to_compare;
        if (!in_to_compare.good()) break;
        // if (evt_number < 50 or evt_number > 10150) {
        //   printf("evt_no_1 = %i, evt_no_2 = %i, energy_deposited_base=%.2f, energy_deposited_to_compare=%.2f\n",evt_number,evt_number_to_compare,energy_deposited_base, energy_deposited_to_compare);
        // }
        if (evt_number != evt_number_to_compare) {
          std::cout << "Problem! Different event numbers?" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        fractional_difference.push_back((energy_deposited_base - energy_deposited_to_compare)/energy_deposited_base);
      }
      in_base.close();
      in_to_compare.close();
      TH1F *h1 = new TH1F(Form("h_fractional_difference_distribution_") + et_portion + Form("_eta2.100_thr%.1f", *thresholds_to_compare_iterator), Form("h_fractional_difference_distribution_") + et_portion + Form("_eta2.100_thr%.1f", *thresholds_to_compare_iterator), 100, 0, 0);
      for (std::vector<Double_t>::iterator fractional_difference_iterator = fractional_difference.begin(); fractional_difference_iterator != fractional_difference.end(); ++fractional_difference_iterator) {
        h1->Fill(*fractional_difference_iterator);
      }

      means.push_back(h1->GetMean());
      means_errors.push_back(h1->GetMeanError());
      rmses.push_back(h1->GetRMS());
      rmses_errors.push_back(h1->GetRMSError());
      
      // f->WriteTObject(h1);
      delete h1;
    } // Ends loop over energies

    ofstream o_means_and_rmses;
    TString o_means_and_rmses_name = Form("analyzed_data_means_and_rmses_version%i_", version_number) + version_name + Form("_thr%.1f",*thresholds_to_compare_iterator);
    o_means_and_rmses.open(o_means_and_rmses_name);
    // o_means_and_rmses << conversion_factor << "    " << conversion_factor_error << "    " << offset << "    " << offset_error << std::endl;
    for (unsigned int counter = 0; counter < means.size(); ++counter) {
      o_means_and_rmses << et_values[counter] << "    " << means[counter] << "    " << means_errors[counter] << "    " << rmses[counter] << "    " << rmses_errors[counter] << std::endl;
    } // Output means and rmses
    o_means_and_rmses.close();
  } // Ends loop over thresholds to compare

  // Next, get maximum values for x-axis
  for (std::vector<Double_t>::iterator et_iterator = et_values.begin(); et_iterator != et_values.end(); ++et_iterator) {
    TString et_portion = Form("et%.0f",*et_iterator);
    TString in_base_for_maxval_filename = Form("data_total_energies_version%i_suppression_study", version_number) + et_portion + Form("_eta2.100_thr%.1f", base_threshold);
    ifstream in_base_for_maxval;
    TString in_to_compare_for_maxval_filename = Form("data_total_energies_version%i_suppression_study", version_number) + et_portion + Form("_eta2.100_thr%.1f", max_threshold);
    ifstream in_to_compare_for_maxval;
    Int_t evt_number;
    Int_t evt_number_to_compare;
    Double_t energy_deposited_base;
    Double_t energy_deposited_to_compare;
    Double_t fdiff_for_maxval;

    Double_t maxval = 0;

    in_base_for_maxval.open(in_base_for_maxval_filename);
    in_to_compare_for_maxval.open(in_to_compare_for_maxval_filename);
    TH1F *histogram_for_maxval = new TH1F("histogram_for_maxval", "histogram_for_maxval", 100, 0, 0);

    while (true) {
      in_base_for_maxval >> evt_number >> energy_deposited_base;
      if (!in_base_for_maxval.good()) break;
      in_to_compare_for_maxval >> evt_number_to_compare >> energy_deposited_to_compare;
      if (!in_to_compare_for_maxval.good()) break;
      // if (evt_number < 50 or evt_number > 10150) {
      //   printf("evt_no_1 = %i, evt_no_2 = %i, energy_deposited_base=%.2f, energy_deposited_to_compare=%.2f\n",evt_number,evt_number_to_compare,energy_deposited_base, energy_deposited_to_compare);
      // }
      if (evt_number != evt_number_to_compare) {
        std::cout << "Problem! Different event numbers?" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      fdiff_for_maxval = (energy_deposited_base - energy_deposited_to_compare)/energy_deposited_base;
      // if (fdiff_for_maxval > maxval) {
      //   maxval = fdiff_for_maxval;
      // }
      histogram_for_maxval->Fill(fdiff_for_maxval);
    }
    // std::cout << "Maximum value for et = " << *et_iterator << ": " << maxval << std::endl;
    Double_t mean_for_maxval = histogram_for_maxval->GetMean();
    // Double_t mean_for_maxval_error = histogram_for_maxval->GetMeanError();
    Double_t sigma_for_maxval = histogram_for_maxval->GetRMS();
    // Double_t sigma_for_maxval_error = histogram_for_maxval->GetRMSError();

    std::cout << "Mean for maxval: " << mean_for_maxval << "   sigma for maxval: " << sigma_for_maxval << std::endl;

    // Double_t sigmas_comparison_array[] = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    // std::vector<Double_t> sigmas_comparison(sigmas_comparison_array, sigmas_comparison_array+sizeof(sigmas_comparison_array)/sizeof(Double_t));
    // for (std::vector<Double_t>::iterator sigmas_comparison_iterator = sigmas_comparison.begin(); sigmas_comparison_iterator != sigmas_comparison.end(); ++sigmas_comparison_iterator) {
    //   std::cout << "mean: " << mean_for_comparison << "  mean + " << *sigmas_comparison_iterator << "*sigma = " << mean_for_comparison + (*sigmas_comparison_iterator)*sigma_for_comparison << std::endl;
    // }
    range_maxvals.push_back(mean_for_maxval + number_of_sigmas_for_max_range*sigma_for_maxval);
    
    delete histogram_for_maxval;
  }

  // Finally, plot the histograms using the max values calculated earlier
  TFile *f = new TFile("fractional_difference_histograms.root","RECREATE");
  
  for (std::vector<Double_t>::iterator thresholds_to_compare_iterator = thresholds_to_compare.begin(); thresholds_to_compare_iterator != thresholds_to_compare.end(); ++thresholds_to_compare_iterator) {
    // std::vector<Double_t> means;
    // std::vector<Double_t> means_errors;
    // std::vector<Double_t> rmses;
    // std::vector<Double_t> rmses_errors;
    // for (std::vector<Double_t>::iterator et_iterator = et_values.begin(); et_iterator != et_values.end(); ++et_iterator) {
    for (unsigned int et_counter = 0; et_counter != et_values.size(); ++et_counter) {
      Double_t et = et_values[et_counter];
      TString et_portion = Form("et%.0f",et);
      // std::cout << "Et = " << et << std::endl;
      TString in_base_filename = Form("data_total_energies_version%i_suppression_study", version_number) + et_portion + Form("_eta2.100_thr%.1f", base_threshold);
      // TCanvas c_fractional_difference_histograms = new TCanvas(Form("c_fractional_difference_histograms_")+et_portion, Form("c_fractional_difference_histograms_")+et_portion, 800, 600);
      ifstream in_base;
      // std::cout << "in_base_filename: " << in_base_filename << std::endl;
      TString in_to_compare_filename = Form("data_total_energies_version%i_suppression_study", version_number) + et_portion + Form("_eta2.100_thr%.1f", *thresholds_to_compare_iterator);
      in_base.open(in_base_filename);
      // in_base.open("data_total_energies_version30_halfet10_eta2.100_thr0.5");
      ifstream in_to_compare;
      // in_to_compare.open(Form("data_total_energies_version%i_half", version_number) + et_portion + Form("_eta2.100_thr%.1f", threshold_to_compare));
      in_to_compare.open(in_to_compare_filename);
      // in_to_compare.open(Form("data_total_energies_version%i_half%s_eta2.100_thr%.1f", version_number, et_portion, threshold_to_compare));
      // in_to_compare.open("data_total_energies_version30_halfet10_eta2.100_thr5.0");
      Int_t evt_number;
      Int_t evt_number_to_compare;
      Double_t energy_deposited_base;
      Double_t energy_deposited_to_compare;
      std::vector<Double_t> fractional_difference;

      while (true) {
        in_base >> evt_number >> energy_deposited_base;
        if (!in_base.good()) break;
        in_to_compare >> evt_number_to_compare >> energy_deposited_to_compare;
        if (!in_to_compare.good()) break;
        // if (evt_number < 50 or evt_number > 10150) {
        //   printf("evt_no_1 = %i, evt_no_2 = %i, energy_deposited_base=%.2f, energy_deposited_to_compare=%.2f\n",evt_number,evt_number_to_compare,energy_deposited_base, energy_deposited_to_compare);
        // }
        if (evt_number != evt_number_to_compare) {
          std::cout << "Problem! Different event numbers?" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        fractional_difference.push_back((energy_deposited_base - energy_deposited_to_compare)/energy_deposited_base);
      }
      in_base.close();
      in_to_compare.close();
      TH1F *h1 = new TH1F(Form("h_fractional_difference_distribution_") + et_portion + Form("_eta2.100_thr%.1f", *thresholds_to_compare_iterator), Form("h_fractional_difference_distribution_") + et_portion + Form("_eta2.100_thr%.1f", *thresholds_to_compare_iterator), 1000, 0, range_maxvals[et_counter]);
      for (std::vector<Double_t>::iterator fractional_difference_iterator = fractional_difference.begin(); fractional_difference_iterator != fractional_difference.end(); ++fractional_difference_iterator) {
        h1->Fill(*fractional_difference_iterator);
      }

      // means.push_back(h1->GetMean());
      // means_errors.push_back(h1->GetMeanError());
      // rmses.push_back(h1->GetRMS());
      // rmses_errors.push_back(h1->GetRMSError());
      
      f->WriteTObject(h1);
      delete h1;
    } // Ends loop over energies

    // ofstream o_means_and_rmses;
    // TString o_means_and_rmses_name = Form("analyzed_data_means_and_rmses_version%i_", version_number) + version_name + Form("_thr%.1f",*thresholds_to_compare_iterator);
    // o_means_and_rmses.open(o_means_and_rmses_name);
    // // o_means_and_rmses << conversion_factor << "    " << conversion_factor_error << "    " << offset << "    " << offset_error << std::endl;
    // for (unsigned int counter = 0; counter < means.size(); ++counter) {
    //   o_means_and_rmses << et_values[counter] << "    " << means[counter] << "    " << means_errors[counter] << "    " << rmses[counter] << "    " << rmses_errors[counter] << std::endl;
    // } // Output means and rmses
    // o_means_and_rmses.close();
  } // Ends loop over thresholds to compare
  
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
