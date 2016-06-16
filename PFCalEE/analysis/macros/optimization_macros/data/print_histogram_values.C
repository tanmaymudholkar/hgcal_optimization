#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<stdlib.h>

#include "TH1F.h"
#include "TFile.h"

void print_histogram_values_to_file(TH1F *input_histogram, TString filename) {
  unsigned nbins_hist = input_histogram->GetNbinsX();
  ofstream o_histogram_contents;
  o_histogram_contents.open(filename);
  for (unsigned counter_bins = 1; counter_bins != 1 + nbins_hist; ++counter_bins) {
    o_histogram_contents << input_histogram->GetXaxis()->GetBinCenter(counter_bins) << "    " << input_histogram->GetBinContent(counter_bins) << std::endl;
  }
  o_histogram_contents.close();
}

void print_histogram_values() {
  Int_t version_number = 30;
  TString version_name = "suppression_study";
  Double_t thresholds_to_compare_array[] = {5.0, 3.5, 2.0, 1.0};
  std::vector<Double_t> thresholds_to_compare(thresholds_to_compare_array, thresholds_to_compare_array+sizeof(thresholds_to_compare_array)/sizeof(Double_t));
  Double_t et_values_array[] = {3,5,10,30,50,70,100,150}; // version 33 and 30, hexagonal geometry
  std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
  TFile *f_histograms_source = new TFile("fractional_difference_histograms.root");

  for (std::vector<Double_t>::iterator threshold_iterator = thresholds_to_compare.begin(); threshold_iterator != thresholds_to_compare.end(); ++threshold_iterator) {
    for (std::vector<Double_t>::iterator et_iterator = et_values.begin(); et_iterator != et_values.end(); ++et_iterator) {
      TString histogram_name = Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", *et_iterator, *threshold_iterator);
      TString filename = Form("analyzed_data_histogram_contents_version%i_", version_number) + version_name + Form("_et%.0f_thr%.1f", *et_iterator, *threshold_iterator);
      std::cout << "histogram name: " << histogram_name << "    filename: " << filename << std::endl;
      TH1F *h_to_print = (TH1F*)f_histograms_source->Get(histogram_name);
      print_histogram_values_to_file(h_to_print, filename);
    }
  }
}
