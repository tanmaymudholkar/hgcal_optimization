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
#include "TColor.h"
#include "TLegend.h"

const Double_t scaling_factor = 1.03;
const Int_t loop_over_these_colors_array[] = {kRed, kGreen, kBlue, kYellow, kCyan};
unsigned int size_of_colors_to_loop_over = sizeof(loop_over_these_colors_array)/sizeof(Int_t);

void generate_overlays_given_et(TFile *f_histograms_source, TCanvas *current_canvas, Int_t pad_position, Double_t base_threshold, Double_t threshold_for_maxval, std::vector<Double_t> thresholds_to_compare, Double_t et) {

  TH1F *histogram_for_maxval = (TH1F*)(f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", et, threshold_for_maxval)));
  Double_t maxval = histogram_for_maxval->GetMaximum();
  maxval = scaling_factor * maxval;

  std::cout << "Maxval: " << maxval << std::endl;

  // TCanvas *c_overlaid_histograms = new TCanvas("overlaid_histograms", "overlaid_histograms", 1024, 768);
  current_canvas->cd(pad_position);
  gStyle->SetOptStat(0);
  TLegend *legend = new TLegend(0.4, 0.5, 0.9, 0.9, "Color Scheme for Thresholds");
  TH1F *base_histogram = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", et, base_threshold));
  legend->AddEntry(base_histogram, Form("%.1f", base_threshold));
  base_histogram->SetTitle(Form("E_{T} = %.0f GeV; fractional difference", et));
  base_histogram->SetMaximum(maxval);
  base_histogram->SetLineColor(kBlack);
  base_histogram->Draw();
  current_canvas->Update();
  // for (std::vector<Double_t>::iterator thresholds_to_compare_iterator = thresholds_to_compare.begin(); thresholds_to_compare_iterator != thresholds_to_compare.end(); ++thresholds_to_compare_iterator) {
  for (unsigned int thresholds_counter = 0; thresholds_counter != thresholds_to_compare.size(); ++thresholds_counter) {
    TH1F *histogram_to_compare = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", et, thresholds_to_compare[thresholds_counter]));
    legend->AddEntry(histogram_to_compare, Form("%.1f", thresholds_to_compare[thresholds_counter]));
    histogram_to_compare->SetMaximum(maxval);
    // base_histogram->SetLineColor(loop_over_these_colors_array[thresholds_counter%size_of_colors_to_loop_over]);
    histogram_to_compare->SetLineColor(loop_over_these_colors_array[thresholds_counter%size_of_colors_to_loop_over]);
    histogram_to_compare->Draw("same");
    current_canvas->Update();
  }
  legend->SetTextSize(0.05);
  legend->Draw();
  current_canvas->Update();
  // c_overlaid_histograms->Print(Form("overlaid_histograms_et%.0f.png", et));
  // delete c_overlaid_histograms;
}

void generate_histogram_overlays() {
  Double_t et_values_array[] = {3,5,10,30,50,70,100,150}; // version 33 and 30, hexagonal geometry
  std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
  Double_t base_threshold = 10.0;
  // Int_t version_number = 30;
  // TString version_name = "suppression_study";
  // Double_t thresholds_to_compare_array[] = {1.0, 2.0, 3.5};
  Double_t thresholds_to_compare_array[] = {7.5, 5.0, 3.5, 2.0, 1.0};
  std::vector<Double_t> thresholds_to_compare(thresholds_to_compare_array, thresholds_to_compare_array+sizeof(thresholds_to_compare_array)/sizeof(Double_t));
  Double_t threshold_for_maxval = thresholds_to_compare[(sizeof(thresholds_to_compare_array)/sizeof(Double_t))-1];
  
  TFile *f_histograms_source = new TFile("fractional_difference_histograms.root");
  TCanvas *c_first_four = new TCanvas("c_first_four", "c_first_four", 1024, 768);
  c_first_four->Divide(2,2);
  generate_overlays_given_et(f_histograms_source, c_first_four, 1, base_threshold, threshold_for_maxval, thresholds_to_compare, 3.0);
  generate_overlays_given_et(f_histograms_source, c_first_four, 2, base_threshold, threshold_for_maxval, thresholds_to_compare, 5.0);
  generate_overlays_given_et(f_histograms_source, c_first_four, 3, base_threshold, threshold_for_maxval, thresholds_to_compare, 10.0);
  generate_overlays_given_et(f_histograms_source, c_first_four, 4, base_threshold, threshold_for_maxval, thresholds_to_compare, 30.0);
  c_first_four->Print("et_3_5_10_30.png");
  delete c_first_four;

  TCanvas *c_next_four = new TCanvas("c_next_four", "c_next_four", 1024, 768);
  c_next_four->Divide(2,2);
  generate_overlays_given_et(f_histograms_source, c_next_four, 1, base_threshold, threshold_for_maxval, thresholds_to_compare, 50.0);
  generate_overlays_given_et(f_histograms_source, c_next_four, 2, base_threshold, threshold_for_maxval, thresholds_to_compare, 70.0);
  generate_overlays_given_et(f_histograms_source, c_next_four, 3, base_threshold, threshold_for_maxval, thresholds_to_compare, 100.0);
  generate_overlays_given_et(f_histograms_source, c_next_four, 4, base_threshold, threshold_for_maxval, thresholds_to_compare, 150.0);
  c_next_four->Print("et_50_70_100_150.png");
  delete c_next_four;
  
  delete f_histograms_source;
  
  
  // for (std::vector<Double_t>::iterator et_iterator = et_values.begin(); et_iterator != et_values.end(); ++et_iterator) {
  //   TH1F *histogram_for_maxval = (TH1F*)(f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", *et_iterator, threshold_for_maxval)));
  //   Double_t maxval = histogram_for_maxval->GetMaximum();
  //   maxval = scaling_factor * maxval;

  //   std::cout << "Maxval: " << maxval << std::endl;

  //   TCanvas *c_overlaid_histograms = new TCanvas("overlaid_histograms", "overlaid_histograms", 1024, 768);
  //   c_overlaid_histograms->cd();
  //   TH1F *base_histogram = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", *et_iterator, base_threshold));
  //   base_histogram->SetMaximum(maxval);
  //   base_histogram->Draw();
  //   c_overlaid_histograms->Update();
  //   for (std::vector<Double_t>::iterator thresholds_to_compare_iterator = thresholds_to_compare.begin(); thresholds_to_compare_iterator != thresholds_to_compare.end(); ++thresholds_to_compare_iterator) {
  //     TH1F *histogram_to_compare = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", *et_iterator, *thresholds_to_compare_iterator));
  //     histogram_to_compare->SetMaximum(maxval);
  //     histogram_to_compare->Draw("same");
  //     c_overlaid_histograms->Update();
  //   }
  //   c_overlaid_histograms->Print(Form("overlaid_histograms_et%.0f.png", *et_iterator));
  //   delete c_overlaid_histograms;
  // }
}



// void generate_overlays_given_et(TFile *f_histograms_source, TCanvas *current_canvas, Int_t pad_position, Double_t base_threshold, std::vector<Double_t> thresholds_to_compare, Double_t et) {
//   // std::vector<TH1F*> histograms;

//   current_canvas->cd(pad_position);

//   TH1F *base_histogram = new TH1F(Form("h_thr%.1f", base_threshold), Form("h_thr%.1f", base_threshold), 100, 0, 0);
//   base_histogram = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", et, base_threshold));
//   Double_t maxval = 0;

//   for (unsigned counter_to_compare = 0; counter_to_compare != thresholds_to_compare.size(); ++counter_to_compare) {
//     TH1F *histogram_to_compare = new TH1F(Form("h_thr%.1f", thresholds_to_compare[counter_to_compare]), Form("h_thr%.1f", thresholds_to_compare[counter_to_compare]), 100, 0, 0);
//     histogram_to_compare = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", et, thresholds_to_compare[counter_to_compare]));
//     if (counter_to_compare == 0) {
//       // maxval = histogram_to_compare->GetMaximum();
//       // base_histogram->SetMaximum(1.1*maxval);
//       base_histogram->Draw();
//     }
//     // histogram_to_compare->SetMaximum(1.1*maxval);
//     histogram_to_compare->Draw("same");
//   }
  
//   // Double_t maxval = histograms[1]->GetMaximum();
//   // histograms[0]->SetMaximum(1.1*maxval);
//   // histograms[0]->Draw();
//   // for (unsigned histograms_counter = 1; histograms_counter != histograms.size(); ++histograms_counter) {
//   //   histograms[histograms_counter]->SetMaximum(1.1*maxval);
//   //   histograms[histograms_counter]->Draw("same");
//   // }
// }

// void generate_histogram_overlays() {
//   Double_t et_values_array[] = {3,5,10,30,50,70,100,150}; // version 33 and 30, hexagonal geometry
//   Double_t base_threshold = 5.0;
//   Int_t version_number = 30;
//   TString version_name = "suppression_study";
//   Double_t thresholds_to_compare_array[] = {0.5, 1.0, 2.0, 3.5};
//   std::vector<Double_t> et_values(et_values_array,et_values_array+sizeof(et_values_array)/sizeof(Double_t));
//   std::vector<Double_t> thresholds_to_compare(thresholds_to_compare_array, thresholds_to_compare_array+sizeof(thresholds_to_compare_array)/sizeof(Double_t));

//   TFile *f_histograms_source = new TFile("fractional_difference_histograms.root");

//   TCanvas *first_four= new TCanvas("Et 3, 5, 10, 30", "Et 3, 5, 10, 30", 1024, 768);
//   first_four->Divide(2,2);
//   // Double_t et = 3.0;

//   first_four->cd();

//   TH1F *base_histogram = new TH1F(Form("h_thr%.1f", base_threshold), Form("h_thr%.1f", base_threshold), 100, 0, 0);
//   base_histogram = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", 3.0, base_threshold));
//   base_histogram->Draw();
//   Double_t maxval = 0;

//   // for (unsigned counter_to_compare = 0; counter_to_compare != thresholds_to_compare.size(); ++counter_to_compare) {
//   //   TH1F *histogram_to_compare = new TH1F(Form("h_thr%.1f", thresholds_to_compare[counter_to_compare]), Form("h_thr%.1f", thresholds_to_compare[counter_to_compare]), 100, 0, 0);
//   //   histogram_to_compare = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", 3.0, thresholds_to_compare[counter_to_compare]));
//   //   if (counter_to_compare == 0) {
//   //     maxval = histogram_to_compare->GetMaximum();
//   //     base_histogram->SetMaximum(1.1*maxval);
//   //     base_histogram->Draw();
//   //   }
//   //   histogram_to_compare->SetMaximum(1.1*maxval);
//   //   histogram_to_compare->Draw("same");
//   //   delete histogram_to_compare;
//   // }

//   // for (unsigned counter_to_compare = 0; counter_to_compare != thresholds_to_compare.size(); ++counter_to_compare) {
//   unsigned counter_to_compare = 0;
//   TH1F *histogram_to_compare = new TH1F(Form("h_thr%.1f", thresholds_to_compare[counter_to_compare]), Form("h_thr%.1f", thresholds_to_compare[counter_to_compare]), 100, 0, 0);
//   histogram_to_compare = (TH1F*)f_histograms_source->Get(Form("h_fractional_difference_distribution_et%.0f_eta2.100_thr%.1f", 3.0, thresholds_to_compare[counter_to_compare]));
//   // if (counter_to_compare == 0) {
//   //   // maxval = histogram_to_compare->GetMaximum();
//   //   // base_histogram->SetMaximum(1.1*maxval);
//   // }
//   // histogram_to_compare->SetMaximum(1.1*maxval);
//   // histogram_to_compare->Draw("same");
//   // delete histogram_to_compare;
//   // }
  
//   // generate_overlays_given_et(f_histograms_source, first_four, 1, base_threshold, thresholds_to_compare, 3.0);
//   // generate_overlays_given_et(f_histograms_source, first_four, 2, base_threshold, thresholds_to_compare, 5.0);
//   // generate_overlays_given_et(f_histograms_source, first_four, 3, base_threshold, thresholds_to_compare, 10.0);
//   // generate_overlays_given_et(f_histograms_source, first_four, 4, base_threshold, thresholds_to_compare, 30.0);
//   first_four->Update();
//   first_four->Print("et_3_5_10_30.png");
//   delete first_four;

//   delete f_histograms_source;

//   // Int_t canvas_index = 0;
//   // TCanvas *c_xy_profile;
//   // for (std::vector<unsigned int>::iterator layer = layers.begin(); layer != layers.end(); ++layer) {
//   //   for (std::vector<TString>::iterator cut_name = cut_names.begin(); cut_name != cut_names.end(); ++cut_name) {
//   //     std::cout << "Threshold 0.5" << std::endl;
//   //     std::cout << "canvas_index : " << canvas_index << "; reading canvas for layer = " << *layer << ", cut_name = " << *cut_name << std::endl;
//   //     c_xy_profile = (TCanvas*)files[0]->Get(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer));
//   //     // xy_profiles.push_back((TCanvas*)files[0]->Get(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer)));
//   //     TH1F *h_xy_profile = (TH1F*)c_xy_profile->GetPrimitive(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer));
//   //     // delete c_xy_profile;
//   //     TCanvas *myc = new TCanvas("myc", "test");
//   //     myc->cd();
//   //     h_xy_profile->Draw();
//   //     myc->Update();
//   //     myc->Print(Form("x_y_positions_") + *cut_name + Form("_layer_%i_thr0.5.png",*layer));
//   //     // delete myc;
//   //     delete h_xy_profile;
//   //     ++canvas_index;

//   //     std::cout << "Threshold 5.0" << std::endl;
//   //     std::cout << "canvas_index : " << canvas_index << "; reading canvas for layer = " << *layer << ", cut_name = " << *cut_name << std::endl;
//   //     c_xy_profile = (TCanvas*)files[4]->Get(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer));
//   //     // xy_profiles.push_back((TCanvas*)files[0]->Get(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer)));
//   //     h_xy_profile = (TH1F*)c_xy_profile->GetPrimitive(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer));
//   //     // delete c_xy_profile;
//   //     myc = new TCanvas("myc", "test");
//   //     myc->cd();
//   //     h_xy_profile->Draw();
//   //     myc->Update();
//   //     myc->Print(Form("x_y_positions_") + *cut_name + Form("_layer_%i_thr5.0.png",*layer));
//   //     // delete myc;
//   //     delete h_xy_profile;
//   //     ++canvas_index;
//   //   }
//   // }

//   // canvas_index = 0;

//   // TCanvas *c_shower_profile;
//   // // for (std::vector<unsigned int>::iterator layer = layers.begin(); layer != layers.end(); ++layer) {
//   // for (std::vector<TString>::iterator cut_name = cut_names.begin(); cut_name != cut_names.end(); ++cut_name) {
//   //   std::cout << "Threshold 0.5" << std::endl;
//   //   std::cout << "canvas_index : " << canvas_index << "; reading canvas for cut_name = " << *cut_name << std::endl;
//   //   c_shower_profile = (TCanvas*)files[0]->Get(Form("shower_profile_")+*cut_name+Form("et100"));
//   //   // xy_profiles.push_back((TCanvas*)files[0]->Get(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer)));
//   //   TH1F *h_shower_profile = (TH1F*)c_shower_profile->GetPrimitive(Form("Shower Profile ")+*cut_name);
//   //   // delete c_xy_profile;
//   //   TCanvas *myc = new TCanvas("myc", "test");
//   //   myc->cd();
//   //   h_shower_profile->Draw();
//   //   myc->Update();
//   //   myc->Print(Form("shower_profile_") + *cut_name + Form("_thr0.5.png"));
//   //   // delete myc;
//   //   delete h_shower_profile;
//   //   // delete myc;
//   //   ++canvas_index;

//   //   std::cout << "Threshold 5.0" << std::endl;
//   //   std::cout << "canvas_index : " << canvas_index << "; reading canvas for cut_name = " << *cut_name << std::endl;
//   //   c_shower_profile = (TCanvas*)files[4]->Get(Form("shower_profile_")+*cut_name+Form("et100"));
//   //   // xy_profiles.push_back((TCanvas*)files[0]->Get(Form("xy_positions_") + *cut_name + Form("_et100_layer_%i",*layer)));
//   //   h_shower_profile = (TH1F*)c_shower_profile->GetPrimitive(Form("Shower Profile ")+*cut_name);
//   //   // delete c_xy_profile;
//   //   myc = new TCanvas("myc", "test");
//   //   myc->cd();
//   //   h_shower_profile->Draw();
//   //   myc->Update();
//   //   myc->Print(Form("shower_profile_") + *cut_name + Form("_thr5.0.png"));
//   //   // delete myc;
//   //   delete h_shower_profile;
//   //   ++canvas_index;
//   // }
//   // // }
  

//   // // Double_t h_2_mips_max = h_2_mips_cut->GetMaximum();
//   // // Double_t h_5_mips_max = h_5_mips_cut->GetMaximum();

//   // // if (h_2_mips_max > h_5_mips_max) {
//   // //   h_2_mips_cut->SetMaximum(h_2_mips_max);
//   // //   h_5_mips_cut->SetMaximum(h_2_mips_max);
//   // // }
//   // // else {
//   // //   h_2_mips_cut->SetMaximum(h_5_mips_max);
//   // //   h_5_mips_cut->SetMaximum(h_5_mips_max);
//   // // }

  

//   // TCanvas *c_2_mips_cut = (TCanvas*)files[2]->Get("total_energy_distributionet100");
//   // TH1F *h_2_mips_cut;
//   // h_2_mips_cut = (TH1F*)c_2_mips_cut->GetPrimitive("Energy Distribution");
//   // TCanvas *c_5_mips_cut = (TCanvas*)files[4]->Get("total_energy_distributionet100");
//   // TH1F *h_5_mips_cut;
//   // h_5_mips_cut = (TH1F*)c_5_mips_cut->GetPrimitive("Energy Distribution");

//   // Double_t maxval = (h_2_mips_cut->GetMaximum() > h_5_mips_cut->GetMaximum()) ? 1.1 * h_2_mips_cut->GetMaximum() : 1.1 * h_5_mips_cut->GetMaximum();

//   // h_2_mips_cut->SetMaximum(maxval);
//   // h_5_mips_cut->SetMaximum(maxval);
  
//   // TCanvas *myc = new TCanvas("myc", "test");
//   // // myc = new TCanvas("myc", "test");
//   // myc->cd();
//   // h_2_mips_cut->GetXaxis()->SetRangeUser(19000, 23000);
//   // h_2_mips_cut->Draw();
//   // myc->Update();
//   // h_5_mips_cut->Draw("same");
//   // myc->Update();
//   // myc->Print("overlay_2_vs_5mipsthreshold.png");
//   // // output_file->WriteTObject(myc);
//   // delete myc;
//   // // delete file;
//   // // delete output_file;
// }
