#include "../effSigmaMacro.C"

#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLine.h"
#include "TFitResult.h"

const TString eta_portion = Form("_eta%.3f", 2.0);
const TString version_name = Form("cracks_study");

const Double_t gapWidth = 6.;

const Double_t gapEvenBegin = 224.;
const Double_t gapEvenEnd = gapEvenBegin + gapWidth;
const Double_t gapOddBegin = 254.;
const Double_t gapOddEnd = gapOddBegin + gapWidth;

const Double_t incomingParticleEnergy = 225.731;

const Int_t nEvents_threshold = 80;
Int_t nBinsBelowThreshold = 0;

Int_t m_nBinsRequired;
Double_t *m_arrayWithZeros;
Double_t *m_generatedParticleX;
Double_t *m_nEvents;
std::vector<Double_t> m_badXPositions;
Int_t m_badBins;
Double_t m_baseEnergy = 0;
Double_t m_baseWidth = 0;

Double_t *m_effectiveSigmas;
Double_t *m_effectiveSigmasErrors;

Double_t *m_ratioEffectiveSigmaToPeakEnergy;
Double_t *m_ratioEffectiveSigmaToPeakEnergyErrors;

Double_t *m_ratioFromFit;
Double_t *m_ratioFromFitErrors;

Double_t *m_fractionalEnergyDifference;
Double_t *m_fractionalEnergyDifferenceErrors;

Double_t *m_ratioEffectiveSigmaToPeakEnergyScaled;
Double_t *m_ratioEffectiveSigmaToPeakEnergyScaledErrors;

Double_t *m_ratioFromFitScaled;
Double_t *m_ratioFromFitScaledErrors;

//The following:
Double_t m_effectiveSigmaMax = -1.;
Double_t m_effectiveSigmaMin = -1.;

Double_t m_ratioEffectiveSigmaToPeakEnergyMax = -1.;
Double_t m_ratioEffectiveSigmaToPeakEnergyMin = -1.;

Double_t m_ratioFromFitMax = -1.;
Double_t m_ratioFromFitMin = -1.;

Double_t m_fractionalEnergyDifferenceMax = -1.;
Double_t m_fractionalEnergyDifferenceMin = -1.;

Double_t m_ratioEffectiveSigmaToPeakEnergyScaledMax = -1.;
Double_t m_ratioEffectiveSigmaToPeakEnergyScaledMin = -1.;

Double_t m_ratioFromFitScaledMax = -1.;
Double_t m_ratioFromFitScaledMin = -1.;
// will be reset

Double_t getRatioError(Double_t numerator, Double_t numeratorError, Double_t denominator, Double_t denominatorError) {
  return (numerator/denominator)*sqrt(pow(numeratorError/numerator,2)+pow(denominatorError/denominator,2));
}

void initializeAxisParameters(TFile *histogramsInputFile) {
  TH1F *h_generatedParticlesX = (TH1F*)histogramsInputFile->Get("h_generatedParticlesX");
  TAxis *h_Xaxis = h_generatedParticlesX->GetXaxis();
  m_nBinsRequired = h_Xaxis->GetNbins();
  m_arrayWithZeros = new Double_t[m_nBinsRequired+2];
  m_generatedParticleX = new Double_t[m_nBinsRequired+2];
  m_nEvents = new Double_t[m_nBinsRequired+2];

  std::cout << "Reading in information about x positions of generated particles:" << std::endl;
    
  for (Int_t binCounter = 0; binCounter <= (1+m_nBinsRequired); ++binCounter) {
    m_arrayWithZeros[binCounter] = 0.;
    m_generatedParticleX[binCounter] = h_Xaxis->GetBinCenter(binCounter);
    m_nEvents[binCounter] = h_generatedParticlesX->GetBinContent(binCounter);
    std::cout << "At binCounter = " << binCounter << ", generatedParticleX = " << m_generatedParticleX[binCounter] << ", nEvents = " << m_nEvents[binCounter] << std::endl;
    if (m_nEvents[binCounter] < nEvents_threshold) {
      m_badXPositions.push_back(m_generatedParticleX[binCounter]);
      ++m_badBins;
    }
  }
}

void printBadBinInfo() {
  std::cout << "There are " << m_badBins << " bad bins, at these locations of x:" << std::endl;
  for (std::vector<Double_t>::iterator m_badXPositionsIterator = m_badXPositions.begin(); m_badXPositionsIterator != m_badXPositions.end(); ++m_badXPositionsIterator) {
    std::cout << *m_badXPositionsIterator << std::endl;
  }
}

void cleanupAxisParameters() {
  delete[] m_generatedParticleX;
  delete[] m_nEvents;
}

void initializeArrayWithCorrectBinning(Double_t*& arrayToInitialize) {
  arrayToInitialize = new Double_t[m_nBinsRequired+2];
  for (Int_t binCounter = 0; binCounter <= (1+m_nBinsRequired); ++binCounter) {
    arrayToInitialize[binCounter] = -1.;
  }
}

void cleanupArray(Double_t *arrayToCleanUp) {
  delete[] arrayToCleanUp;
}

Double_t getEffectiveSigma(TH1F *inputHistogram) {
  return effSigmaMacro(inputHistogram);
}

void calculateBaseParameters(TFile *histogramsInputFile) {
  m_baseEnergy = 0;
  m_baseWidth = 0;
  
  TH1F *h_generatedParticlesX = (TH1F*)histogramsInputFile->Get("h_generatedParticlesX");
  TAxis *h_Xaxis = h_generatedParticlesX->GetXaxis();
  Double_t firstBinX = h_Xaxis->GetBinCenter(1);
  Int_t minBinForBaseEnergy = h_Xaxis->FindBin(firstBinX + 0.1*(gapEvenBegin - firstBinX));
  Int_t maxBinForBaseEnergy = h_Xaxis->FindBin(firstBinX + 0.7*(gapEvenBegin - firstBinX));
  std::cout << "minBinForBaseEnergy = " << minBinForBaseEnergy << "; maxBinForBaseEnergy = " << maxBinForBaseEnergy << std::endl;
    
  Int_t binsToCount = 0;
  for (Int_t binCounter = minBinForBaseEnergy; binCounter <= maxBinForBaseEnergy; ++binCounter) {
    if (m_nEvents[binCounter] < nEvents_threshold) {
      continue;
    }
    TH1F *h_energyHistogram = (TH1F*)histogramsInputFile->Get(Form("h_energyHistogram_binNumber_%i",binCounter));
    TFitResultPtr gaussianFit = h_energyHistogram->Fit("gaus", "IMES");
    Double_t fitMean = gaussianFit->Parameter(1);
    // Double_t fitWidth = gaussianFit->Parameter(2);
    Double_t fitWidth = getEffectiveSigma(h_energyHistogram);
    m_baseEnergy += fitMean;
    m_baseWidth += fitWidth;
    // std::cout << "fitMean = " << fitMean << std::endl;
    ++binsToCount;
  }
  m_baseEnergy = m_baseEnergy/binsToCount;
  m_baseWidth = m_baseWidth/binsToCount;
}

Double_t getElementToFill(TH1F *h_energyHistogram, Int_t analysisType) {
  if (analysisType == 1) {
    Double_t effectiveSigma = getEffectiveSigma(h_energyHistogram);
    return effectiveSigma;
  }
  else if (analysisType == 2) {
    Double_t effectiveSigma = getEffectiveSigma(h_energyHistogram);
    Double_t peakEnergy = h_energyHistogram->GetXaxis()->GetBinCenter(h_energyHistogram->GetMaximumBin());
    return effectiveSigma/peakEnergy;
  }
  return -1.;
}

std::pair<Double_t, Double_t> getElementAndErrorToFill(TH1F *h_energyHistogram, Int_t analysisType, Int_t binCounter) {
  std::pair<Double_t, Double_t> elementAndError;
  if (analysisType == 1) {
    Double_t effectiveSigma = getEffectiveSigma(h_energyHistogram);
    elementAndError.first = effectiveSigma;
    elementAndError.second = effectiveSigma/(sqrt(2*m_nEvents[binCounter]));
  }
  else if (analysisType == 2) {
    Double_t effectiveSigma = getEffectiveSigma(h_energyHistogram);
    Double_t peakEnergy = h_energyHistogram->GetXaxis()->GetBinCenter(h_energyHistogram->GetMaximumBin());
    elementAndError.first = effectiveSigma/peakEnergy;
    elementAndError.second = effectiveSigma/(peakEnergy*sqrt(2*m_nEvents[binCounter]));
  }
  else if (analysisType == 3) {
    TFitResultPtr gaussianFit = h_energyHistogram->Fit("gaus", "IMES");
    Double_t fitMean = gaussianFit->Parameter(1);
    Double_t fitSigma = gaussianFit->Parameter(2);
    Double_t fitMeanError = gaussianFit->ParError(1);
    Double_t fitSigmaError = gaussianFit->ParError(2);
    elementAndError.first = fitSigma/fitMean;
    elementAndError.second = getRatioError(fitSigma, fitSigmaError, fitMean, fitMeanError);
  }
  else if (analysisType == 4) {
    TFitResultPtr gaussianFit = h_energyHistogram->Fit("gaus", "IMES");
    Double_t fitMean = gaussianFit->Parameter(1);
    Double_t fitMeanError = gaussianFit->ParError(1);
    elementAndError.first = (m_baseEnergy-fitMean)/m_baseEnergy;
    elementAndError.second = fitMeanError/m_baseEnergy;
  }
  else if (analysisType == 5) {
    Double_t effectiveSigma = getEffectiveSigma(h_energyHistogram);
    Double_t peakEnergy = h_energyHistogram->GetXaxis()->GetBinCenter(h_energyHistogram->GetMaximumBin());
    elementAndError.first = effectiveSigma*m_baseEnergy/(peakEnergy*peakEnergy);
    elementAndError.second = effectiveSigma*m_baseEnergy/(peakEnergy*peakEnergy*sqrt(2*m_nEvents[binCounter]));
  }
  else if (analysisType == 6) {
    TFitResultPtr gaussianFit = h_energyHistogram->Fit("gaus", "IMES");
    Double_t fitMean = gaussianFit->Parameter(1);
    Double_t fitSigma = gaussianFit->Parameter(2);
    Double_t fitMeanError = gaussianFit->ParError(1);
    Double_t fitSigmaError = gaussianFit->ParError(2);
    elementAndError.first = fitSigma*m_baseEnergy/(fitMean*fitMean);
    elementAndError.second = (m_baseEnergy/fitMean)*getRatioError(fitSigma, fitSigmaError, fitMean, fitMeanError);
  }
  else {
    elementAndError.first = -1.;
    elementAndError.second = -1.;
  }
  return elementAndError;
}

void fillArray(TFile *histogramsInputFile, Double_t *arrayToFill, Double_t &arrayMax, Double_t &arrayMin, Int_t analysisType) {
  for (Int_t binCounter = 0; binCounter <= (1+m_nBinsRequired); ++binCounter) {
    if (m_nEvents[binCounter] < nEvents_threshold) {
      continue;
    }
    
    TH1F *h_energyHistogram = (TH1F*)histogramsInputFile->Get(Form("h_energyHistogram_binNumber_%i",binCounter));
    Double_t elementToFill = getElementToFill(h_energyHistogram, analysisType);
    arrayToFill[binCounter] = elementToFill;
    if (arrayMin == -1 || arrayMax == -1) {
      arrayMax = elementToFill;
      arrayMin = elementToFill;
    }
    else {
      if (elementToFill > arrayMax) arrayMax = elementToFill;
      if (elementToFill < arrayMin) arrayMin = elementToFill;
    }
    // std::cout << "For analysis type " << analysisType << ", at binCounter = " << binCounter << ", nEvents = " << m_nEvents[binCounter] << ", x = " << m_generatedParticleX[binCounter] << ", elementToFill = " << elementToFill << ", arrayMax = " << arrayMax << ", arrayMin = " << arrayMin << std::endl;
  }
}

void fillArrayWithError(TFile *histogramsInputFile, Double_t *arrayToFill, Double_t *arrayErrorToFill, Double_t &arrayMax, Double_t &arrayMin, Int_t analysisType) {
  for (Int_t binCounter = 0; binCounter <= (1+m_nBinsRequired); ++binCounter) {
    if (m_nEvents[binCounter] < nEvents_threshold) {
      continue;
    }
    TH1F *h_energyHistogram = (TH1F*)histogramsInputFile->Get(Form("h_energyHistogram_binNumber_%i",binCounter));
    std::pair<Double_t, Double_t> elementAndErrorToFill = getElementAndErrorToFill(h_energyHistogram, analysisType, binCounter);
        
    arrayToFill[binCounter] = elementAndErrorToFill.first;
    arrayErrorToFill[binCounter] = elementAndErrorToFill.second;
    
    if (arrayMin == -1 || arrayMax == -1) {
      arrayMax = elementAndErrorToFill.first;
      arrayMin = elementAndErrorToFill.first;
    }
    else {
      if (elementAndErrorToFill.first > arrayMax) arrayMax = elementAndErrorToFill.first;
      if (elementAndErrorToFill.first < arrayMin) arrayMin = elementAndErrorToFill.first;
    }
    // std::cout << "For analysis type " << analysisType << ", at binCounter = " << binCounter << ", nEvents = " << m_nEvents[binCounter] << ", x = " << m_generatedParticleX[binCounter] << ", elementToFill = " << elementAndErrorToFill.first << ", elementError = " << elementAndErrorToFill.second << ", arrayMax = " << arrayMax << ", arrayMin = " << arrayMin << std::endl;
  }
}

void plotVerticalGapLines(TCanvas *outputCanvas, Double_t ymax, Double_t ymin) {
  outputCanvas->cd();
  TLine *l_gapEvenBegin = new TLine(gapEvenBegin, ymin - 0.085*(ymax - ymin), gapEvenBegin, ymax + 0.085*(ymax - ymin));
  TLine *l_gapEvenEnd = new TLine(gapEvenEnd, ymin - 0.085*(ymax - ymin), gapEvenEnd, ymax + 0.085*(ymax - ymin));
  TLine *l_gapOddBegin = new TLine(gapOddBegin, ymin - 0.085*(ymax - ymin), gapOddBegin, ymax + 0.085*(ymax - ymin));
  TLine *l_gapOddEnd = new TLine(gapOddEnd, ymin - 0.085*(ymax - ymin), gapOddEnd, ymax + 0.085*(ymax - ymin));

  l_gapEvenBegin->SetLineColor(kRed);
  l_gapEvenEnd->SetLineColor(kRed);
  l_gapOddBegin->SetLineColor(kBlue);
  l_gapOddEnd->SetLineColor(kBlue);
  
  l_gapEvenBegin->Draw("SAME");
  l_gapEvenEnd->Draw("SAME");
  l_gapOddBegin->Draw("SAME");
  l_gapOddEnd->Draw("SAME");
}

void plotGraph(TFile *graphsOutputFile, TString title, Double_t *toPlot, Double_t ymax, Double_t ymin) {
  TCanvas *outputCanvas = new TCanvas(title, title, 1024,768);
  TGraph *graphToPlot = new TGraph(m_nBinsRequired+2, m_generatedParticleX, toPlot);
  graphToPlot->SetTitle(title);
  graphToPlot->GetYaxis()->SetRangeUser(ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin));
  graphToPlot->Draw("A*");
  plotVerticalGapLines(outputCanvas, ymax, ymin);
  graphsOutputFile->WriteTObject(outputCanvas);
  delete graphToPlot;
  delete outputCanvas;
}

void plotGraphWithError(TFile *graphsOutputFile, TString title, Double_t *toPlot, Double_t *errorsToPlot, Double_t ymax, Double_t ymin) {
  TCanvas *outputCanvas = new TCanvas(title, title, 1024,768);
  TGraphErrors *graphToPlot = new TGraphErrors(m_nBinsRequired+2, m_generatedParticleX, toPlot, m_arrayWithZeros, errorsToPlot);
  graphToPlot->SetTitle(title);
  graphToPlot->GetYaxis()->SetRangeUser(ymin - 0.1*(ymax-ymin), ymax + 0.1*(ymax-ymin));
  graphToPlot->Draw("A*");
  plotVerticalGapLines(outputCanvas, ymax, ymin);
  graphsOutputFile->WriteTObject(outputCanvas);
  delete graphToPlot;
  delete outputCanvas;
}

// Double_t getOverallEffectiveResolution(TFile *histogramsInputFile) {
//   TH1F *h_allEventsEnergyHistogram = (TH1F*)histogramsInputFile->Get(Form("allEventsEnergyHistogram"));
//   Double_t effectiveSigma = getEffectiveSigma(h_allEventsEnergyHistogram);
//   std::cout << "Effective sigma in chosen units = " << effectiveSigma << std::endl;
//   return effectiveSigma/m_baseEnergy;
// }

Double_t printEffectiveSigma(TFile *histogramsInputFile) {
  TH1F *h_allEventsEnergyHistogram = (TH1F*)histogramsInputFile->Get(Form("allEventsEnergyHistogram"));
  Double_t effectiveSigma = getEffectiveSigma(h_allEventsEnergyHistogram)*incomingParticleEnergy/m_baseEnergy;
  std::cout << "Overall effective sigma = " << effectiveSigma << " GeV" << std::endl;
  std::cout << "Overall effective resolution = " << effectiveSigma/incomingParticleEnergy << std::endl;

  Double_t baseEffectiveSigma = m_baseWidth*incomingParticleEnergy/m_baseEnergy;
  std::cout << "Base effective sigma = " << baseEffectiveSigma << " GeV" << std::endl;
  std::cout << "Base effective resolution = " << baseEffectiveSigma/incomingParticleEnergy << std::endl;
  // return effectiveSigma/m_baseEnergy;
}

void analyzeHistograms() {
  TFile *histogramsInputFile = new TFile(Form("root_histograms/histograms_")+version_name+eta_portion+Form("_results.root"),"READ");
  TFile *graphsOutputFile = new TFile(Form("root_histograms/width_analysis_")+version_name+eta_portion+Form("_results.root"),"RECREATE");
  initializeAxisParameters(histogramsInputFile);

  initializeArrayWithCorrectBinning(m_effectiveSigmas);
  initializeArrayWithCorrectBinning(m_effectiveSigmasErrors);
  // fillArray(histogramsInputFile, m_effectiveSigmas, m_effectiveSigmaMax, m_effectiveSigmaMin, 1);
  // plotGraph(graphsOutputFile, "Effective sigma versus x-position of generated photon", m_effectiveSigmas, m_effectiveSigmaMax, m_effectiveSigmaMin);
  fillArrayWithError(histogramsInputFile, m_effectiveSigmas, m_effectiveSigmasErrors, m_effectiveSigmaMax, m_effectiveSigmaMin, 1);
  plotGraphWithError(graphsOutputFile, "Effective sigma versus x-position of generated photon", m_effectiveSigmas, m_effectiveSigmasErrors, m_effectiveSigmaMax, m_effectiveSigmaMin);
  cleanupArray(m_effectiveSigmas);
  cleanupArray(m_effectiveSigmasErrors);

  initializeArrayWithCorrectBinning(m_ratioEffectiveSigmaToPeakEnergy);
  initializeArrayWithCorrectBinning(m_ratioEffectiveSigmaToPeakEnergyErrors);
  // fillArray(histogramsInputFile, m_ratioEffectiveSigmaToPeakEnergy, m_ratioEffectiveSigmaToPeakEnergyMax, m_ratioEffectiveSigmaToPeakEnergyMin, 2);
  // plotGraph(graphsOutputFile, "Ratio of effective sigma to peak energy versus x-position of generated photon", m_ratioEffectiveSigmaToPeakEnergy, m_ratioEffectiveSigmaToPeakEnergyMax, m_ratioEffectiveSigmaToPeakEnergyMin);
  fillArrayWithError(histogramsInputFile, m_ratioEffectiveSigmaToPeakEnergy, m_ratioEffectiveSigmaToPeakEnergyErrors, m_ratioEffectiveSigmaToPeakEnergyMax, m_ratioEffectiveSigmaToPeakEnergyMin, 2);
  plotGraphWithError(graphsOutputFile, "Ratio of effective sigma to peak energy versus x-position of generated photon", m_ratioEffectiveSigmaToPeakEnergy, m_ratioEffectiveSigmaToPeakEnergyErrors, m_ratioEffectiveSigmaToPeakEnergyMax, m_ratioEffectiveSigmaToPeakEnergyMin);
  cleanupArray(m_ratioEffectiveSigmaToPeakEnergy);
  cleanupArray(m_ratioEffectiveSigmaToPeakEnergyErrors);

  initializeArrayWithCorrectBinning(m_ratioFromFit);
  initializeArrayWithCorrectBinning(m_ratioFromFitErrors);
  fillArrayWithError(histogramsInputFile, m_ratioFromFit, m_ratioFromFitErrors, m_ratioFromFitMax, m_ratioFromFitMin, 3);
  plotGraphWithError(graphsOutputFile, "Ratio of sigma to mean of Gaussian fit versus x-position of generated photon", m_ratioFromFit, m_ratioFromFitErrors, m_ratioFromFitMax, m_ratioFromFitMin);
  cleanupArray(m_ratioFromFit);
  cleanupArray(m_ratioFromFitErrors);

  calculateBaseParameters(histogramsInputFile);
  
  initializeArrayWithCorrectBinning(m_fractionalEnergyDifference);
  initializeArrayWithCorrectBinning(m_fractionalEnergyDifferenceErrors);
  fillArrayWithError(histogramsInputFile, m_fractionalEnergyDifference, m_fractionalEnergyDifferenceErrors, m_fractionalEnergyDifferenceMax, m_fractionalEnergyDifferenceMin, 4);
  plotGraphWithError(graphsOutputFile, "Fractional energy difference versus x-position of generated photon", m_fractionalEnergyDifference, m_fractionalEnergyDifferenceErrors, m_fractionalEnergyDifferenceMax, m_fractionalEnergyDifferenceMin);
  cleanupArray(m_fractionalEnergyDifference);
  cleanupArray(m_fractionalEnergyDifferenceErrors);

  initializeArrayWithCorrectBinning(m_ratioEffectiveSigmaToPeakEnergyScaled);
  initializeArrayWithCorrectBinning(m_ratioEffectiveSigmaToPeakEnergyScaledErrors);
  // fillArray(histogramsInputFile, m_ratioEffectiveSigmaToPeakEnergyScaled, m_ratioEffectiveSigmaToPeakEnergyScaledMax, m_ratioEffectiveSigmaToPeakEnergyScaledMin, 2);
  // plotGraph(graphsOutputFile, "Ratio of effective sigma to peak energy versus x-position of generated photon", m_ratioEffectiveSigmaToPeakEnergyScaled, m_ratioEffectiveSigmaToPeakEnergyScaledMax, m_ratioEffectiveSigmaToPeakEnergyScaledMin);
  fillArrayWithError(histogramsInputFile, m_ratioEffectiveSigmaToPeakEnergyScaled, m_ratioEffectiveSigmaToPeakEnergyScaledErrors, m_ratioEffectiveSigmaToPeakEnergyScaledMax, m_ratioEffectiveSigmaToPeakEnergyScaledMin, 5);
  plotGraphWithError(graphsOutputFile, "Ratio of scaled effective sigma to peak energy versus x-position of generated photon", m_ratioEffectiveSigmaToPeakEnergyScaled, m_ratioEffectiveSigmaToPeakEnergyScaledErrors, m_ratioEffectiveSigmaToPeakEnergyScaledMax, m_ratioEffectiveSigmaToPeakEnergyScaledMin);
  cleanupArray(m_ratioEffectiveSigmaToPeakEnergyScaled);
  cleanupArray(m_ratioEffectiveSigmaToPeakEnergyScaledErrors);

  initializeArrayWithCorrectBinning(m_ratioFromFitScaled);
  initializeArrayWithCorrectBinning(m_ratioFromFitScaledErrors);
  fillArrayWithError(histogramsInputFile, m_ratioFromFitScaled, m_ratioFromFitScaledErrors, m_ratioFromFitScaledMax, m_ratioFromFitScaledMin, 6);
  plotGraphWithError(graphsOutputFile, "Ratio of scaled sigma to mean of Gaussian fit versus x-position of generated photon", m_ratioFromFitScaled, m_ratioFromFitScaledErrors, m_ratioFromFitScaledMax, m_ratioFromFitScaledMin);
  cleanupArray(m_ratioFromFitScaled);
  cleanupArray(m_ratioFromFitScaledErrors);

  // Double_t overallEffectiveResolution = getOverallEffectiveResolution(histogramsInputFile);

  printEffectiveSigma(histogramsInputFile);

  // std::cout << "Overall effective resolution: " << overallEffectiveResolution << std::endl;
  // // std::cout << "... as compared to base resolution: " << m_baseWidth/m_baseEnergy << std::endl;
  // std::cout << "... as compared to base width: " << m_baseWidth << std::endl;
  
  cleanupAxisParameters();
  delete histogramsInputFile;
  delete graphsOutputFile;
  // std::cout << "Done analyzing! " << nBinsBelowThreshold << " bins below threshold." << std::endl;
  printBadBinInfo();
}
