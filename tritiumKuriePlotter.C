#include "TROOT.h"
#include "TFile.h"
#include "Riostream.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TH2F.h"
#include "TVirtualPad.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TVector2.h"
#include "TProfile.h"

using std::cout;
using std::endl;

void   divideByNumberEntries(TH1F* histo);
void   divideByBinWidth(TH1F* histo);
double integral(TH1F* histo, double min, double max);
double average(TH1F* histo);
void   resetToFunctionRatio(TH1F* histo, TF1* func);
TH1F*  efficiency(TH1F* numerator, TH1F* denominator);
void   prepareCanvas(TCanvas*& c1,TPad*& p1, int nDivisionsX, int nDivisionsY);
void   prettify(TH1F* h);
void   prettify(TH2F* h);
void   prettify(TF1* h);
void   prettify(TGraph* h);
void   prettify(TProfile* h);
TH2F*  getRangeHisto(double xmin, double xmax, double ymin, double ymax, char* title);



// MAIN ANALYSIS FUNCTION :
void tritiumKuriePlotter() {

  // General macro to analyse testbeam results
  //
  // To run this macro :
  // root [] gSystem->CompileMacro("tritiumKuriePlotter.C");
  // root [] tritiumKuriePlotter()

  // Aesthetics :
  gStyle->SetOptStat(00000000);
  // gStyle->SetFillColor(0);
  gROOT->SetStyle("Plain");
  
  TCanvas* c1;
  TPad*    p1;
  
  
  /*
  // ---------------------------------------------------------------------------------------------- //
  // These plots are for minimum neutrino mass and 1, 10 meV resolution :
  //
  // Read in the spectra (this is all hard-wired and not smart) :
  TFile* root_file_1meV           = new TFile("kurie_1meV_resolution.root","READ");
  TF1*   kurie_smeared_1meV       = (TF1*)(root_file_1meV->Get("T Kurie (smeared)"));
  TF1*   kurie_smeared_1meV_null  = (TF1*)(root_file_1meV->Get("T Kurie (smeared m0)"));
  TH1F*  kurie_smeared_1meV_data  = (TH1F*)(root_file_1meV->Get("kurie_data"));
  TFile* root_file_10meV          = new TFile("kurie_10meV_resolution.root","READ");
  TF1*   kurie_smeared_10meV      = (TF1*)(root_file_10meV->Get("T Kurie (smeared)"));
  TF1*   kurie_smeared_10meV_null = (TF1*)(root_file_10meV->Get("T Kurie (smeared m0)"));
  TH1F*  kurie_smeared_10meV_data = (TH1F*)(root_file_10meV->Get("kurie_data"));
  //
  TCanvas* ch1 = new TCanvas("H-3","H-3",600,600);
  ch1->SetGridx();
  ch1->SetGridy();
  kurie_smeared_1meV->SetTitle(";Energy [keV]; K(T)");
  kurie_smeared_1meV->GetXaxis()->SetLabelSize(0.03);
  kurie_smeared_1meV->GetYaxis()->SetLabelSize(0.03);
  kurie_smeared_1meV->Draw();
  kurie_smeared_1meV_null->SetLineStyle(2);
  kurie_smeared_1meV_null->Draw("same");
  kurie_smeared_1meV_data->Draw("same");
  //
  kurie_smeared_10meV->SetLineColor(4);
  kurie_smeared_10meV->Draw("same");
  kurie_smeared_10meV_null->SetLineColor(4);
  kurie_smeared_10meV_null->SetLineStyle(2);
  kurie_smeared_10meV_null->Draw("same");
  // kurie_smeared_10meV_data->Draw("same");
  //
  // ---------------------------------------------------------------------------------------------- //
  */
  
  // ---------------------------------------------------------------------------------------------- //
  // These plots are for 100 meV neutrino mass and 18.6 meV, 100 meV resolution :
  //
  // kurie_100meV_100meV_resolution.root
  // kurie_100meV_18.6meV_resolution.root
  // kurie_100meV_neg_resolution.root
  //
  // Read in the spectra (this is all hard-wired and not smart) :
  TFile* root_file_true                = new TFile("kurie_100meV_neg_resolution.root","READ");
  TF1*   kurie_100meV_unsmeared        = (TF1*)(root_file_true->Get("T Kurie (smeared)"));
  TF1*   kurie_0mev_unsmeared          = (TF1*)(root_file_true->Get("T Kurie (smeared m0)"));
  TH1F*  kurie_100meV_unsmeared_data   = (TH1F*)(root_file_true->Get("kurie_data"));
  TFile* root_file_100meV_18_6meV_res  = new TFile("kurie_100meV_18.6meV_resolution.root","READ");
  TF1*   kurie_100meV_18_6meV_res      = (TF1*)(root_file_100meV_18_6meV_res->Get("T Kurie (smeared)"));
  TH1F*  kurie_100meV_18_6meV_res_data = (TH1F*)(root_file_100meV_18_6meV_res->Get("kurie_data"));
  TFile* root_file_100meV_100meV_res   = new TFile("kurie_100meV_100meV_resolution.root","READ");
  TF1*   kurie_100meV_100meV_res       = (TF1*)(root_file_100meV_100meV_res->Get("T Kurie (smeared)"));
  TH1F*  kurie_100meV_100meV_res_data  = (TH1F*)(root_file_100meV_100meV_res->Get("kurie_data"));
  //
  TCanvas* ch1 = new TCanvas("H-3","H-3",600,600);
  ch1->SetGridx();
  ch1->SetGridy();
  kurie_0mev_unsmeared->SetTitle(";Energy [keV]; K(T)");
  kurie_0mev_unsmeared->GetXaxis()->SetLabelSize(0.03);
  kurie_0mev_unsmeared->GetYaxis()->SetLabelSize(0.03);
  kurie_0mev_unsmeared->SetLineColor(1);
  kurie_0mev_unsmeared->SetLineStyle(2);
  kurie_0mev_unsmeared->Draw();
  
  kurie_100meV_unsmeared->SetLineColor(1);
  kurie_100meV_unsmeared->SetLineStyle(0);
  kurie_100meV_unsmeared->Draw("same");
  kurie_100meV_unsmeared_data->SetLineColor(1);
  kurie_100meV_unsmeared_data->SetMarkerStyle(kOpenSquare);
  kurie_100meV_unsmeared_data->SetMarkerSize(0.5);
  kurie_100meV_unsmeared_data->Draw("same");
  
  kurie_100meV_18_6meV_res->SetLineColor(2);
  kurie_100meV_18_6meV_res->SetLineStyle(0);
  kurie_100meV_18_6meV_res->Draw("same");
  
  kurie_100meV_100meV_res->SetLineColor(4);
  kurie_100meV_100meV_res->SetLineStyle(0);
  kurie_100meV_100meV_res->Draw("same");
  //
  // ---------------------------------------------------------------------------------------------- //
  
  
  
}

void divideByNumberEntries(TH1F* histo) {

  // Divide each bin content by the number of entries in the histogram
  // Useful for converting frequency distributions into probability distributions
  double integralByHand = 0;
  for (int ibin = 1; ibin <= histo->GetNbinsX(); ++ibin) {
    integralByHand += histo->GetBinContent(ibin);
    // No protection against zero-width bins :
    histo->SetBinContent(ibin,histo->GetBinContent(ibin)/histo->GetEntries());
  }
  cout << "integral =  " << integralByHand << endl;

}

void divideByBinWidth(TH1F* histo) {

  // Divide each bin content by the bin-width.
  // Useful for converting frequency distributions into differential cross-sections.
  for (int ibin = 1; ibin <= histo->GetNbinsX(); ++ibin) {
    // No protection against zero-width bins :
    histo->SetBinContent(ibin,histo->GetBinContent(ibin)/histo->GetBinWidth(ibin));
  }
  
}

double integral(TH1F* histo, double min, double max) {
  
  double integralValue = 0.0;
  for (int ibin = 1; ibin <= histo->GetNbinsX(); ++ibin) {
    double binCenter = histo->GetBinCenter(ibin);
    if (binCenter >= min &&
	binCenter <= max) integralValue += histo->GetBinContent(ibin)*histo->GetBinWidth(ibin);
  }
 
  return integralValue;
 
}

double average(TH1F* histo) {
  
  double sumxfx = 0.0;
  double sumfx  = 0.0;
  double mean   = 0.0;
  for (int ibin = 1; ibin <= histo->GetNbinsX(); ++ibin) {
    sumxfx += histo->GetBinContent(ibin)*histo->GetBinCenter(ibin);
    sumfx  += histo->GetBinContent(ibin);
  }
  if (sumfx > 0.0) mean = sumxfx/sumfx;
 
  return mean;
 
}

void resetToFunctionRatio(TH1F* histo, TF1* func) {

  // Reset histogram contents to ratio with provided function :
  for (int ibin = 1; ibin <= histo->GetNbinsX(); ++ibin) {
    // Integrate the function over the bin width :
    double rangeLo = histo->GetBinLowEdge(ibin);
    double rangeHi = histo->GetBinLowEdge(ibin)+histo->GetBinWidth(ibin);
    double denominator = func->Integral(rangeLo,rangeHi)/histo->GetBinWidth(ibin);
    // cout << "rangeLo, rangeHi, width = " << rangeLo << " , " << rangeHi << " , " << histo->GetBinWidth(ibin) << endl;
    // histo->SetBinContent(ibin,histo->GetBinContent(ibin)/func->Eval(histo->GetBinCenter(ibin)));
    histo->SetBinContent(ibin,histo->GetBinContent(ibin)/denominator);
  }  
  
}

TH1F* efficiency(TH1F* numerator, TH1F* denominator) {

  // Calculate efficiency assuming numerator is a subset of denominator and integer numbers of events.
  TH1F* efficiencyResult = new TH1F(*numerator);
  efficiencyResult->SetTitle("efficiency");
  for (int ibin=0; ibin <= denominator->GetNbinsX()+1; ++ibin) {
    double num = numerator->GetBinContent(ibin);
    double den = denominator->GetBinContent(ibin);
    double eff = 0.0;
    double effErr = 0.0;
    if (den > 0.0) {
      eff = num/den;
      effErr = sqrt(eff*(1.0-eff)/den);
    }
    efficiencyResult->SetBinContent(ibin,eff);
    efficiencyResult->SetBinError(ibin,effErr);
  }
  return efficiencyResult;
}



void prepareCanvas(TCanvas*& c1,TPad*& p1, int nDivisionsX, int nDivisionsY) {
  int nPixelX = 400;
  int nPixelY = 400;
  if (nDivisionsY == 1) nPixelY = 400;
  c1 = new TCanvas("c1","plot_output",nPixelX,nPixelY);
  c1->GetFrame()->SetFillColor(0);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  c1->GetFrame()->SetLineColor(2);
  p1 = new TPad("plot_output","plot_output",0.00,0.00,1.00,1.00);
  p1->Divide(nDivisionsX,nDivisionsY,0.01,0.01,0);
  p1->Draw(); p1->Range(0,0,1,1);
}

void prettify(TH2F* h) {
  
  h->SetTitle("");

  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.8);  
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetTitleSize(0.05);
  h->SetLineWidth(1);

}

void prettify(TH1F* h) {
  
  h->SetTitle("");
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetTitleOffset(0.9);  
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleSize(0.05);
  h->SetLineWidth(1);

}

void prettify(TF1* f) {
  
  f->SetTitle("");
  f->GetXaxis()->CenterTitle();
  f->GetXaxis()->SetTitleOffset(1.2);
  f->GetYaxis()->CenterTitle();
  f->SetLineWidth(1);

}


void prettify(TGraph* g) {
  
  g->SetTitle("");
  g->GetXaxis()->CenterTitle();
  g->GetXaxis()->SetTitleOffset(1.2);
  g->GetYaxis()->CenterTitle();
  g->SetLineWidth(1);

}

void prettify(TProfile* h) {
  
  h->SetTitle("");
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(1.2);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleOffset(0.7);  
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleSize(0.05);
  h->SetLineWidth(1);

}


TH2F*  getRangeHisto(double xmin, double xmax, double ymin, double ymax, char* title) {

  TH2F* plotRange = new TH2F(title,title,100,xmin,xmax,100,ymin,ymax);

  prettify(plotRange);
  plotRange->SetTitle(title);

  return plotRange;
}

