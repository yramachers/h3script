#include <iostream>
// ******
// Usage: start ROOT
// .L h3script.C
// plotcurve(0.0)
// ** for a zero lightest neutrino mass eigenstate
// ** or
// plotcurve(1.0e-2))
// ** for a 10 meV lightest neutrino mass eigenstate
// ******
// simple Kurie plot script for Tritium
// includes energy resolution smearing
// and neutrino mixing effects with
// mixings and mass splittings hard-coded, see
// betashape function below.
// ******************************************************
// *** run plotcurve(neutrino mass in eV) function at the end
// *** edit hard-coded energy resolution value in plotcurve(mn)
// *** Output: - integral spectrum fraction of range of interest
// *** and - plot of derivative of Kurie curve to see
// *** shape changes clearly.
// *** INFO: Tritium has the specific activity of 357 TBq/g,
// *** of interest for integral fraction scaling.
// *******************************************************

// ------------------------------------------------------- //
// Tritium Data :
Double_t Qvalue = 18.574; // units are [keV]
Double_t T_half = 12.32;  // units are [yr]
// ------------------------------------------------------- //


TCanvas* ch3;
TF1* fh3;
TF1* fh3n;
TF1* fh3n_smear;
TF1* kernel;
TF1* kh3_smear;
TF1* kh3;
TF1* fh3_mn0;
TF1* fh3n_mn0;
TF1* fh3n_smear_mn0;
TF1* kernel_mn0;
TF1* kh3_smear_mn0;
TF1* kh3_mn0;
Double_t neutrinoMass[3]; // ** in [eV] **
Double_t PMNS_Ue[3];
Double_t effectiveNeutrinoMass(Double_t mneV);

Double_t fermiFunction(Double_t T);

// from Schmitz p.222-223, eqn.6.69
// citing J.L. Vuilleumier, Rep. Progr. Phys. 49 (1986) 1293
Double_t betashape(Double_t *x, Double_t *par)
{
  // xx is in energy units [keV] and represents kinetic energy!
  Double_t xx = x[0];
  Double_t Phasespacef;
  Double_t result;
  Double_t sum;
  Double_t eps;
  
  // constants
  Double_t electronmass = 511.0; // [keV]
  Double_t pi = TMath::Pi();
  Double_t finestructure = 1.0/137.0;
  Double_t constant;

  Double_t mneutrino_keV = std::max(0.0,par[3])/1000.; // in [keV]

  if (xx >= (Qvalue-mneutrino_keV)) return 0.0;
  if (xx <= 1.0e-6) xx = 1.0e-6;

  Double_t momentum = TMath::Sqrt(xx * (xx + 2.0*electronmass));
  Double_t Fermif = fermiFunction(xx);

  effectiveNeutrinoMass(par[3]);
  
  sum = 0.0;
  eps = (Qvalue - xx)*(Qvalue - xx); // [keV]
  for (Int_t i=0;i<3;i++) 
  {
    if (eps >= (neutrinoMass[i]*neutrinoMass[i]*1.0e-6))
      sum += PMNS_Ue[i]*PMNS_Ue[i] * TMath::Sqrt(eps - neutrinoMass[i]*neutrinoMass[i]*1.0e-6);
  }
  Phasespacef = (xx + electronmass)*momentum * (Qvalue - xx) * sum;

  result = Phasespacef * Fermif;
  return result;
}

Double_t effectiveNeutrinoMass(Double_t mneV)
{
  // The units here are ** [eV] **

  // cout << "mneV = " << mneV << endl;
  
  // nu levels; pick mean values of ranges
  // measured weights: Ue1=0.8-0.84, Ue2=0.52-0.58, Ue3=0.14-0.16
  // levels: mn=me, dm12^2=7.5e-5 eV^2, dm23^2=2.4e-3 eV^2
  // m2_min = 9 meV, m3 = sqrt(m2^2+2.4e-3) = 50 meV
  if (mneV>=0.0) {
    // Normal hierarchy :
    neutrinoMass[0] = mneV;
    neutrinoMass[1] = TMath::Sqrt(mneV * mneV + 7.5e-5);
    neutrinoMass[2] = TMath::Sqrt(neutrinoMass[1]*neutrinoMass[1] + 2.4e-3);
    // Inverted hierarchy :
    // neutrinoMass[2] = mneV;
    // neutrinoMass[0] = TMath::Sqrt(mneV*mneV + 2.4e-3);
    // neutrinoMass[1] = TMath::Sqrt(neutrinoMass[0]*neutrinoMass[0] + 7.5e-5);
  } else {
    neutrinoMass[0] = 0.0;
    neutrinoMass[1] = 0.0;
    neutrinoMass[2] = 0.0;
  }

  // cout << "neutrinoMass[2] = " << neutrinoMass[2] << endl;

  PMNS_Ue[0] = 0.82;
  PMNS_Ue[1] = 0.55;
  PMNS_Ue[2] = 0.15;
  
  Double_t mv_effective2 = 0;
  for (int im=0; im<3; ++im){
    mv_effective2 += pow(neutrinoMass[im],2)*pow(PMNS_Ue[im],2);
  }
  return sqrt(mv_effective2);
}

Double_t normalisedbeta(Double_t *x, Double_t *par)
{
// global Normfactor from fh3->Integral(10^-6,Qvalue)
  
  Double_t result = betashape(x,par) * 1.0/par[2]; 
  return result;
}

Double_t kuriefunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t full = normalisedbeta(x,par);
  Double_t mneutrino_keV = std::max(0.0,par[3])/1000.;
  if (xx >= (Qvalue-mneutrino_keV)) return 0.0;

  Double_t electronmass = 511.0; // [keV]
  Double_t momentum = TMath::Sqrt(xx * (xx + 2.0*electronmass));
  Double_t Fermif = fermiFunction(xx);
  
  Double_t Kurief = TMath::Sqrt(full / ((xx + electronmass)*momentum*Fermif));
  return Kurief;
}

Double_t simplesmooth(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t Eresolution = par[0];
  Double_t weight = TMath::Gaus(xx,par[1],Eresolution,kTRUE);
  //  return kuriefunction(x,par) * weight;
  return normalisedbeta(x,par) * weight;
}

Double_t kurie_smeared(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t left = xx - 5.0*par[0];
  if (left < 15.0) left = 15.0;
  Double_t right = xx + 5.0*par[0];
  Int_t np = 200;
  Double_t* w = new Double_t [np];
  Double_t* y = new Double_t [np];
  
  // Make the kernel dynamically, for any resolution & neutrino mass :
  Double_t Eres = par[0];
  TF1* kernel = new TF1("corefunc",simplesmooth,18.5,Qvalue+5*Eres,4);
  kernel->SetParameter(0,Eres);
  kernel->SetParameter(2,par[2]); // normalisation factor
  kernel->SetParameter(3,par[3]); // neutrino mass
  
  kernel->SetParameter(1,xx);
  kernel->CalcGaussLegendreSamplingPoints(np,y,w,1.e-8);
  Double_t full = kernel->IntegralFast(np,y,w,left,right);

  Double_t electronmass = 511.0; // [keV]
  Double_t momentum = TMath::Sqrt(xx * (xx + 2.0*electronmass));
  Double_t Fermif = fermiFunction(xx);

  Double_t Kurief = TMath::Sqrt(full / ((xx + electronmass)*momentum*Fermif));
  delete [] w;
  delete [] y;
  delete kernel;
  return Kurief;
}

Double_t normalisedbeta_smeared(Double_t *x, Double_t *par)
{
  
  Double_t xx = x[0];
  Double_t left = xx - 5.0*par[0];
  if (left < 0.0) left = 0.0;
  Double_t right = xx + 5.0*par[0];
  Int_t np = 200;
  Double_t* w = new Double_t [np];
  Double_t* y = new Double_t [np];
  
  // Make the kernel dynamically, for any resolution & neutrino mass :
  Double_t Eres = par[0];
  TF1* kernel = new TF1("corefunc",simplesmooth,18.5,Qvalue+5*Eres,4);
  kernel->SetParameter(0,Eres);
  kernel->SetParameter(2,par[2]); // normalisation factor
  kernel->SetParameter(3,par[3]); // neutrino mass
  
  kernel->SetParameter(1,xx);
  kernel->CalcGaussLegendreSamplingPoints(np,y,w,1.e-8);
  Double_t full = kernel->IntegralFast(np,y,w,left,right);
  
  delete [] w;
  delete [] y;
  delete kernel;
  return full;
}


Double_t fermiFunction(Double_t T)
{
  // constants
  Double_t electronmass = 511.0; // [keV]
  Double_t pi = TMath::Pi();
  Double_t finestructure = 1.0/137.0;
  // For tritium :
  Double_t Z = 1.0;
  
  // for eta, get v in units of c from p/m (Schmitz p.222, this is fully relativistic):
  Double_t momentum = TMath::Sqrt(T * (T + 2.0*electronmass));
  Double_t eta = (T + electronmass)/momentum * (Z+1) *finestructure;
  return 2*pi*eta/(1.0-TMath::Exp(-2.0*pi*eta));

}

//-----------------------------------------------------
// run plotcurve function
// .L h3script.C
// plotcurve(1.0e-2)
// for a 10 meV lightest neutrino mass eigenstate
// uses neutrino mass in units of [eV]
// for script convenience, not user convenience
//-----------------------------------------------------
void plotcurve(Double_t mn) // mn in [eV]
{
  cout.precision(10);
  
  
  // Define some common choices for sensitivity studies :

  
  // Here is a general attempt to set the left-limit based on the mass, but it may be overriden below:
  Double_t range_left;
  if (mn > 0.0)
    range_left = Qvalue - 5*(mn/1000.); // [keV]
  else
    range_left = Qvalue - 6.0e-5; // [keV]

  
  // --------------------------------------------- //
  // For 100 meV (NH) :
  Double_t exposure = 5.0E18;
  // Double_t Eres = 18.6e-6;  // (note : cannot be 0)
  // Double_t Eres = 100.0e-6;  // (note : cannot be 0)
  Double_t Eres = 100.0e-8;  // (note : cannot be 0)
  // --------------------------------------------- //
  // For minimum m_v (NH) :
  // Double_t exposure = 5.0E21;
  // Double_t Eres = 1.0e-6;
  // --------------------------------------------- //
  // For minimum m_v (IH) :
  // Double_t exposure = 7.0E18;
  // Double_t Eres = 1.0e-8
  // --------------------------------------------- //
  
  
  // Double_t Eres = 1.0e-5; // Gauss std [keV]
  // Double_t Eres = 1.0e-8;
  
  // This also sets up the neutrino mass array, and must be called :
  cout << " Effective neutrino mass = " << effectiveNeutrinoMass(mn) << " eV " << endl;
  cout << " Energy resolution       = " << Eres << " keV " << endl;
  
  
  // Hard-wire this. For small neutrino masses a 0.1 eV range works :
  // range_left = Qvalue - 3.0e-4;
  // For larger neutrino masses :
  // range_left = Qvalue - 0.25e-3;
  
  // Double_t range_right = Qvalue + 5.0*Eres;
  Double_t range_right = Qvalue + 10.0E-06;
  std::cout << " range_left, range_right: " << range_left << " , " << range_right << std::endl;

  fh3 = new TF1("T beta spectrum",betashape,1.0e-6,18.6,4);
  fh3->SetParameter(0,0.0);
  fh3->SetParameter(1,0.0);
  fh3->SetParameter(2,0.0);
  fh3->SetParameter(3,mn);
  Double_t norm = fh3->Integral(1.0e-6,(Qvalue-mn));

  // UNSMEARED BETA-SPECTRUM :
  fh3n = new TF1("T beta spectrum",normalisedbeta,range_left,range_right,4);
  fh3n->SetParameter(0,Eres);
  fh3n->SetParameter(1,0.0);
  fh3n->SetParameter(2,norm);
  fh3n->SetParameter(3,mn);

  // SMEARED BETA-SPECTRUM :
  fh3n_smear = new TF1("T beta spectrum (smeared)",normalisedbeta_smeared,range_left,range_right,4);
  fh3n_smear->SetParameter(0,Eres);
  fh3n_smear->SetParameter(1,0.0);
  fh3n_smear->SetParameter(2,norm);
  fh3n_smear->SetParameter(3,mn);
  
  // SMEARED KURIE FUNCTION :
  kh3_smear = new TF1("T Kurie (smeared)",kurie_smeared,range_left,range_right,4);
  kh3_smear->SetParameter(0,Eres);
  kh3_smear->SetParameter(1,0.0);
  kh3_smear->SetParameter(2,norm);
  kh3_smear->SetParameter(3,mn);
  kh3_smear->SetRange(range_left, range_right);

  fh3->SetNpx(1000);
  fh3n->SetNpx(1000);
  fh3n_smear->SetNpx(1000);
  kh3_smear->SetNpx(1000);

  kh3 = new TF1("T Kurie",kuriefunction,range_left,range_right,4);
  kh3->SetParameter(0,Eres);
  kh3->SetParameter(1,0.0);
  kh3->SetParameter(2,norm);
  kh3->SetParameter(3,mn);
  kh3->SetNpx(10000);

  // Null hypothesis comparison plots:
  // ---------------------------------
  // Use negative mass
  
  fh3_mn0 = new TF1("T beta spectrum",betashape,1.0e-6,18.6,4);
  fh3_mn0->SetParameter(0,0.0);
  fh3_mn0->SetParameter(1,0.0);
  fh3_mn0->SetParameter(2,0.0);
  fh3_mn0->SetParameter(3,-1.0);
  Double_t norm_mn0 = fh3_mn0->Integral(1.0e-6,Qvalue);
  
  // UNSMEARED BETA-SPECTRUM (mv=0) :
  fh3n_mn0 = new TF1("T beta spectrum",normalisedbeta,range_left,range_right,4);
  fh3n_mn0->SetParameter(0,Eres);
  fh3n_mn0->SetParameter(1,0.0);
  fh3n_mn0->SetParameter(2,norm_mn0);
  fh3n_mn0->SetParameter(3,-1.0);
  
  // SMEARED BETA-SPECTRUM (mv=0) :
  fh3n_smear_mn0 = new TF1("T beta spectrum (smeared)",normalisedbeta_smeared,range_left,range_right,4);
  fh3n_smear_mn0->SetParameter(0,Eres);
  fh3n_smear_mn0->SetParameter(1,0.0);
  fh3n_smear_mn0->SetParameter(2,norm_mn0);
  fh3n_smear_mn0->SetParameter(3,-1.0);
  
  // SMEARED KURIE FUNCTION :
  kh3_smear_mn0 = new TF1("T Kurie (smeared m0)",kurie_smeared,range_left,range_right,4);
  kh3_smear_mn0->SetParameter(0,Eres);
  kh3_smear_mn0->SetParameter(1,0.0);
  kh3_smear_mn0->SetParameter(2,norm_mn0);
  kh3_smear_mn0->SetParameter(3,-1.0);
  kh3_smear_mn0->SetRange(range_left, range_right);
  
  fh3_mn0->SetNpx(1000);
  fh3n_mn0->SetNpx(1000);
  fh3n_smear_mn0->SetNpx(1000);
  kh3_smear_mn0->SetNpx(1000);
  
  kh3_mn0 = new TF1("T Kurie",kuriefunction,range_left,range_right,4);
  kh3_mn0->SetParameter(0,Eres);
  kh3_mn0->SetParameter(1,0.0);
  kh3_mn0->SetParameter(2,norm_mn0);
  kh3_mn0->SetParameter(3,-1.0);
  kh3_mn0->SetNpx(10000);

  // Make pseudo-data histograms:
  // ----------------------------
  Double_t lambda_T = TMath::Log(2.0)/T_half;
  Double_t n_decays = exposure*lambda_T;
  
  // NOT RIGHT YET. CHECK WHICH FUNCTION AND STATISTICAL ERRORS.
  
  // Bin by some factor of the resolution :
  // int nbins = int((Qvalue+5*Eres-range_left)/(1.0*Eres));
  // int nbins = int((range_right-range_left)/(2.0*Eres));
  int nbins = 43;
  // 1 meV bins :
  // nbins = int((range_right-range_left)*1.0E06);
  
  cout << " Kurie pseudo-data plot : nbins = " << nbins << endl;
  TH1F* kurie_hist = new TH1F("kurie_data","kurie_data",nbins,range_left,range_right);
  Double_t nEventsBinnedTotal = 0;
  Double_t spectrumFraction = 0;
  Double_t chisq = 0;
  Int_t nbins_chisq = 0;
  for (int ibin = 1; ibin <= nbins; ibin++){
    Double_t T_bin = kurie_hist->GetBinCenter(ibin);
    Double_t T_bin_width = kurie_hist->GetBinWidth(ibin);
    // Take the number of events from the actual spectrum, not the Kurie function :
    Double_t N_events = n_decays*T_bin_width*fh3n_smear->Eval(T_bin);
    spectrumFraction   += T_bin_width*fh3n->Eval(T_bin);
    nEventsBinnedTotal += N_events;
    // cout << " Kurie psuedo-data plot : bin center, width, N_events = " << T_bin << " , " << T_bin_width << " , " << N_events << endl;
    kurie_hist->SetBinContent(ibin,kh3_smear->Eval(T_bin));
    if (N_events>0.01) {
      kurie_hist->SetBinError(ibin,(1.0/sqrt(N_events))*kh3_smear->Eval(T_bin));
    } else {
      kurie_hist->SetBinError(ibin,0.0);
    }
    // Calculate the chi-squared relative to null hypothesis:
    // -- only include bins containing a minimum number of events --
    if (N_events > 1.0) {
      Double_t chisq_bin = pow((kurie_hist->GetBinContent(ibin) - kh3_smear_mn0->Eval(T_bin)),2.)/pow(kurie_hist->GetBinError(ibin),2.);
      std::cout << " T, chisq_bin = " << T_bin << " , " << chisq_bin << std::endl;
      chisq += chisq_bin;
      nbins_chisq++;
    }
  }
  
  // Do the plotting:
  // ----------------
  // Aesthetics :
  gStyle->SetOptStat(00000000);
  // gStyle->SetFillColor(0);
  gROOT->SetStyle("Plain");

  ch3 = new TCanvas("H-3","H-3",1200,600);
  ch3->SetGridx();
  ch3->SetGridy();
  // fh3smooth->SetTitle("Kurie function derivative;Energy [keV];dN(E)/dE");
  // fh3smooth->DrawDerivative("AL");
  // fh3smooth->Draw();
  //
  // fh3n->Draw();
  // fh3n_smear->SetLineStyle(2);
  // fh3n_smear->Draw("same");
  
  // fh3n_mn0->SetLineStyle(2);
  // fh3n_mn0->Draw("same");
  
  // kh3_smear->Draw();
  // kh3_smear_mn0->SetLineStyle(2);
  // kh3_smear_mn0->Draw("same");
  
  kh3_smear->SetTitle(";Energy [keV]; K(T)");
  kh3_smear->Draw();
  kh3_smear_mn0->SetLineStyle(2);
  kh3_smear_mn0->Draw("same");
  kurie_hist->Draw("same");

  std::cout << " exposure: " << exposure << std::endl;
  std::cout << " n_decays: " << n_decays << std::endl;
  std::cout << " integral spectrum fraction: " << fh3n->Integral(range_left, Qvalue) << std::endl;
  std::cout << " my spectrum fraction: " << spectrumFraction << std::endl;
  std::cout << " total binned events: " << nEventsBinnedTotal << std::endl;
  std::cout << " chi2 w.r.t. null hypothesis: " << chisq << std::endl;
  std::cout << " chi2 probability: " << TMath::Prob(chisq,1) << std::endl;
  
  
  // --------------------------------------------------------------------------- //
  // Save output for combined plotting :
  TFile* histFile = new TFile("kurie_100meV_neg_resolution.root","RECREATE");
  kh3_smear->Write();
  kh3_smear_mn0->Write();
  kurie_hist->Write();
  histFile->Close();
  // --------------------------------------------------------------------------- //

  
}





