#include <iostream>
// ******
// Usage: start ROOT
// .L h3script.C
// plotcurve(0.0)
// ** for a zero lightest neutrino mass eigenstate
// ** or
// plotcurve(1.0e-5))
// ** for a 10 meV lightest neutrino mass eigenstate
// ******
// simple Kurie plot script for Tritium
// includes energy resolution smearing
// and neutrino mixing effects with
// mixings and mass splittings hard-coded, see
// betashape function below.
// ******************************************************
// *** run plotcurve(neutrino mass in keV) function at the end
// *** edit hard-coded energy resolution value in plotcurve(mn)
// *** Output: - integral spectrum fraction of range of interest
// *** and - plot of derivative of Kurie curve to see
// *** shape changes clearly.
// *** INFO: Tritium has the specific activity of 357 TBq/g,
// *** of interest for integral fraction scaling.
// *******************************************************

Double_t Qvalue = 18.574; // units in [keV]
TCanvas* ch3;
TF1* fh3;
TF1* fh3n;
TF1* fkurie;
TF1* kernel;
TF1* fh3smooth;
TF1* kurieh3;

// from Schmitz p.222-223, eqn.6.69
// citing J.L. Vuilleumier, Rep. Progr. Phys. 49 (1986) 1293
Double_t betashape(Double_t *x, Double_t *par)
{
  // xx is in energy units [keV] and represents kinetic energy!
  Double_t xx = x[0];
  Double_t Fermif;
  Double_t eta,Z;
  Double_t Phasespacef;  
  Double_t momentum;
  Double_t result;
  Double_t sum;
  Double_t eps;
  
  // constants
  Double_t electronmass = 511.0; // [keV]
  Double_t pi = TMath::Pi();
  Double_t finestructure = 1.0/137.0;
  Double_t constant;
  Double_t mneutrino = par[3]; // in [keV]
  Double_t mneV = mneutrino * 1.0e3; // in [eV]
  // H specific, daughter Z, Helium
  Z = 2.0;

  // nu levels; pick mean values of ranges
  // measured weights: Ue1=0.8-0.84, Ue2=0.52-0.58, Ue3=0.14-0.16
  // levels: mn=me, dm12^2=7.5e-5 eV^2, dm23^2=2.4e-3 eV^2
  // m2_min = 9 meV, m3 = sqrt(m2^2+2.4e-3) = 50 meV
  Double_t level2 = TMath::Sqrt(mneV * mneV + 7.5e-5);
  Double_t level3 = TMath::Sqrt(level2*level2 + 2.4e-3);
  Double_t level[3] = {mneV, level2, level3}; // eV
  Double_t weights[3] = {0.82, 0.55, 0.15}; // U_ei measured mixing

  if (xx >= (Qvalue-mneutrino)) return 0.0;
  if (xx <= 1.0e-6) xx = 1.0e-6;

  // for eta, get v in units of c from p/m (non-rel.)
  momentum = TMath::Sqrt(xx * (xx + 2.0*electronmass));
  eta = (xx + electronmass)/momentum * Z*finestructure;
  Fermif = 2*pi*eta/(1.0-TMath::Exp(-2.0*pi*eta));

  sum = 0.0;
  eps = (Qvalue - xx)*(Qvalue - xx); // [keV]
  for (Int_t i=0;i<3;i++) 
  {
    if (eps >= (level[i]*level[i]*1.0e-6)) 
      sum += weights[i]*weights[i] * TMath::Sqrt(eps - level[i]*level[i]*1.0e-6);
  }
  Phasespacef = (xx + electronmass)*momentum * (Qvalue - xx) * sum;

  result = Phasespacef * Fermif;
  return result;
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
  Double_t Kurief,Fermif;
  Double_t eta;
  Double_t momentum;
  Double_t mneutrino = par[3];
  if (xx >= (Qvalue-mneutrino)) return 0.0;

  // constants
  Double_t electronmass = 511.0;
  Double_t pi = TMath::Pi();
  Double_t finestructure = 1.0/137.0;
  // H-3 specific
  Double_t Z = 2.0;

  momentum = TMath::Sqrt(xx * (xx + 2.0*electronmass));
  eta = electronmass/momentum * Z*finestructure;
  Fermif = 2*pi*eta/(1.0-TMath::Exp(-2.0*pi*eta));

  Kurief = TMath::Sqrt(full / ((xx + electronmass)*momentum*Fermif));
  return Kurief;
}

Double_t smoothed(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t left = xx - 5.0*par[0];
  if (left < 15.0) left = 15.0;
  Double_t right = xx + 5.0*par[0];
  Int_t np = 200;
  Double_t* w = new Double_t [np];
  Double_t* y = new Double_t [np];
  kernel->SetParameter(1,xx);
  kernel->CalcGaussLegendreSamplingPoints(np,y,w,1.e-8);
  Double_t full = kernel->IntegralFast(np,y,w,left,right);
  Double_t Kurief,Fermif;
  Double_t eta;
  Double_t momentum;
  Double_t mneutrino = par[3];

  // constants
  Double_t electronmass = 511.0;
  Double_t pi = TMath::Pi();
  Double_t finestructure = 1.0/137.0;
  // H-3 specific
  Double_t Z = 2.0;

  momentum = TMath::Sqrt(xx * (xx + 2.0*electronmass));
  eta = electronmass/momentum * Z*finestructure;
  Fermif = 2*pi*eta/(1.0-TMath::Exp(-2.0*pi*eta));

  Kurief = TMath::Sqrt(full / ((xx + electronmass)*momentum*Fermif));
  delete [] w;
  delete [] y;
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

//-----------------------------------------------------
// run plotcurve function
// .L h3script.C
// plotcurve(1.0e-5) 
// for a 10 meV lightest neutrino mass eigenstate
// uses neutrino mass in units of keV 
// for script convenience, not user convenience
//-----------------------------------------------------
void plotcurve(Double_t mn) // mn in [keV]
{
  Double_t Eres = 1.0e-6; // 1 meV Gauss std [keV]
  Double_t range_left;
  if (mn > 0.0)
    range_left = Qvalue - 10*mn; // [keV]
  else
    range_left = Qvalue - 1.0e-4; // [keV]


  fh3 = new TF1("h3func",betashape,1.0e-6,18.6,4);
  fh3->SetParameter(0,mn);
  fh3->SetParameter(1,mn);
  fh3->SetParameter(2,mn);
  fh3->SetParameter(3,mn);
  Double_t norm = fh3->Integral(1.0e-6,(Qvalue-mn));

  fh3n = new TF1("h3norm",normalisedbeta,range_left,Qvalue,4);
  fh3n->SetParameter(0,Eres);
  fh3n->SetParameter(1,0.0);
  fh3n->SetParameter(2,norm);
  fh3n->SetParameter(3,mn);

  kernel = new TF1("corefunc",simplesmooth,18.5,Qvalue+5*Eres,4);
  kernel->SetParameter(0,Eres);
  kernel->SetParameter(1,0.0);
  kernel->SetParameter(2,norm);
  kernel->SetParameter(3,mn);

  fh3smooth = new TF1("h3smooth",smoothed,range_left,Qvalue+5*Eres,4);
  fh3smooth->SetParameter(0,Eres);
  fh3smooth->SetParameter(1,0.0);
  fh3smooth->SetParameter(2,norm);
  fh3smooth->SetParameter(3,mn);
  fh3smooth->SetRange(range_left, Qvalue); 

  fh3->SetNpx(10000);
  fh3n->SetNpx(10000);
  fh3smooth->SetNpx(10000);
  kernel->SetNpx(10000);

  kurieh3 = new TF1("kurieplot",kuriefunction,range_left,Qvalue+5*Eres,4);
  kurieh3->SetParameter(0,Eres);
  kurieh3->SetParameter(1,0.0);
  kurieh3->SetParameter(2,norm);
  kurieh3->SetParameter(3,mn);
  kurieh3->SetNpx(10000);

  ch3 = new TCanvas("H-3","H-3",800,600);
  fh3smooth->SetTitle("Kurie function derivative;Energy [keV];dN(E)/dE");
  fh3smooth->DrawDerivative("AL");

  std::cout << " integral spectrum fraction: " << fh3n->Integral(range_left, Qvalue) << std::endl;
}
