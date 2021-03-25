#include <iostream>
#include "TH2D.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"

// Directives and declarations for namespaces and namespace members
using std::string, std::stoi, std::stod, std::cout, std::abs;

// Program-specific helper function declarations
/**
 * Identical string comparison function
 * param a: first string to be compared
 * param b: second string to be compared
 * return 1 if a & b are identical character-wise and in length, 0 otherwise
 */
bool strsame(string a, string b);

// Begin primary program function
int main(int argc, char** argv){

  // Command line option parsing variables
  bool argexists = 0;
  bool printhelp = 0;

  // Command line option storage variables
  int mode = 0;
  int Nmeas = 1;
  double par_experiment = 0.0;

  // Parse and process command line options
  for(int i = 1; i < argc; ++i){
    if(strsame(argv[i],"--gaus")){
      argexists = 1;
      mode = 0;
    }
    if(strsame(argv[i],"--exp")){
      argexists = 1;
      mode = 1;
    }
    if(strsame(argv[i],"--meas")){
      argexists = 1;
      Nmeas = stoi(argv[++i]);
    }
    if(strsame(argv[i],"--slice")){
      argexists = 1;
      par_experiment = stod(argv[++i]);
    }
    if(!argexists){
      printhelp = 1;
      cout << "Undefined option: " << argv[i] << "\n";
    }
  }
  /* Print the executable usage instructions if the user adds -h or --help,
     doesn't provide required input, or provides an undefined option */
  if(printhelp){
  cout << "\nUsage: " << argv[0] << " [options]\n"
       << "  options and descriptions:\n"
       << "   --exp             generate neyman const. for exp rate param\n"
       << "   --gaus            generate neyman const. for gaus mean (default)\n"
       << "   --meas [int]      number of measurements per experiment (1)\n"
       << "   --slice [number]  value of paramter for 1d \"slice\" hist (0)\n";
  return 0;
  }

  // Mode variables
  int a = 0;
  int b = 0;
  double c = 0.0;

  // Other helper variables
  int Nexp  = 100000;
  double par_true;
  double par_best;
  double sigma = 2;
  double x;

  if(mode == 0){
    a = -100;
    b = 100;
    c = 10.0;
  }
  else if(mode == 1){
    a = 1;
    b = 100;
    c = 100.0;
  }
  int bins = abs(b-a)+1;

  TH2D* hist2D = new TH2D("hist2D","hist2D",
        bins, double(a)/c, double(b)/c,
        bins, double(a)/c, double(b)/c);

  TH1D* hist1D = new TH1D("hist1D","hist1D",
        bins, double(a)/c, double(b)/c);

  // loop through different par_true
  for(int i = a; i <= b; i++){
    par_true = double(i)/c;

    for(int e = 0; e < Nexp; e++){
      par_best = 0.;
      // loop through measurements in experiment
      for(int m = 0; m < Nmeas; m++){
        if(mode == 0){
          x = gRandom->Gaus(par_true, sigma);
          par_best += x;
        }
        else if(mode == 1){
          x = gRandom->Exp(par_true);
          par_best += x;
        }
      }

      // our "measurement" for par best fit
      par_best = par_best / double(Nmeas);

      hist2D->Fill(par_true, par_best);
    }
  }

  // assume we do experiment and measure par_experiment

  // which bin in 2D histogram on y-axis matches measurement?
  int ibin = hist2D->GetYaxis()->FindBin(par_experiment);

  // Fill 1D histogram with slice from 2D corresponding to measurement
  for(int i = 0; i < bins; i++){
    hist1D->SetBinContent(i+1, hist2D->GetBinContent(i+1, ibin));
  }

  // normalize as probability distribution
  hist1D->Scale(1./hist1D->Integral());


  // some formating settings
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  TCanvas* can0 = (TCanvas*) new TCanvas("c0","c0",
					 450, 400);
  double hlo = 0.15;
  double hhi = 0.2;
  double hbo = 0.15;
  double hto = 0.07;
  can0->SetLeftMargin(hlo);
  can0->SetRightMargin(hhi);
  can0->SetBottomMargin(hbo);
  can0->SetTopMargin(hto);
  can0->SetGridx();
  can0->SetGridy();
  can0->cd();
  can0->Draw();

  hist2D->Draw("colz");
  hist2D->GetXaxis()->CenterTitle();
  hist2D->GetXaxis()->SetTitleSize(0.05);
  hist2D->GetYaxis()->CenterTitle();
  hist2D->GetYaxis()->SetTitleSize(0.05);
  hist2D->GetXaxis()->SetTitleOffset(1.1);
  hist2D->GetYaxis()->SetTitleOffset(1.2);
  if(mode == 0){
    hist2D->GetXaxis()->SetTitle("#mu true");
    hist2D->GetYaxis()->SetTitle("#mu meas");
  } else if(mode == 1){
    hist2D->GetXaxis()->SetTitle("#lambda true");
    hist2D->GetYaxis()->SetTitle("#lambda meas");
  }



  can0->SaveAs("2DNeyman.png");

  //////////////////////////////////////////

  TCanvas* can1 = (TCanvas*) new TCanvas("c1","c1",
					 450, 400);

  can1->SetLeftMargin(0.15);
  can1->SetRightMargin(0.1);
  can1->SetBottomMargin(hbo);
  can1->SetTopMargin(hto);
  can1->SetGridx();
  can1->SetGridy();
  can1->cd();
  can1->Draw();

  hist1D->SetLineColor(kBlue+2);
  hist1D->SetLineWidth(3);
  hist1D->SetFillColor(kBlue);
  hist1D->SetFillStyle(3004);
  hist1D->Draw("hist");

  hist1D->GetXaxis()->CenterTitle();
  hist1D->GetXaxis()->SetTitleSize(0.05);
  hist1D->GetYaxis()->CenterTitle();
  hist1D->GetYaxis()->SetTitleSize(0.05);
  hist1D->GetXaxis()->SetTitleOffset(1.1);
  hist1D->GetYaxis()->SetTitleOffset(1.3);
  if(mode == 0){
    hist1D->GetXaxis()->SetTitle("#mu");
    hist1D->GetYaxis()->SetTitle("P( #mu | #mu_{meas}, #sigma)");
  } else if(mode == 1){
    hist1D->GetXaxis()->SetTitle("#lambda");
    hist1D->GetYaxis()->SetTitle("P( #lambda | #lambda_{meas})");
  }

  can1->SaveAs("1DNeyman.png");

  return 0;
}

// Program-specific helper function definitions
bool strsame(string a, string b){
  if(a.length()==b.length()){
    int n = a.length();
    for(int i = 0; i < n; i++){
      if(a.at(i)!=b.at(i)){
        return 0;
      }
    }
    return 1;
  }
  return 0;
}
