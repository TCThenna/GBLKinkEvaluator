/******************************************************************************************
 * Original code by Hendrik Jansen                                                        *
 * 2016/7                                                                                 *
 *                                                                                        *
 *                                                                                        *
 ******************************************************************************************/

	
#include "AnaTel.h"
#include "TF1.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <utility>
#include <map>
#include <sstream>

//Root headers
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TStopwatch.h"
#include "TLatex.h"
#include "TPaveStats.h"
#include "TColor.h"
#include "TLegendEntry.h"

std::string inputDir20 = "../../analysis-20mm/output/histograms/";
int runmode = 0;
std::string inputFile = "run00";
TFile* _outputFile;
ofstream pulls20;
ofstream pulls150;

std::string whichfitter = "GBLKinkEstimator_kappa100";

TH1D *delta0 = new TH1D("temp0","dXY 0",1000, -1, 5);
TH1D *delta1 = new TH1D("temp1","dXY 1",1000, -1, 5);
TH1D *delta2 = new TH1D("temp2","dXY 2",1000, -1, 5);
TH1D *delta3 = new TH1D("temp3","dXY 3",1000, -1, 5);
TH1D *delta4 = new TH1D("temp4","dXY 4",1000, -1, 5);
TH1D *delta5 = new TH1D("temp5","dXY 5",1000, -1, 5);
TH1D *delta6 = new TH1D("temp6","dX 0 and 5",1000, -1, 5);
TH1D *delta7 = new TH1D("temp7","dX 1 and 4",1000, -1, 5);
TH1D *delta8 = new TH1D("temp8","dX 2 and 3",1000, -1, 5);
TH1D *delta9 = new TH1D("temp9","dY 0 and 5",1000, -1, 5);
TH1D *delta10 = new TH1D("temp10","dY 1 and 4",1000, -1, 5);
TH1D *delta11 = new TH1D("temp11","dY 2 and 3",1000, -1, 5);


const Int_t alucount = 8;
double aluthick[alucount] = {0.0, 0.013, 0.025, 0.05, 0.1, 0.2, 1.0, 10.0};
const Int_t planescount = 12; // planes *2

TGraphErrors* g_mean[planescount]; // for 2*6 planes
TGraphErrors* g_sigma[planescount];

TGraphErrors* g_mean_res_offset;
TGraphErrors* g_mean_res_offset_inner;

TMultiGraph* mg_obsresX = new TMultiGraph("mg_obsresX","obs_resX");
TMultiGraph* mg_obsresY = new TMultiGraph("mg_obsresY","obs_resY");

TMultiGraph* mg_effi = new TMultiGraph("mg_effi3","effi3");
// 10 grapahs for 2 set-ups * 5 energies
TGraphErrors* effi0;
TGraphErrors* effi1;
TGraphErrors* effi2;
TGraphErrors* effi3;
TGraphErrors* effi4;
TGraphErrors* effi5;
TGraphErrors* effi6;
TGraphErrors* effi7;
TGraphErrors* effi8;
TGraphErrors* effi9;


std::vector<double> v_runnumber;
std::vector<double> v_erunnumber;
std::vector<double> v_mean;
std::vector<double> v_emean;
std::vector<double> v_sigma;
std::vector<double> v_esigma;
std::vector<double> v_residual_offset;
std::vector<double> v_mean_res_offset;
std::vector<double> v_mean_res_offset_inner;

std::vector<std::vector<double>> clustersizes(10,     std::vector<double>());
std::vector<std::vector<double>> rms_clustersizes(10, std::vector<double>());

std::vector<std::vector<double>> effis(10,     std::vector<double>());
std::vector<std::vector<double>> rms_effis(10, std::vector<double>());

// Intrinsic resolution
Double_t m26_resolution =0.;
Double_t m26_res_error =0.;
Double_t global_plot_error = 0;


// Beam energy
Double_t global_beam = 0.;
Double_t global_spread = 0.;
Double_t global_thickness = 0.;

// Average noise and efficiency
Double_t avgnoise = 0.;
Double_t avgeffi = 0.;
Double_t avgnoise_error = 0.;
Double_t avgeffi_error = 0.;

// Average clustersize
Double_t avgclustersize =0;
Double_t avgclustersize_error =0;

Double_t avgmeas =0;
Double_t avgmeas_error =0;

// Plot options
Bool_t plot_residuals = false ;

// Verbosity
Bool_t verbose0 = false;
Bool_t verbose1 = false;

// Biased or unbiased
bool IsBiased = true;


// Prediction graphs for the smilie plot
#define ngraphs 6

// Total plane count
#define nplanes 6

// Plane position for the smilie plot
Double_t posx[nplanes] = { 0.0, 150.0, 300.0, 450.0, 600.0, 750.0 };
Double_t posx_error[nplanes] = { 2.5,2.5,2.5,2.5,2.5,2.5};


//char telescopebuild[50];
std::string telescopebuild;
int planedistance;

// The observed resolution & error
Double_t obsresol_x[nplanes];
Double_t obsresol_y[nplanes];
Double_t obsresol_error_x[nplanes] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
Double_t obsresol_error_y[nplanes] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
//Double_t obsresol_x_biased[nplanes];
//Double_t obsresol_y_biased[nplanes];
//Double_t obsresol_error_x_biased[nplanes] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
//Double_t obsresol_error_y_biased[nplanes] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };


// Put a cut on the threshold, 0 for all
Int_t global_thresh = 0;

// color array
Int_t color[10] = {1,2,3,4,6,7,8,9,11,12};

// What subsensor information do we want? Format: A-D for submatrix, 1-7 for clustersize, empty for all.
// e.g.: submask = "A1" or submask = "3" or submask = "C"
TString submask="";





void DoGausFit(TH1D *histo, double &mean, double &sigma)
{
  TF1 *f1 = new TF1("f1","gaus");

  float range = 2.;

  float tmp_mean= histo->GetMean();
  float tmp_rms = histo->GetRMS();

  f1->SetParameter(1,tmp_mean);
  f1->SetParameter(2,tmp_rms);
  f1->SetParLimits(2,tmp_rms*0.5,tmp_rms*2.);

  //limits
  float upper_limit = tmp_mean + range*tmp_rms;
  float lower_limit = tmp_mean - range*tmp_rms;

  histo->Fit(f1,"BQ","",lower_limit,upper_limit);

  mean = f1->GetParameter(1);
  sigma= f1->GetParameter(2);

}

void DoConstFit(TH1 *histo, double &mean, double &sigma, double low, double high)
{
  TF1 *f1 = new TF1("f1","[0]");

  float tmp_mean= histo->GetMean();

  f1->SetParameter(0,tmp_mean);

  histo->Fit(f1,"BQ","",low,high);

  mean = f1->GetParameter(0);
  sigma= f1->GetParError(0);

}

void filehelper(std::string &help, Int_t RN)
{
  help = inputDir20 + inputFile;
  //std::cout << "help: " << help << ", runnumber = " << RN << std::endl;
  //std::cout << "run number = " << RN << std::endl;
  //TString filename(help.c_str());
  if (RN <= 999)
    help+="0";
  if (RN <= 99)
    help+="0";
  if (RN <= 9)
    help+="0";
  help+= std::to_string(RN);
  help+="-";
  help+=whichfitter;
  //if(RN < 113 || RN > 700) help += "_probcut01_kappa12";
  help+=".root";
}

void CanvasSetter(TCanvas *canvas)
{
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  canvas->SetFillColor(0);
}

void DoMean(std::vector<double> vec, double &mean)
{
  mean = .0;
  for( double entry : vec)
    mean += entry;
  mean /= (float)vec.size(); 
}

void DoRMS(std::vector<double> vec, double &rms)
{
  double rms_local = .0;
  double mean_local = .0;
  DoMean(vec, mean_local);
  for( double entry : vec)
    rms_local += (entry - mean_local)*(entry-mean_local)/(float)vec.size();
  rms = sqrt(rms_local); 
}

void CSchecker(Int_t runnumber, Double_t ebeam, Int_t j)
{

  std::cout << "Start CSchecker for run " << runnumber << " with E = " << ebeam << std::endl;

  std::string help;

  filehelper(help, runnumber);

  std::cout << "Reading file: " << help << std::endl;
  TFile *tfile = new TFile(help.c_str());
  if(tfile == NULL) {
    std::cout << "   --- File does not exist, return! ---   " << std::endl;
    return;
  } else std::cout << "   --- File found ---   " << std::endl;

  // canvas
  TCanvas* canvas;
  std::string canv_name = "CSchecker_";
  canv_name += std::to_string(runnumber);
  canvas = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canvas);


  std::vector<double> cs_perPlane;
  TH1D *h_clustersize[6];

  h_clustersize[0] = (TH1D*) (tfile->Get("Fitter06/Tracks/clustersize0"))->Clone();
  h_clustersize[1] = (TH1D*) (tfile->Get("Fitter06/Tracks/clustersize1"))->Clone();
  h_clustersize[2] = (TH1D*) (tfile->Get("Fitter06/Tracks/clustersize2"))->Clone();
  h_clustersize[3] = (TH1D*) (tfile->Get("Fitter06/Tracks/clustersize3"))->Clone();
  h_clustersize[4] = (TH1D*) (tfile->Get("Fitter06/Tracks/clustersize4"))->Clone();
  h_clustersize[5] = (TH1D*) (tfile->Get("Fitter06/Tracks/clustersize5"))->Clone();

  h_clustersize[0]->GetXaxis()->SetRangeUser(1, 4.9);
  h_clustersize[1]->GetXaxis()->SetRangeUser(1, 4.9);
  h_clustersize[2]->GetXaxis()->SetRangeUser(1, 4.9);
  h_clustersize[3]->GetXaxis()->SetRangeUser(1, 4.9);
  h_clustersize[4]->GetXaxis()->SetRangeUser(1, 4.9);
  h_clustersize[5]->GetXaxis()->SetRangeUser(1, 4.9);

  for(int i = 0; i < 6; i++) cs_perPlane.push_back((h_clustersize[i]->GetMean()-0.5)); // -0.5 for correction of how ROOT calculated mean of histogramm for which the Range was changed. ROOT uses bin centres - i.e. 1.5, 2.5, 3.5, ...), but actuall clustersizes are 1, 2, 3, ...
  for(int i = 0; i < 6; i++) std::cout << "CS " << i << " = " << cs_perPlane.at(i) << std::endl;

  double mean_cs = .0;
  double rms_cs = .0;
  DoMean(cs_perPlane, mean_cs);
  DoRMS(cs_perPlane, rms_cs);

  std::cout << " Mean CS = " << mean_cs << "  RMS CS = " << rms_cs << std::endl;
  std::cout << " - run finished - " << std::endl;

  //std::vector<double> test;
  //clustersizes.push_back(test);
  (clustersizes.at(j)).push_back(mean_cs);
  //rms_clustersizes.push_back(test);
  (rms_clustersizes.at(j)).push_back(rms_cs);
}

//read effi from root files
void effi_reader(Int_t runnumber, Double_t ebeam, Int_t j)
{
  std::cout << "\nStart effi_reader for RN " << runnumber << " with E = " << ebeam << std::endl;

  //  TString filename("histograms/run00");
  std::string help;
  filehelper(help, runnumber);

  std::cout << "Reading file: " << help << std::endl;
  TFile *tfile = new TFile(help.c_str());
  if(tfile == NULL) {
    std::cout << "   --- File does not exist, return! ---   " << std::endl;
    return;
  } else  std::cout << "   --- File found ---   " << std::endl;

  // canvas
  TCanvas *canv;
  std::string canv_name = "m26effi_";
  canv_name += std::to_string(runnumber);
  canv = new TCanvas(canv_name.c_str(),canv_name.c_str(),20,20,800,800); 
  CanvasSetter(canv);
  canv->Divide(1,2);

  TProfile *p_m26effi[2];

  tfile->cd();

  std::cout << "clone histos";

  p_m26effi[0] = (TProfile*) (tfile->Get("Fitter06/Effi/effix3"))->Clone();
  p_m26effi[1] = (TProfile*) (tfile->Get("Fitter06/Effi/effiy3"))->Clone();

  double mean_eff[2];
  double rms_eff[2];
  mean_eff[0] = p_m26effi[0]->GetMean(2);
  rms_eff[0]  = p_m26effi[0]->GetRMS(2);

  p_m26effi[0]->GetXaxis()->SetRangeUser(-5,5);
  mean_eff[1] = p_m26effi[0]->GetMean(2);
  rms_eff[1]  = p_m26effi[0]->GetRMS(2);

  std::cout << " mean eff = " << mean_eff[0] << " pm " << rms_eff[0] << std::endl;

  (effis.at(j)).push_back(mean_eff[0]);

  DoConstFit(p_m26effi[0], mean_eff[0], rms_eff[0], -9, 9);
  (rms_effis.at(j)).push_back(rms_eff[0]);

  canv->cd(1);
  p_m26effi[0]->Draw(); 
  canv->cd(2);
  p_m26effi[1]->Draw(); 


  _outputFile->cd();
  canv->Modified();
  canv->Update();
  canv->Write();
  canv->Close();

  std::cout << "  -> done." << std::endl;
}


// Fitting of each file -> resolution
void fitter(Int_t runnumber, Double_t ebeam)
{
  std::cout << "\nStart fitter for RN " << runnumber << " with E = " << ebeam << std::endl;

  //  TString filename("histograms/run00");
  std::string help;
  filehelper(help, runnumber);

  std::cout << "Reading file: " << help << std::endl;
  TFile *tfile = new TFile(help.c_str());
  if(tfile == NULL) {
    std::cout << "   --- File does not exist, return! ---   " << std::endl;
    return;
  } else{
    std::cout << "   --- File found ---   " << std::endl;

  }

  // canvas
  TCanvas *canv;
  std::string canv_name = "m26fitter_";
  canv_name += std::to_string(runnumber);
  canv = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canv);
  canv->Divide(2,6);

  TCanvas *canv_pulls;
  canv_name = "m26fitter_pulls_";
  canv_name += std::to_string(runnumber);
  canv_pulls = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canv_pulls);
  canv_pulls->Divide(2,6);

  TCanvas *canv_kink_pulls;
  canv_name = "m26fitter_kink_pulls_";
  canv_name += std::to_string(runnumber);
  canv_kink_pulls = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canv_kink_pulls);
  canv_kink_pulls->Divide(2,4);





  // histos
  TH1D *h_residual_biased[planescount];
  TH1D *h_pull_biased[planescount];
  TH1D *h_kink_pull[8];
  tfile->cd();


  std::cout << "clone histos";
  // Load file
  //


  h_residual_biased[0] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx0"))->Clone();
  h_residual_biased[1] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry0"))->Clone();
  h_residual_biased[2] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx1"))->Clone();
  h_residual_biased[3] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry1"))->Clone();
  h_residual_biased[4] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx2"))->Clone();
  h_residual_biased[5] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry2"))->Clone();
  h_residual_biased[6] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx3"))->Clone();
  h_residual_biased[7] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry3"))->Clone();
  h_residual_biased[8] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx4"))->Clone();
  h_residual_biased[9] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry4"))->Clone();
  h_residual_biased[10] = (TH1D*)(tfile->Get("Fitter06/GBL/gblrx5"))->Clone();
  h_residual_biased[11] = (TH1D*)(tfile->Get("Fitter06/GBL/gblry5"))->Clone();

  h_pull_biased[0] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpx0"))->Clone();
  h_pull_biased[1] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpy0"))->Clone();
  h_pull_biased[2] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpx1"))->Clone();
  h_pull_biased[3] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpy1"))->Clone();
  h_pull_biased[4] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpx2"))->Clone();
  h_pull_biased[5] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpy2"))->Clone();
  h_pull_biased[6] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpx3"))->Clone();
  h_pull_biased[7] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpy3"))->Clone();
  h_pull_biased[8] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpx4"))->Clone();
  h_pull_biased[9] = (TH1D*) (tfile->Get("Fitter06/GBL/gblpy4"))->Clone();
  h_pull_biased[10] = (TH1D*)(tfile->Get("Fitter06/GBL/gblpx5"))->Clone();
  h_pull_biased[11] = (TH1D*)(tfile->Get("Fitter06/GBL/gblpy5"))->Clone();

  h_kink_pull[0] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx1"))->Clone();
  //h_kink_pull[1] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty1"))->Clone();
  h_kink_pull[2] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx2"))->Clone();
  //h_kink_pull[3] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty2"))->Clone();
  h_kink_pull[4] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx3"))->Clone();
  //h_kink_pull[5] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty3"))->Clone();
  h_kink_pull[6] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx4"))->Clone();
  //h_kink_pull[7] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty4"))->Clone();


  std::cout << "  - done" << std::endl;
  // Add the clustersizes 2-4 into the histos
  //
  

  //  range
  for(int i=0;i<planescount;i++){
    h_residual_biased[i]->GetXaxis()->SetRangeUser(-100.,100);
  }

  // TCanvas *canv2;
  // canv2 = new TCanvas("m26fitter2","m21fitter(",900,10,600,800); 
  // canv2->SetFillColor(0);
  // canv2->Divide(2,3);

  std::cout << "\n  Reading residuals from file...\n" << std::endl;

  std::vector<double> vec_res_sigma_biasedX;
  std::vector<double> vec_res_sigma_biasedY;
  std::vector<double> vec_pull_meanX;
  std::vector<double> vec_pull_meanY;
  std::vector<double> vec_pull_sigmaX;
  std::vector<double> vec_pull_sigmaY;
  std::vector<double> vec_kink_pull_meanX;
  std::vector<double> vec_kink_pull_sigmaX;

  // Sort into X or Y
  for(Int_t i = 0; i < nplanes*2; i++ )
  {
    if(verbose1) std::cout << "  --- i = " << i << std::endl;
    canv->cd(i+1);
    h_residual_biased[i]->Draw();
    gPad->Update();
    TPaveStats * ps = (TPaveStats *)h_residual_biased[i]->FindObject("stats");
    ps->SetX1NDC(0.65);
    ps->SetX2NDC(0.88);
    ps->SetY1NDC(0.30);
    ps->SetY2NDC(0.98);
    canv->Modified();
    canv->Update();

    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);


    int j = 6;
    if(i == 0 || i == 1)
      j = 0;
    if(i == 2 || i == 3)
      j = 1;
    if(i == 4 || i == 5)
      j = 2;
    if(i == 6 || i == 7)
      j = 3;
    if(i == 8 || i == 9)
      j = 4;
    if(i == 10 || i == 11)
      j = 5;

    TString st = "";
    char tmpstring[500];
    sprintf(tmpstring, "Sensor %1.1i", j);
    st = tmpstring;


    // Fit the residuals
    //TF1 *f1 = new TF1("f1","gaus");
    TF1 *f1_biased = new TF1("f1_biased","gaus");

    //float mean1 = h_m26[i]->GetMean();
    //float rms1 = h_m26[i]->GetRMS();

    float mean_biased = h_residual_biased[i]->GetMean();
    float rms_biased  = h_residual_biased[i]->GetRMS();

    //std::cout << "Mean = " << mean1 << "   RMS = " << rms1 << std::endl;

    //f1->SetParameter(1,mean1);
    //f1->SetParameter(2,rms1);
    //f1->SetParLimits(2,rms1*0.5,rms1*2.);
    //f1->FixParameter(2,16.);
    f1_biased->SetParameter(1,mean_biased);
    f1_biased->SetParameter(2,rms_biased);

    //h_m26[i]->Fit(f1,"QB");
    h_residual_biased[i]->Fit(f1_biased,"EMQI0");

    //Double_t mean_1 = f1->GetParameter(1);
    //Double_t sigma_1 = f1->GetParameter(2);
    //Double_t sigma_error_1 = f1->GetParError(2);


    //std::cout << " --- simga from first fit = " << sigma_1 << std::endl;

    Double_t mean_1_biased = f1_biased->GetParameter(1);
    Double_t sigma_1_biased = f1_biased->GetParameter(2);
    Double_t sigma_error_1_biased = f1_biased->GetParError(2);

    // Repeat within 2 sigma

    float nsig = 2.;
    //float low_lim = mean_1        - nsig*sigma_1;
    //float high_lim = mean_1        + nsig*sigma_1;
    //if(low_lim < -25.) low_lim = -24.9;
    //if(high_lim > 25.) high_lim = 24.9;
    float low_lim_biased = mean_1_biased        - nsig*sigma_1_biased;
    float high_lim_biased = mean_1_biased        + nsig*sigma_1_biased;
    //TF1 *f2 = new TF1("f2","gaus", low_lim, high_lim);
    TF1 *f2_biased = new TF1("f2_biased","gaus", low_lim_biased, high_lim_biased);
    //f2->SetLineWidth(2);
    //f2->SetLineStyle(1);
    //f2->SetLineColor(kBlack);
    f2_biased->SetLineWidth(2);
    f2_biased->SetLineStyle(1);
    f2_biased->SetLineColor(kBlack);
    //h_m26[i]->       Fit(f2,"EQMIR", "");
    h_residual_biased[i]->Fit(f2_biased,"EQMIR","");
    //Double_t mean_2 = f2->GetParameter(1);
    //Double_t emean_2 = f2->GetParError(1);
    Double_t mean_2_biased = f2_biased->GetParameter(1);
    Double_t emean_2_biased = f2_biased->GetParError(1);


    //Double_t sigma_2 = f2->GetParameter(2);
    //Double_t esigma_2 = f2->GetParError(2);
    Double_t sigma_2_biased = f2_biased->GetParameter(2);
    Double_t esigma_2_biased = f2_biased->GetParError(2);


    // Create telescope
    //
    //telescopebuild = "AnaTel_wide.geom";
    //AnaTel *tscope = new AnaTel(telescopebuild.c_str(), ebeam);
    //double sigma_int = 3.42e-3;
    //tscope->SetResolution(sigma_int);

    //std::cout << " measured UNBIASED width in um = " << sigma_2 << "   expected from GBL = " << tscope->GetWidth((int)i/2,0)*1e3 << std::endl; // removed 1e3
    //std::cout << " measured   BIASED width in um = " << sigma_2_biased << "   expected from GBL = " << tscope->GetWidth((int)i/2,1)*1e3 << std::endl; // removed 1e3
    //delete tscope;

    //Double_t sigma_error_2 = f2->GetParError(2) * 1000.0;
    //Double_t sigma_error_2 = 0.05 * sigma_2;
    Double_t sigma_error_2_biased = 0.05 * sigma_2_biased;

    // Now get resolution estimate
    //Double_t tel_resol = resol_estimate(sigma_2,j);

    //if (verbose0)
    //cout << " Measured residual width plane " << j << " is " <<  tel_resol << " mu m" << std::endl;

    // Fill fitted resolution into array
    if(i == 0 || i == 2 || i == 4 ||  i == 6 || i == 8 || i == 10)
    {
	obsresol_x[j] = sigma_2_biased*1.e-3;
	obsresol_error_x[j] = sigma_error_2_biased*1.e-3;
    }
    else
    {
	obsresol_y[j] = sigma_2_biased*1.e-3;
	obsresol_error_y[j] = sigma_error_2_biased*1.e-3;
    }

    if(verbose0) std::cout << " perform fits " << std::endl;
    double mean = 0.;
    double sigma = 0.;

    canv_pulls->cd(i+1);
    h_pull_biased[i]->Draw();
    gPad->Update();
    ps = (TPaveStats *)h_pull_biased[i]->FindObject("stats");
    ps->SetX1NDC(0.65);
    ps->SetX2NDC(0.88);
    ps->SetY1NDC(0.30);
    ps->SetY2NDC(0.98);
    DoGausFit(h_pull_biased[i], mean, sigma);
    canv_pulls->Modified();
    canv_pulls->Update();

    if(i%2 == 0) {
      //vec_res_sigma_unbiasedX.push_back(sigma_2);
      vec_res_sigma_biasedX.push_back(sigma_2_biased);
      vec_pull_meanX.push_back(mean);
      vec_pull_sigmaX.push_back(sigma);
    }
    if(i%2 == 1) {
      //vec_res_sigma_unbiasedY.push_back(sigma_2);
      vec_res_sigma_biasedY.push_back(sigma_2_biased);
      vec_pull_meanY.push_back(mean);
      vec_pull_sigmaY.push_back(sigma);
    }

    mean = 0;
    sigma = 0;
    if( i > 1 && i < 10){
      if(i%2 == 0) {
        canv_kink_pulls->cd(i+1-2);    
        h_kink_pull[i]->Draw();
        gPad->Update();
        ps = (TPaveStats *)h_kink_pull[i]->FindObject("stats");
        ps->SetX1NDC(0.65);
        ps->SetX2NDC(0.88);
        ps->SetY1NDC(0.30);
        ps->SetY2NDC(0.98);

	DoGausFit(h_kink_pull[i-2], mean, sigma);
	vec_kink_pull_meanX.push_back(mean);
	vec_kink_pull_sigmaX.push_back(sigma);
      }
      canv_kink_pulls->Modified();
      canv_kink_pulls->Update();
    }

    //if(i%2 == 1) {
    //  vec_kinks_meanY.push_back(mean);
    //  vec_kinks_sigmaY.push_back(sigma);
    //}

    mean = 0;
    sigma = 0;


    //if(i%2 == 1) {
    //  vec_kinks_meanY.push_back(mean);
    //  vec_kinks_sigmaY.push_back(sigma);
    //}



    if (i == 0) {
      std::string s_runnumber = "Run number ";
      s_runnumber += std::to_string(runnumber);
      TLatex* tex = new TLatex(-0.04,51,s_runnumber.c_str());
      tex->SetTextSize(0.15);
      tex->SetLineWidth(2);
      tex->Draw();
    }

    //Fill mean vector
    //std::cout << "PB, j = " << j << " mean = " << mean_2 << " emean = " << emean_2 << std::endl;
      v_mean.push_back(mean_2_biased);
      v_emean.push_back(emean_2_biased);
      v_sigma.push_back(sigma_2_biased);
      v_esigma.push_back(esigma_2_biased);

    delete f1_biased;
    delete f2_biased;




  }

  // cout pulls

  double avg_pull_mean = .0; 
  DoMean(vec_pull_meanX,avg_pull_mean);
  std::cout << " vec_pull_meanX = ";
  for(unsigned int i = 0; i<6; i++) std::cout << vec_pull_meanX.at(i) << " ";
  std::cout  << "   avgX = " << avg_pull_mean  << std::endl;

  avg_pull_mean = .0;
  std::cout << " vec_pull_meanY = ";
  for(unsigned int i = 0; i<6; i++){
    std::cout << vec_pull_meanY.at(i) << " ";
    avg_pull_mean += vec_pull_meanY.at(i)/6.;
  }
  std::cout  << "   avgY = " << avg_pull_mean << std::endl;

  double avg_pull_sigma = .0;
  double rms_pull_sigma = .0;
  std::cout << " vec_pull_sigmaX = ";
  for(unsigned int i = 0; i<6; i++){
    std::cout << vec_pull_sigmaX.at(i) << " ";
    avg_pull_sigma += vec_pull_sigmaX.at(i)/6.;
  }
  std::cout  << "   avgX = " << avg_pull_sigma;
  DoRMS(vec_pull_sigmaX,rms_pull_sigma);
  std::cout  << "   rmsX = " << rms_pull_sigma << std::endl;

  avg_pull_sigma = .0;
  rms_pull_sigma = .0;
  std::cout << " vec_pull_sigmaY = ";
  for(unsigned int i = 0; i<6; i++){
    std::cout << vec_pull_sigmaY.at(i) << " ";
    avg_pull_sigma += vec_pull_sigmaY.at(i)/6.;
  }
  std::cout  << "   avgY = " << avg_pull_sigma;
  DoRMS(vec_pull_sigmaY,rms_pull_sigma);
  std::cout  << "   rmsY = " << rms_pull_sigma << std::endl;

  double avg_res_sigma = .0;
  double rms_res_sigma = .0;
  std::cout << " vec_res_sigma_biasedX = ";
  for(unsigned int i = 0; i<6; i++){
    std::cout << vec_res_sigma_biasedX.at(i) << " ";
    avg_res_sigma += vec_res_sigma_biasedX.at(i)/6.;
  }
  //std::cout  << "   avgX = " << avg_res_sigma;
  //DoRMS(vec_res_sigma_biasedX,rms_res_sigma);
  //std::cout  << "   rmsX = " << rms_res_sigma << std::endl;


  float avg_kink_pull_sigma = .0;
  std::cout << " vec_kink_pull_sigmaX = ";
  for(unsigned int i = 0; i<4; i++){
    std::cout << vec_kink_pull_sigmaX.at(i) << " ";
    avg_kink_pull_sigma += vec_kink_pull_sigmaX.at(i)/4.;
  }
  std::cout  << " avg kink pull X = " << avg_kink_pull_sigma << std::endl;

  /*avg_kinks_sigma = .0;
    std::cout << " vec_kinks_sigmaY = ";
    for(unsigned int i = 0; i<6; i++){
    std::cout << vec_kinks_sigmaY.at(i) << " ";
    avg_kinks_sigma += vec_kinks_sigmaY.at(i)/6.;
    }
    std::cout  << " avg kink Y = " << avg_kinks_sigma << std::endl;
    */



  v_runnumber.push_back((double)runnumber);
  v_erunnumber.push_back(0.0);

  TString outputname="run_";
  outputname+=runnumber;
  outputname+="_";
  outputname+=submask;
  //canv->Print("pics/"+outputname+"_measured_residuals.eps");
  _outputFile->cd();
  canv->Write();
  canv->Close();
  canv_pulls->Write();
  canv_pulls->Close();
  canv_kink_pulls->Write();
  canv_kink_pulls->Close();

  // Run mean offset calculator with energy and measured resolutions
  //run_res_offset( ebeam, obsresol_x, obsresol_error_x, obsresol_y, obsresol_error_y );

  // Run minimization with energy and measured resolutions
  //run_global( ebeam, obsresol_x, obsresol_error_x, obsresol_y, obsresol_error_y );
  //ebeam = global_beam;

  Double_t positions[6] = {0.,1.,2.,3.,4.,5.};
  TGraph* g_obsresidualX = new TGraph(6,positions, obsresol_x);
  TGraph* g_obsresidualY = new TGraph(6,positions, obsresol_y);

  //std::cout << " runnumber%20 = " << runnumber%10 << std::endl;
  g_obsresidualX->SetMarkerStyle(runnumber%10 + 20);
  g_obsresidualX->SetMarkerColor(runnumber%9+1);
  g_obsresidualX->SetMarkerSize(1.2);
  g_obsresidualY->SetMarkerStyle(runnumber%10 + 20);
  g_obsresidualY->SetMarkerColor(runnumber%9+1);
  g_obsresidualY->SetMarkerSize(1.2);


  mg_obsresX->Add(g_obsresidualX,"p");
  mg_obsresY->Add(g_obsresidualY,"p");

  delete tfile;

  delete canv;
  delete canv_pulls;
  delete canv_kink_pulls;


}

// Fitting of each file -> noise and efficiency
void noise(Int_t runnumber)
{
  std::string help;
  help = inputDir20 + inputFile;
  TString filename(help.c_str());
  if (runnumber <= 999)
    filename+="0";
  if (runnumber <= 99)
    filename+="0";
  filename+=runnumber;
  filename+="-fitter.root";
  TFile *tfile = new TFile(filename);
  TH1D *h_m26[planescount];
  TH1D *h_m26_eff[planescount];
  tfile->cd();

  std::cout << " - This is noise() for run " << runnumber << "." << std::endl;

  // Load file
  h_m26[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTnoiseX"))->Clone();
  h_m26[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTnoiseY"))->Clone();
  h_m26[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTnoiseX"))->Clone();
  h_m26[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTnoiseY"))->Clone();
  h_m26[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTnoiseX"))->Clone();
  h_m26[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTnoiseY"))->Clone();
  h_m26[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTnoiseX"))->Clone();
  h_m26[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTnoiseY"))->Clone();
  h_m26[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTnoiseX"))->Clone();
  h_m26[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTnoiseY"))->Clone();
  h_m26[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTnoiseX"))->Clone();
  h_m26[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTnoiseY"))->Clone();

  h_m26_eff[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTeffiX"))->Clone();
  h_m26_eff[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTeffiY"))->Clone();
  h_m26_eff[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTeffiX"))->Clone();
  h_m26_eff[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTeffiY"))->Clone();
  h_m26_eff[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTeffiX"))->Clone();
  h_m26_eff[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTeffiY"))->Clone();
  h_m26_eff[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTeffiX"))->Clone();
  h_m26_eff[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTeffiY"))->Clone();
  h_m26_eff[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTeffiX"))->Clone();
  h_m26_eff[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTeffiY"))->Clone();
  h_m26_eff[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTeffiX"))->Clone();
  h_m26_eff[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTeffiY"))->Clone();



  // Add the values for all 6 sensor planes, in x and y and divide by 12

  Double_t noisevalue = 0;
  Double_t noisevalue_error = 0;
  Double_t efficiency = 0;
  Double_t efficiency_error = 0;
  for(int i=4;i<8;i++)
  {
    noisevalue += h_m26[i]->GetMean(2);
    noisevalue_error += h_m26[i]->GetMeanError(2);
    efficiency += h_m26_eff[i]->GetMean(2);
    efficiency_error += h_m26_eff[i]->GetMeanError(2);
  }
  avgnoise = noisevalue / 4.0;
  avgnoise_error = noisevalue_error / 4.0;
  avgeffi = efficiency / 4.0;
  avgeffi_error = efficiency_error / 4.0;



  std::cout << "   - run done" << std::endl;

}



void FillGraph(){
  int size = (v_mean.size()) / (2*nplanes);
  std::cout << " size of vect = " << size << std::endl;


  TGraphErrors* g_mean0 = new TGraphErrors(size);
  TGraphErrors* g_mean1 = new TGraphErrors(size);
  TGraphErrors* g_mean2 = new TGraphErrors(size);
  TGraphErrors* g_mean3 = new TGraphErrors(size);
  TGraphErrors* g_mean4 = new TGraphErrors(size);
  TGraphErrors* g_mean5 = new TGraphErrors(size);
  TGraphErrors* g_mean6 = new TGraphErrors(size);
  TGraphErrors* g_mean7 = new TGraphErrors(size);
  TGraphErrors* g_mean8 = new TGraphErrors(size);
  TGraphErrors* g_mean9 = new TGraphErrors(size);
  TGraphErrors* g_mean10 = new TGraphErrors(size);
  TGraphErrors* g_mean11 = new TGraphErrors(size);

  g_mean[0] = g_mean0;
  g_mean[1] = g_mean1;
  g_mean[2] = g_mean2;
  g_mean[3] = g_mean3;
  g_mean[4] = g_mean4;
  g_mean[5] = g_mean5;
  g_mean[6] = g_mean6;
  g_mean[7] = g_mean7;
  g_mean[8] = g_mean8;
  g_mean[9] = g_mean9;
  g_mean[10]= g_mean10;
  g_mean[11]= g_mean11;

  std::string name;
  for(int i = 0; i<2*nplanes; i++)
  {
    name = "Residual Mean" + std::to_string(i)+" per RN";
    g_mean[i]->SetNameTitle(name.c_str(), name.c_str());
  }

  g_mean0->SetLineColor(kRed);
  g_mean1->SetLineColor(kRed);
  g_mean10->SetLineColor(kRed);
  g_mean11->SetLineColor(kRed);


  TGraphErrors* g_sigma0 = new TGraphErrors(size);
  TGraphErrors* g_sigma1 = new TGraphErrors(size);
  TGraphErrors* g_sigma2 = new TGraphErrors(size);
  TGraphErrors* g_sigma3 = new TGraphErrors(size);
  TGraphErrors* g_sigma4 = new TGraphErrors(size);
  TGraphErrors* g_sigma5 = new TGraphErrors(size);
  TGraphErrors* g_sigma6 = new TGraphErrors(size);
  TGraphErrors* g_sigma7 = new TGraphErrors(size);
  TGraphErrors* g_sigma8 = new TGraphErrors(size);
  TGraphErrors* g_sigma9 = new TGraphErrors(size);
  TGraphErrors* g_sigma10 = new TGraphErrors(size);
  TGraphErrors* g_sigma11 = new TGraphErrors(size);


  g_sigma[0] = g_sigma0;
  g_sigma[1] = g_sigma1;
  g_sigma[2] = g_sigma2;
  g_sigma[3] = g_sigma3;
  g_sigma[4] = g_sigma4;
  g_sigma[5] = g_sigma5;
  g_sigma[6] = g_sigma6;
  g_sigma[7] = g_sigma7;
  g_sigma[8] = g_sigma8;
  g_sigma[9] = g_sigma9;
  g_sigma[10]= g_sigma10;
  g_sigma[11]= g_sigma11;


  g_sigma0->SetLineColor(kRed);
  g_sigma1->SetLineColor(kRed);
  g_sigma10->SetLineColor(kRed);
  g_sigma11->SetLineColor(kRed);

  g_mean_res_offset = new TGraphErrors(size);
  g_mean_res_offset->SetNameTitle("MeanResOffset","MeanResOffset");
  g_mean_res_offset_inner = new TGraphErrors(size);
  g_mean_res_offset_inner->SetNameTitle("MeanResOffset_inner","MeanResOffset_inner");

  for(int i = 0; i<2*nplanes; i++)
  {
    name = "Residual width" + std::to_string(i) + " per RN";
    g_sigma[i]->SetNameTitle(name.c_str(), name.c_str());
    g_sigma[i]->SetMarkerStyle(i%10 + 20);
    g_sigma[i]->SetMarkerColor(i%9+1);
    g_sigma[i]->SetMarkerSize(1.2);
    g_sigma[i]->SetMarkerStyle(i%10 + 20);
    g_sigma[i]->SetMarkerColor(i%9+1);
    g_sigma[i]->SetMarkerSize(1.2);
  }

  std::cout << " SetPoints and Erros" << std::endl;

  for(int i = 0; i < size; i++){
    //std::cout << "i = " << i << std::endl;
    for(int j = 0; j < (2*nplanes); j++){
      //if(j==0) std::cout << i << " " << j << " " << v_runnumber.at(i) << " " << v_mean.at(i*2*nplanes +j) << " " << v_erunnumber.at(i) << " " <<  v_emean.at(i*2*nplanes+j) << std::endl; 
      g_mean[j]->SetPoint     (i,  v_runnumber[i],  v_mean[i*2*nplanes+j]); 
      g_mean[j]->SetPointError(i, v_erunnumber[i], v_emean[i*2*nplanes+j]); 

      g_sigma[j]->SetPoint     (i,  v_runnumber[i],  v_sigma[i*2*nplanes+j]); 
      g_sigma[j]->SetPointError(i, v_erunnumber[i], v_esigma[i*2*nplanes+j]); 
    }

    //g_mean_res_offset->SetPoint     (i,  v_runnumber[i], v_mean_res_offset.at(i)); 
    //g_mean_res_offset->SetPointError(i, v_erunnumber[i], 0.0); 

    //g_mean_res_offset_inner->SetPoint     (i,  v_runnumber[i], v_mean_res_offset_inner.at(i)); 
    //g_mean_res_offset_inner->SetPointError(i, v_erunnumber[i], 0.0); 

  }

  std::cout << "Write to File" << std::endl;
  _outputFile->cd();

  for(int i = 0; i<2*nplanes; i++)
  {
    g_mean[i]->Draw("AP");
    g_mean[i]->Write();
  }


  for(int i = 0; i<2*nplanes; i++)
  {
    g_sigma[i]->Draw("AP");
    g_sigma[i]->Write();
  }

  //g_mean_res_offset->Write();
  //g_mean_res_offset_inner->Write();

  std::cout << "MGs" << std::endl;
  TMultiGraph* mg_mean = new TMultiGraph("mg_mean","mg_mean");

  for(int i = 0; i<2*nplanes; i++)
  {
    mg_mean->Add(g_mean[i],"p");
  }

  mg_mean->Draw("ap");
  mg_mean->GetXaxis()->SetTitle("runnumber");
  mg_mean->GetYaxis()->SetTitle("measured mean X, Y");
  mg_mean->SetTitle("measured mean X, Y");

  mg_mean->Write();

  TMultiGraph* mg_sigma = new TMultiGraph("mg_sigma","mg_sigma");

  for(int i = 0; i<2*nplanes; i++)
  {
    mg_sigma->Add(g_sigma[i],"p");
  }

  mg_sigma->Draw("ap");
  mg_sigma->GetXaxis()->SetTitle("runnumber");
  mg_sigma->SetName("measured residual X, Y");
  mg_sigma->GetYaxis()->SetTitle("measured residual X, Y");

  mg_sigma->Write();

  mg_obsresX->Write();
  mg_obsresY->Write();

  delete mg_obsresX;
  delete mg_obsresY;

}

// Let's go!
int main(int argc, char *argv[0])
{
  std::cout << "Usage: " << argv[0] << "<runmode> <data input path>" << std::endl;
  std::cout << "0 for all" << std::endl;
  std::cout << "2 for threshold" << std::endl;
  std::cout << "9 for global_thres=6" << std::endl;
  std::cout << "11 for read pulls and write to file" << std::endl;
  std::cout << "111 mean number of fired pixels" << std::endl;
  std::cout << "16 Plot predicted and measured residual width versus z position (smiley)" << std::endl;
  std::cout << "17 read effi plots and create plots" << std::endl;
  std::cout << "18 make reso plots a.f.o. dz_DUT" << std::endl;
  std::cout << "19 make reso plots a.f.o. eps_DUT" << std::endl;
  std::cout << "20 DIY tscope\n" << std::endl;

  //gROOT->SetBatch();
  //gSystem->Load("libMinuit");
  //gSystem->Load("lib/libGBL.so");
  //gSystem->AddIncludePath("include/");
  gROOT->SetStyle("Plain");


  // arguments

  if( argc == 3 ) {
    runmode = atoi(argv[1]);
    inputDir20 = argv[2];
  }
  else if( argc == 1 ) {
    std::cout << "Using default arguments." << std::endl;
  }
  else if( argc > 3 ) {
    std::cout << "Too many arguments supplied." << std::endl;
    return 1;
  }
  else if( argc < 3 ) {
    std::cout << "Too less arguments supplied." << std::endl;
    return 1;
  }
  std::cout << "Run mode is " << runmode << std::endl;
  std::cout << "Data input path is " << inputDir20 << std::endl;

  _outputFile = new TFile("../bin/output.root", "RECREATE");
  _outputFile->cd();




  // Initial telescope constructor name, AnaTel_wide.geom for 150mm, AnaTel_narrow.geom for 20mm data
  //char telescopebuild[50];
  //int tempint = 0;

  //telescopebuild = "AnaTel_wide.geom";


  // Initialises the run vectors
  std::vector<int> run_5_20;
  std::vector<int> run_4_20;
  std::vector<int> run_3_20;
  std::vector<int> run_2_20;
  std::vector<int> run_1_20;

  // Fills the run vectors with runnumbers
  if(true)
  {

    // 20 mm
    // energy 5 GeV:
    run_5_20.push_back(5);	// alu     0
    run_5_20.push_back(57);	// alu    13
    run_5_20.push_back(7);	// alu    25
    run_5_20.push_back(20);	// alu    50
    run_5_20.push_back(22);	// alu   100
    run_5_20.push_back(36);	// alu   200
    run_5_20.push_back(38);	// alu  1000
    run_5_20.push_back(50);	// alu 10000

    // energy 4 GeV:
    run_4_20.push_back(4);	// alu     0
    run_4_20.push_back(56);	// alu    13
    run_4_20.push_back(9);	// alu    25
    run_4_20.push_back(19);	// alu    50
    run_4_20.push_back(23);	// alu   100
    run_4_20.push_back(35);	// alu   200
    run_4_20.push_back(39);	// alu  1000
    run_4_20.push_back(49);	// alu 10000

    // energy 3 GeV:
    run_3_20.push_back(3);	// alu     0
    run_3_20.push_back(55);	// alu    13
    run_3_20.push_back(11);	// alu    25
    run_3_20.push_back(18);	// alu    50
    run_3_20.push_back(24);	// alu   100
    run_3_20.push_back(33);	// alu   200
    run_3_20.push_back(41);	// alu  1000
    run_3_20.push_back(48);	// alu 10000

    // energy 2 GeV:
    run_2_20.push_back(1);	// alu     0
    run_2_20.push_back(54);	// alu    13
    run_2_20.push_back(12);	// alu    25
    run_2_20.push_back(17);	// alu    50
    run_2_20.push_back(26);	// alu   100
    run_2_20.push_back(32);	// alu   200
    run_2_20.push_back(42);	// alu  1000
    run_2_20.push_back(47);	// alu 10000

    // energy 1 GeV:
    run_1_20.push_back(2);	// alu     0
    run_1_20.push_back(52);	// alu    13
    run_1_20.push_back(13);	// alu    25
    run_1_20.push_back(15);	// alu    50
    run_1_20.push_back(28);	// alu   100
    run_1_20.push_back(30);	// alu   200
    run_1_20.push_back(44);	// alu  1000
    run_1_20.push_back(46);	// alu 10000

  }



  // Open config file
  //std::string s_config = "../conf/config.txt";
  //ifstream conf(s_config.c_str());
  //if(!conf)
  //{
  //  std::cout << "../conf/config.txt file missing!" << std::endl;
  //  return 1;
  //} else
  //{
  //  conf >> runmode;
  //}

  if (runmode == 20)
  {

    std::cout << " DIY tscope" << std::endl;

    double p = 5.; // GeV

    std::cout << "p " << p << " GeV" << std::endl;

    telescopebuild = "AnaTel_narrow_DUT.geom";

    AnaTel *tscope = new AnaTel(telescopebuild.c_str(), p);
    tscope->SetPlane(0,   0, 0.055,          93.66, 0.00324);
    tscope->SetPlane(1,  20, 0.055,          93.66, 0.00324);
    tscope->SetPlane(2,  40, 0.055,          93.66, 0.00324);
    tscope->SetPlane(3,  60,    -1,          0.01, -1);
    tscope->SetPlane(4,  80, 0.055,          93.66, 0.00324);
    tscope->SetPlane(5, 100, 0.055,          93.66, 0.00324);
    tscope->SetPlane(6, 120, 0.055,          93.66, 0.00324);
    tscope->SetKappa(.753);
    //std::cout << " - Tscope created, w/ DUT" << std::endl;

    double tmp; 
    for(int i = 0; i < 7; i++){
      //if(i != 3) {
      //tmp = tscope->GetWidthGBL(i, 1); // at DUT
      //std::cout << " biased width at " << i << " = " << tmp << std::endl;
      //} else {
      //tmp = tscope->GetPointingResGBL(i, 0); // at DUT
      //std::cout << " unbiased track res at plane " << i << " = " << tmp << std::endl;
      tmp = tscope->GetPointingResGBL(i, 1); // at DUT
      std::cout << " biased track res at plane " << i << " = " << tmp << std::endl;
      //}
    }


  }



  if (runmode == 19)
  {

    std::cout << " make reso plots a.f.o. eps_DUT" << std::endl;

    double p = 120.; // GeV

    const int msteps = 3;
    double dz_DUT[msteps] = {20, 60, 100};

    //const int nsteps = 300; // number of steps for eps_DUT
    const int nsteps = 894; // number of steps for eps_DUT 120
    //double eps_DUT      = 0.0005; //  
    double eps_DUT      = 0.005; //  
    //double delta_epsDUT = 0.001; //  per step
    double delta_epsDUT = 0.001; //  per step 120

    int colorr[3] = {1,2,4};

    //double reso[nsteps][msteps];
    std::vector< std::vector<double> > reso;

    std::cout << "p " << p << " GeV" << std::endl;

    for( int h = 0; h < 2; h++){
      telescopebuild = "AnaTel_narrow_DUT.geom"; // 20
      double dz = 20.;
      if(h == 1) {
	telescopebuild = "AnaTel_wide_DUT.geom"; // 150
	dz = 150.;
      }


      for( int i = 0; i < msteps; i++){



	std::vector<double> reso_n;
	//eps_DUT = 0.0005;
	eps_DUT = 0.006; // 120
	for( int j = 0;  j < nsteps; j++){


	  AnaTel *tscope = new AnaTel(telescopebuild.c_str(), p);
	  double position = 0.;
	  for( int k = 0; k < 7; k++){
	    if (k != 3) tscope->SetPlane(k, position, 0.055,   93.66, 0.00324);
	    else        tscope->SetPlane(k, position,    -1, eps_DUT, 0.0);

	    if (k == 2 || k == 3) position += dz_DUT[i];
	    else position += dz;
	  }
	  tscope->SetKappa(.864);
	  std::cout << " - Tscope created, w/ DUT" << std::endl;

	  double tmp = tscope->GetPointingResGBL(3, 0); // at DUT
	  std::cout << " reso = " << tmp << std::endl;
	  reso_n.push_back(tmp); // [um]

	  eps_DUT += delta_epsDUT;

	} // for eps_DUT
	reso.push_back(reso_n);

      } // for dz_DUT
    } // for dz = 20 , 150

    std::cout << " --- Start plotting --- " << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1",800,800);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->cd();
    c1->cd();
    gPad->SetLogx();
    //gPad->SetLogy(); // 120

    TMultiGraph *mg = new TMultiGraph();
    mg->SetName("Reso");
    mg->SetTitle(""); 

    TGraphErrors *gr[6];
    for(unsigned int i = 0; i < reso.size(); i++) {
      std::cout << "i = " << i << std::endl;
      //gr->SetPoint(0,0,0);
      gr[i] = new TGraphErrors(nsteps);
      if (i<3) gr[i]->SetLineStyle(7);
      gr[i]->SetLineWidth(3);
      gr[i]->SetLineColor(colorr[i%3]);

      //eps_DUT = 0.0005;
      eps_DUT = 0.006; // 120

      for(unsigned int j = 0; j < reso.at(i).size(); j++) {

	gr[i]->SetPoint(j, eps_DUT, reso.at(i).at(j)*1e3);
	eps_DUT += delta_epsDUT;

      }

      gr[i]->GetXaxis()->SetRangeUser(0,200);
      mg->Add(gr[i], "l");
    }

    mg->Draw("al");

    mg->GetXaxis()->SetTitle("#epsilon_{DUT}");
    mg->GetYaxis()->SetTitle("track resolution at DUT [#mum]");
    mg->GetXaxis()->SetLimits(0.005,1.0);
    mg->GetXaxis()->SetTitleOffset(1.14);
    mg->GetYaxis()->SetTitleOffset(1.44);
    //mg->SetMinimum(1.01);
    //mg->SetMaximum(10);  
    mg->SetMinimum(0.0); // 120
    mg->SetMaximum(3.5);   // 120
    mg->GetXaxis()->SetMoreLogLabels();
    mg->GetYaxis()->SetMoreLogLabels();
    mg->GetXaxis()->SetNoExponent();

    TLegend *leg = new TLegend(0.54,0.15,0.92,0.35);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetShadowColor(0);
    leg->SetLineColor(0);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(2);
    leg->SetMargin(0.4);

    //leg->SetHeader(" dz =  20 mm      dz = 150 mm");
    //leg->AddEntry("NULL","dz =  20 mm","");
    //leg->AddEntry("NULL","dz = 150 mm","");
    leg->AddEntry("null",  "dz_{DUT}","");
    leg->AddEntry("null",  "dz_{DUT}","");
    leg->AddEntry(gr[0], "  20 mm","l");
    leg->AddEntry(gr[3], "  20 mm","l");
    leg->AddEntry(gr[1], "  60 mm","l");
    leg->AddEntry(gr[4], "  60 mm","l");
    leg->AddEntry(gr[2], "100 mm","l");
    leg->AddEntry(gr[5], "100 mm","l");

    TLegend *leg3 = new TLegend(0.53,0.35,1.18,0.4);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetShadowColor(0);
    leg3->SetLineColor(0);
    leg3->SetFillColor(kWhite);
    leg3->AddEntry("NULL","dz =  20 mm    dz = 150 mm","h");

    leg->Draw();
    leg3->Draw();

    TLine *linee = new TLine(0.00071, 1.2, 0.00071, 7.5);
    linee->SetLineWidth(2);
    linee->SetLineStyle(5);
    linee->Draw();

    TLatex *tex = new TLatex(0.0013,7.0,"p = 5 GeV");
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(0.00062,7.987388,"#epsilon_{M26}");
    tex->SetTextSize(0.04);
    tex->SetLineWidth(2);
    tex->Draw();  

    tex = new TLatex(0.1768824,2.518734,"#circ");
    tex->SetTextColor(4);
    tex->SetTextSize(0.1);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.3166828,2.286983,"#circ");
    tex->SetTextColor(2);
    tex->SetTextSize(0.1);
    tex->SetLineWidth(2);
    tex->Draw(); 
    tex = new TLatex(0.3,2.13057,"#circ");
    tex->SetTextColor(1);
    tex->SetTextSize(0.1);
    tex->SetLineWidth(2);
    //tex->Draw(); 

    // ------------>Primitives in pad: Canvas_1_n2_1
    TPad *Canvas_1_n2_1 = new TPad("Canvas_1_n2_1", "newpad2",0.04396985,0.02971576,0.2135678,0.09431525);
    Canvas_1_n2_1->Draw();
    Canvas_1_n2_1->cd();
    Canvas_1_n2_1->Range(0,0,1,1);
    Canvas_1_n2_1->SetFillColor(0);
    Canvas_1_n2_1->SetBorderMode(0);
    Canvas_1_n2_1->SetBorderSize(2);
    Canvas_1_n2_1->SetFrameBorderMode(0);
    Canvas_1_n2_1->Modified();
    // ------------>Primitives in pad: Canvas_1_n2_2
    TPad *Canvas_1_n2_2 = new TPad("Canvas_1_n2_2", "newpad",0.1155779,0.0749354,0.1457286,0.124031);
    Canvas_1_n2_2->Draw();
    Canvas_1_n2_2->cd();
    Canvas_1_n2_2->Range(0,0,1,1);
    Canvas_1_n2_2->SetFillColor(0);
    Canvas_1_n2_2->SetBorderMode(0);
    Canvas_1_n2_2->SetBorderSize(2);
    Canvas_1_n2_2->SetFrameBorderMode(0);
    Canvas_1_n2_2->Modified();
    c1->cd();
    tex = new TLatex(0.01,0.0231331,"0");
    tex->SetTextSize(0.04005168);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.013777148,3.033107,"p = 120 GeV");
    tex->SetLineWidth(2);
    tex->Draw();

    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
    c1->ToggleToolBar();

    c1->Update();
    c1->Modified();

    c1->Write();
    mg->Write();

    _outputFile->Close();

    std::cout << "Written to file: " << std::endl;
    std::cout << " Done." << std::endl;


  }

  if (runmode == 18)
  {

    std::cout << " make reso plots a.f.o. dz_DUT" << std::endl;

    double p = 5.; // GeV

    const int msteps = 4;
    double X0_DUT_frac[msteps] = {0.03, 0.01, 0.003, 0.001};

    const int nsteps = 200; // number of steps for dz_DUT
    double dz_DUT = 1; // [mm] 
    double delta_dzDUT = 1; // [mm] per step

    //double reso[nsteps][msteps];
    std::vector< std::vector<double> > reso;

    std::cout << "p " << p << " GeV" << std::endl;

    for( int h = 0; h < 2; h++){
      telescopebuild = "AnaTel_narrow_DUT.geom"; // 20
      double dz = 20.;
      if(h==0) std::cout << " dz = 20" << std::endl;
      if(h==1) std::cout << " dz = 150" << std::endl;

      if(h == 1) {
	telescopebuild = "AnaTel_wide_DUT.geom"; // 150
	dz = 150.;
      }


      for( int i = 0; i < msteps; i++){



	std::vector<double> reso_n;
	dz_DUT = 1;
	for( int j = 0;  j < nsteps; j++){
	  if(i==3) std::cout << "j = " << j << " :   ";


	  AnaTel *tscope = new AnaTel(telescopebuild.c_str(), p);
	  double position = 0.;
	  for( int k = 0; k < 7; k++){
	    if (k != 3) tscope->SetPlane(k, position, 0.055,          93.66, 0.00324);
	    //if (k != 3) tscope->SetPlane(k, position, 0.055,          93.66, 0.003324);
	    //if (k != 3) tscope->SetPlane(k, position, 0.055,          93.66, 0.003156);
	    else        tscope->SetPlane(k, position,    -1, X0_DUT_frac[i], 0.0);

	    if (k == 2 || k == 3) position += dz_DUT;
	    else position += dz;
	  }
	  tscope->SetKappa(.864);
	  //std::cout << " - Tscope created, w/ DUT" << std::endl;

	  double tmp = tscope->GetPointingResGBL(3, 0); // at DUT
	  if(i==3) std::cout << " reso = " << tmp << std::endl;
	  reso_n.push_back(tmp); // [um]

	  dz_DUT += delta_dzDUT;

	} // for dz_DUT
	reso.push_back(reso_n);

      } // for eps_DUT
    } // for dz = 20 , 150

    std::cout << " --- Start plotting --- " << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1",800,800);
    c1->Range(0,0,1,1);
    c1->SetFillColor(0);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameBorderMode(0);
    c1->SetRightMargin(0.05);
    c1->SetLeftMargin(0.15);
    c1->SetTickx(1);
    c1->SetTicky(1);
    c1->cd();
    c1->cd();

    TMultiGraph *mg = new TMultiGraph();
    mg->SetName("Reso");
    mg->SetTitle(""); 

    const int line[msteps] = {2,7,9,1};

    TGraphErrors *gr[8];
    for(unsigned int i = 0; i < reso.size(); i++) {
      std::cout << "i = " << i << std::endl;
      //gr->SetPoint(0,0,0);
      gr[i] = new TGraphErrors(nsteps);
      gr[i]->SetLineStyle(line[i%4]);
      gr[i]->SetLineWidth(3);
      if (i<4) gr[i]->SetLineColor(2);

      dz_DUT = 1.;

      for(unsigned int j = 0; j < reso.at(i).size(); j++) {

	gr[i]->SetPoint(j, dz_DUT, reso.at(i).at(j)*1e3);
	dz_DUT += delta_dzDUT;

      }

      gr[i]->GetXaxis()->SetRangeUser(0,200);
      mg->Add(gr[i], "l");
    }

    mg->Draw("al");

    mg->GetXaxis()->SetTitle("dz_{DUT} [mm]");
    mg->GetYaxis()->SetTitle("track resolution at z_{DUT} [#mum]");
    mg->GetXaxis()->SetLimits(0,200);
    mg->GetXaxis()->SetTitleOffset(1.14);
    mg->GetYaxis()->SetTitleOffset(1.44);
    mg->SetMinimum(0);
    mg->SetMaximum(10);
    //mg->SetMaximum(3); // 120
    TLegend *leg = new TLegend(0.54,0.15,0.92,0.35);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetShadowColor(0);
    leg->SetLineColor(0);
    leg->SetFillColor(kWhite);
    leg->SetNColumns(2);
    leg->SetMargin(0.4);

    //leg->SetHeader(" dz =  20 mm      dz = 150 mm");
    //leg->AddEntry("NULL","dz =  20 mm","");
    //leg->AddEntry("NULL","dz = 150 mm","");
    leg->AddEntry("null",  "#epsilon_{DUT}","");
    leg->AddEntry("null",  "#epsilon_{DUT}","");
    leg->AddEntry(gr[0], "0.03","l");
    leg->AddEntry(gr[4], "0.03","l");
    leg->AddEntry(gr[1], "0.01","l");
    leg->AddEntry(gr[5], "0.01","l");
    leg->AddEntry(gr[2], "0.003","l");
    leg->AddEntry(gr[6], "0.003","l");
    leg->AddEntry(gr[3], "0.001","l");
    leg->AddEntry(gr[7], "0.001","l");

    TLegend *leg3 = new TLegend(0.53,0.33,1.18,0.38);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->SetShadowColor(0);
    leg3->SetLineColor(0);
    leg3->SetFillColor(kWhite);
    leg3->AddEntry("NULL","dz =  20 mm    dz = 150 mm","h");

    leg->Draw();
    leg3->Draw();

    TLatex * tex = new TLatex(69.91206,2.52907,"#circ");
    tex->SetTextSize(0.1);
    tex->SetLineWidth(2);
    tex->Draw(); // 120

    tex = new TLatex(13.00,9.0,"p = 5 GeV"); // 9.00 for 5 GeV, 2.7 for 120 GeV
    //tex = new TLatex(13.00,2.7,"p = 120 GeV"); // 9.00 for 5 GeV, 2.7 fo120 GeV
    tex->SetLineWidth(2);
    tex->Draw();
    c1->Update();
    c1->Modified();

    c1->Write();
    mg->Write();

    _outputFile->Close();

    std::cout << "Written to file: " << std::endl;
    std::cout << " Done." << std::endl;


  }

  if (runmode == 17)
  {

    std::cout << " read effi plots and create plots" << std::endl;

    std::vector<double> thres;
    std::vector<double> ethres;
    // Get histograms to be checked
    //
    for (Int_t j=0; j<alucount; j++)
    {
      std::cout << " j = " << j << std::endl;


      thres.push_back(aluthick[j]);
      ethres.push_back(.0);

      effi_reader( run_1_20[j], 1.0, j );
      effi_reader( run_2_20[j], 2.0, j );
      effi_reader( run_3_20[j], 3.0, j );
      effi_reader( run_4_20[j], 4.0, j );
      effi_reader( run_5_20[j], 5.0, j );

    }

    std::cout << "size of vector = " << effis.at(0).size() << std::endl;
    std::cout << "size of vector = " << effis.at(1).size()<< std::endl;
    std::cout << "size of vector = " << effis.at(9).size()<< std::endl;

    double E[5] = {2., 3., 4., 5., 6.};

    std::cout << " --------- effi read done ---------- " << std::endl;

    _outputFile->cd();
    TCanvas *cmg150 = new TCanvas("cmg150","cmg150",10,10,800,600);
    cmg150->cd();


    int marker_style[10] = {20,21,22,33,34, 24,25,26,27,28};
    int colorr[5] = {1,2,3,4,6};

    TGraphErrors *h_effi_E_150[5];
    for(int i = 0; i < 5; i++) {
      //std:: cout << i << std::endl;
      h_effi_E_150[i]= new TGraphErrors(alucount);
      h_effi_E_150[i]->SetMarkerSize(2.5);
      h_effi_E_150[i]->SetName("Graph");
      h_effi_E_150[i]->SetTitle("");
      h_effi_E_150[i]->SetMarkerStyle(marker_style[9 - i]);
      h_effi_E_150[i]->SetMarkerColor(colorr[4-i]);
    }

    for(int i = 0; i < 5; i++){
      //std:: cout << i;
      for(int j = 0; j < 10; j++){
	//std:: cout << " " << j << std::endl;
	h_effi_E_150[i]->SetPoint(j, thres.at(j)+(i-2)*0.05, (effis.at(j)).at(i));
	h_effi_E_150[i]->SetPointError(j, 0., (rms_effis.at(j)).at(i));
      }
    }

    TMultiGraph *mg150 = new TMultiGraph();
    mg150->SetName("effi150");

    for(int i = 0; i < 5; i++){
      //std:: cout << i << std::endl;
      //h_effi_E_150[i]->SetLineColor(1);
      //h_effi_E_150[i]->SetMarkerStyle(20);
      //h_effi_E_150[i]->SetMarkerSize(1.5);
      //h_effi_E_150[i]->SetMarkerColor(color[i]);
      mg150->Add(h_effi_E_150[i],"p");
    }

    cmg150->cd();
    mg150->Draw("ap");
    mg150->SetMinimum(.0);
    gPad->SetLogx();
    mg150->GetXaxis()->SetTitle("M26 threshold #xi");
    mg150->GetYaxis()->SetTitle("telescope efficiency");
    mg150->SetTitle("Mean effi vs energy");

    TLegend *leg150 = new TLegend(0.1,0.7,0.48,0.9);
    leg150->SetName("my_leg");
    leg150->SetHeader("Energy");
    for(int i = 0; i < 5; i++) {
      //std:: cout << i << std::endl;
      leg150->AddEntry(h_effi_E_150[i], std::to_string(E[i]).c_str(),"lep");   
    }

    leg150->Draw();

    std::cout << "Write to File" << std::endl;

    cmg150->Modified();
    cmg150->Update();

    mg150->Write();
    cmg150->Write();


    TCanvas *cmg20 = new TCanvas("cmg20","cmg20",10,10,1100,1100);
    cmg20->cd();

    TGraphErrors *h_effi_E_20[5];
    for(int i = 0; i < 5; i++) {
      std:: cout << i << std::endl;
      h_effi_E_20[i]= new TGraphErrors(alucount);
      h_effi_E_20[i]->SetMarkerSize(2.5);
      h_effi_E_20[i]->SetName("Graph");
      h_effi_E_20[i]->SetTitle("");
      h_effi_E_20[i]->SetMarkerStyle(marker_style[4-i]);
      h_effi_E_20[i]->SetMarkerColor(colorr[4-i]);
    }

    for(int i = 0; i < 5; i++){
      //std:: cout << i << std::endl;
      for(int j = 0; j < 10; j++){
	//std:: cout << j << std::endl;
	h_effi_E_20[i]->SetPoint(j, thres.at(j)+(i-2)*0.05, (effis.at(j)).at(i+5));
	h_effi_E_20[i]->SetPointError(j, 0., (rms_effis.at(j)).at(i));
      }
    }


    TMultiGraph *mg20 = new TMultiGraph();
    mg20->SetName("effi20");

    for(int i = 0; i < 5; i++){
      //std:: cout << i << std::endl;
      //h_effi_E_20[i]->SetLineColor(1);
      //h_effi_E_20[i]->SetMarkerStyle(25);
      //h_effi_E_20[i]->SetMarkerSize(1.5);
      //h_effi_E_20[i]->SetMarkerColor(color[i]);
      mg20->Add(h_effi_E_20[i],"p");
    }

    cmg20->cd();
    mg20->Draw("ap");
    mg20->SetMinimum(.0);
    mg20->GetXaxis()->SetTitle("M26 threshold #xi");
    mg20->GetYaxis()->SetTitle("#epsilon_{trip}");
    mg20->SetTitle("Mean effi vs energy");

    TLegend *leg20 = new TLegend(0.1,0.7,0.48,0.9);
    leg20->SetName("my_leg");
    leg20->SetHeader("Thresholds");
    leg20->SetFillColor(kWhite);
    for(int i = 0; i < 5; i++) {
      //std:: cout << i << std::endl;
      leg20->AddEntry(h_effi_E_20[i], std::to_string((int)thres.at(i)).c_str(),"lep");   
    }

    leg20->Draw();

    std::cout << "Write to File" << std::endl;

    cmg20->Modified();
    cmg20->Update();

    mg20->Write();
    cmg20->Write();



    TCanvas *cmgCombo = new TCanvas("cmgCombo","cmgcombo",20, 1920, 900,900);
    cmgCombo->Range(0,0,1,1);
    cmgCombo->SetFillColor(0);
    cmgCombo->SetBorderMode(0);
    cmgCombo->SetBorderSize(2);
    cmgCombo->SetFrameBorderMode(0);
    cmgCombo->SetRightMargin(0.05);
    cmgCombo->SetLeftMargin(0.15);
    cmgCombo->SetTickx(1);
    cmgCombo->SetTicky(1);
    cmgCombo->cd();

    TMultiGraph *mgCombo = new TMultiGraph("mgCombo","mgCombo");
    mgCombo->SetName("effiCombo");
    mgCombo->SetTitle("");

    mgCombo->Add(mg20,"p");
    mgCombo->Add(mg150,"p");

    mgCombo->Draw("ap");
    mgCombo->SetMinimum(.0);
    mgCombo->GetXaxis()->SetTitle("M26 threshold #xi_{n}");
    mgCombo->GetYaxis()->SetTitle("triplet efficiency #epsilon_{trip}");
    mgCombo->GetXaxis()->SetLimits(2,13);
    mgCombo->SetMinimum(0.8);
    mgCombo->SetMaximum(1.0);
    mgCombo->GetXaxis()->SetTitleOffset(1.14);
    mgCombo->GetYaxis()->SetTitleOffset(1.74);

    TLegend *legCombo = new TLegend(0.22,0.17,0.52,0.52,NULL,"brNDC");
    legCombo->SetBorderSize(0);
    legCombo->SetTextSize(0.03);
    legCombo->SetFillStyle(1001);
    legCombo->SetShadowColor(0);
    legCombo->SetLineColor(0);
    legCombo->SetFillColor(kWhite);
    //for(int i = 0; i < 5; i++)
    //  legCombo->AddEntry(h_effi_E_150[i], std::to_string((int)thres.at(i)).c_str(),"lep");   

    legCombo->AddEntry("NULL","dz = 20 mm","h");
    legCombo->AddEntry(h_effi_E_20[4],"p =  6.0 GeV","p");
    legCombo->AddEntry(h_effi_E_20[3],"p =  5.0 GeV","p");
    legCombo->AddEntry(h_effi_E_20[2],"p =  4.0 GeV","p");
    legCombo->AddEntry(h_effi_E_20[1],"p =  3.0 GeV","p");
    legCombo->AddEntry(h_effi_E_20[0],"p =  2.0 GeV","p");
    legCombo->AddEntry("NULL","","");
    legCombo->AddEntry("NULL","dz = 150 mm","h");
    legCombo->AddEntry(h_effi_E_150[4],"p =  6.0 GeV","p");
    legCombo->AddEntry(h_effi_E_150[3],"p =  5.0 GeV","p");
    legCombo->AddEntry(h_effi_E_150[2],"p =  4.0 GeV","p");
    legCombo->AddEntry(h_effi_E_150[1],"p =  3.0 GeV","p");
    legCombo->AddEntry(h_effi_E_150[0],"p =  2.0 GeV","p");
    legCombo->Draw();

    cmgCombo->Modified();
    cmgCombo->Update();

    mgCombo->Write();
    cmgCombo->Write();

    _outputFile->Close();

  }

  if (runmode == 16)
  {
    std::cout << " Plot predicted and measured residual width versus z position (smiley) " << std::endl;

    // need various graphs
    // 1) with points for X
    // 2) green band for uncertainty
    // 3) predictions for sigma int = 2, 3, 4, 5 um (need all?)

    TMultiGraph *mg = new TMultiGraph();

    bool Is20 = 0;
    // First fill measurements 
    //double z[6] = {0., 20., 40., 60., 80., 100.};
    //double resXY[6]  = {2.17179, 2.62271, 2.87623, 2.87197, 2.63783, 2.17014};
    //double eresXY[6] = {0.};
    //for(int i = 0; i< 6; i++) eresXY[i] = resXY[i]*0.024;

    double z[6] = {0., 150., 300., 450., 600., 750.}; // 150
    double resXY[6]  = {1.09139, 2.04053, 2.06699, 2.02072, 2.00273, 1.07391}; // 150
    double eresXY[6] = {0.}; 
    for(int i = 0; i< 6; i++) eresXY[i] = resXY[i]*0.079;


    // mean of X and Y per Z position
    TGraphErrors *res_meas_XY = new TGraphErrors(6);
    res_meas_XY->SetName("");
    res_meas_XY->SetTitle("");
    res_meas_XY->SetLineWidth(1);
    res_meas_XY->SetLineColor(kGreen);
    res_meas_XY->SetMarkerStyle(22);
    res_meas_XY->SetMarkerSize(2.5);
    res_meas_XY->SetFillColor(kGreen);    

    for (int i = 0; i< 6; i++){
      res_meas_XY->SetPoint(i, z[i] ,  resXY[i]);
      res_meas_XY->SetPointError(i,0, eresXY[i]);
    }

    mg->Add(res_meas_XY,"3");

    TGraph *res_meas_XY2 = new TGraph(6);
    res_meas_XY2->SetName("");
    res_meas_XY2->SetTitle("");
    res_meas_XY2->SetLineWidth(1);
    res_meas_XY2->SetLineColor(1);
    res_meas_XY2->SetMarkerColor(2);
    res_meas_XY2->SetMarkerStyle(22);
    res_meas_XY2->SetMarkerSize(2.5);
    res_meas_XY2->SetFillColor(kGreen);    

    for (int i = 0; i< 6; i++){
      res_meas_XY2->SetPoint(i, z[i] ,  resXY[i]);
    }

    mg->Add(res_meas_XY2,"p");

    // Second add estimates for other values of sigma int
    double ebeam = 6.;

    // Create telescope
    //
    telescopebuild = "AnaTel_narrow.geom"; // 20
    if(!Is20) telescopebuild = "AnaTel_wide.geom"; // 150
    AnaTel *tscope20 = new AnaTel(telescopebuild.c_str(), ebeam);
    tscope20->SetKappa(0.86);
    std::cout << " - Tscope created, NO DUT" << std::endl;

    // vary z
    const int steps_z = 6;
    double step_width_z = 20; //   20
    if(!Is20) step_width_z = 150; // 150
    double start_z = 0;
    std::vector<double> zs;

    for(unsigned int ii = 0; ii<steps_z; ii++) zs.push_back(start_z + ii*step_width_z);

    double sigma_int = 2.5e-3;
    tscope20->SetResolution(sigma_int);
    TGraph* gr_rhat_2 = new TGraph(6);
    gr_rhat_2->SetNameTitle("","");
    double rhat = -1.;
    for (int ii = 0; ii < steps_z; ii++){
      rhat = tscope20->GetRhat(ii,1);
      //std::cout << " rhat is " << rhat << std::endl;
      gr_rhat_2->SetPoint(ii, zs.at(ii), rhat*1e3);
      rhat = -1.;
    }
    gr_rhat_2->SetLineWidth(2);
    gr_rhat_2->SetLineStyle(2);
    mg->Add(gr_rhat_2,"l");

    tscope20->SetResolution(3.0e-3);
    TGraph* gr_rhat_3 = new TGraph(6);
    gr_rhat_3->SetNameTitle("","");
    for (int ii = 0; ii < steps_z; ii++){
      rhat = tscope20->GetRhat(ii,1);
      gr_rhat_3->SetPoint(ii, zs.at(ii), rhat*1e3);
      rhat = -1.;
    }
    gr_rhat_3->SetLineWidth(2);
    gr_rhat_3->SetLineStyle(3);
    mg->Add(gr_rhat_3,"l");

    tscope20->SetResolution(3.5e-3);
    TGraph* gr_rhat_4 = new TGraph(6);
    gr_rhat_4->SetNameTitle("","");
    for (int ii = 0; ii < steps_z; ii++){
      rhat = tscope20->GetRhat(ii,1);
      gr_rhat_4->SetPoint(ii, zs.at(ii), rhat*1e3);
      rhat = -1.;
    }
    gr_rhat_4->SetLineWidth(2);
    gr_rhat_4->SetLineStyle(4);
    mg->Add(gr_rhat_4,"l");

    tscope20->SetResolution(4.0e-3);
    TGraph* gr_rhat_5 = new TGraph(6);
    gr_rhat_5->SetNameTitle("","");
    for (int ii = 0; ii < steps_z; ii++){
      rhat = tscope20->GetRhat(ii,1);
      gr_rhat_5->SetPoint(ii, zs.at(ii), rhat*1e3);
      rhat = -1.;
    }
    gr_rhat_5->SetLineWidth(2);
    gr_rhat_5->SetLineStyle(5);
    mg->Add(gr_rhat_5,"l");

    tscope20->SetResolution(3.236e-3);
    TGraph* gr_rhat_true = new TGraph(6);
    gr_rhat_true->SetNameTitle("","");
    for (int ii = 0; ii < steps_z; ii++){
      rhat = tscope20->GetRhat(ii,1);
      gr_rhat_true->SetPoint(ii, zs.at(ii), rhat*1e3);
      rhat = -1.;
    }
    gr_rhat_true->SetLineWidth(2);
    mg->Add(gr_rhat_true,"l");

    TCanvas *smiley = new TCanvas("smiley","smiley",20,1920,800,800);
    smiley->Range(0,0,1,1);
    smiley->SetFillColor(0);
    smiley->SetBorderMode(0);
    smiley->SetBorderSize(2);
    smiley->SetFrameBorderMode(0);
    smiley->SetRightMargin(0.05);
    smiley->SetLeftMargin(0.15);
    smiley->SetTickx(1);
    smiley->SetTicky(1);
    smiley->cd();

    mg->Draw("apl");

    mg->SetMinimum(0);
    mg->SetMaximum(8);
    mg->GetXaxis()->SetLimits(-10,110); // 20
    if(!Is20) mg->GetXaxis()->SetLimits(-75,825); // 150
    mg->SetTitle("");
    mg->GetXaxis()->SetTitle("z [mm]");
    mg->GetYaxis()->SetTitle("biased residual width [#mum]");
    mg->GetXaxis()->SetTitleOffset(1.14);
    mg->GetYaxis()->SetTitleOffset(1.14);

    TLegend *leg = new TLegend(.22,0.55,0.47,0.85,NULL,"brNDC");
    leg->SetBorderSize(0);
    leg->SetTextSize(0.03);
    leg->SetFillStyle(4000);
    leg->SetShadowColor(0);
    leg->SetLineColor(0);
    leg->SetFillColor(kWhite);
    TLegendEntry *entry=
      leg->AddEntry("NULL","Predictions:","h");
    leg->AddEntry(gr_rhat_5,"#sigma_{M26} = 4.0 #mum","l");
    leg->AddEntry(gr_rhat_4,"#sigma_{M26} = 3.5 #mum","l");
    leg->AddEntry(gr_rhat_true,"#sigma_{M26} = 3.24 #mum","l");
    leg->AddEntry(gr_rhat_3,"#sigma_{M26} = 3.0 #mum","l");
    leg->AddEntry(gr_rhat_2,"#sigma_{M26} = 2.5 #mum","l");
    //entry->SetLineWidth(3);

    leg->Draw();

    TLegend *leg2 = new TLegend(0.57,0.55,0.82,0.85,NULL,"brNDC");
    leg2->SetBorderSize(0);
    leg2->SetTextSize(0.03);
    leg2->SetLineColor(1);
    leg2->SetLineStyle(1);
    leg2->SetLineWidth(1);
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry("NULL","DESY 6 GeV e^{-}","");
    leg2->AddEntry("NULL","M26 threshold: 6","");
    if(Is20)leg2->AddEntry("NULL","dz = 20 mm","");
    else leg2->AddEntry("NULL","dz = 150 mm","");
    leg2->AddEntry(res_meas_XY2,"meas. biased residual","p");
    leg2->AddEntry(res_meas_XY, "1 #sigma uncertainty","f");
    leg2->Draw();

    smiley->Write();

  }



  if (runmode == 12)
  {
    std::cout << " This checks the mean number of fired pixel in a cluster (a.k.a. cluster charge in case of digital chip)." <<
      "\n - Attention! You can choose if average over ALL clusters or just those that make it into a track! Default is those from tracks" <<  std::endl;

    whichfitter = "GBLfitter-CS-CSdepSigmaInt";


    std::vector<double> thres;
    std::vector<double> ethres;
    // Get histograms to be checked
    //
    for (Int_t j=0; j<alucount; j++)
    {
      std::cout << " j = " << j << std::endl;


      thres.push_back(aluthick[j]);
      ethres.push_back(.0);

      CSchecker( run_1_20[j], 1.0, j );
      CSchecker( run_2_20[j], 2.0, j );
      CSchecker( run_3_20[j], 3.0, j );
      CSchecker( run_4_20[j], 4.0, j );
      CSchecker( run_5_20[j], 5.0, j );

    }

    _outputFile->cd();
    TCanvas *cmg150 = new TCanvas("cmg150","cmg150",10,10,800,600);
    cmg150->cd();

    TGraphErrors *h_clustersize_E_150[10];
    for(int i = 0; i < 10; i++)
      h_clustersize_E_150[i]= new TGraphErrors(5);

    double E[5] = {2., 3., 4., 5., 6.};
    double eE[5];
    for(int i = 0; i < 5; i++) 
      eE[i] = 0.05*E[i];

    for(int i = 0; i < 10; i++){
      for(int j = 0; j < 5; j++){
	h_clustersize_E_150[i]->SetPoint(j, E[j], (clustersizes.at(i)).at(j));
	h_clustersize_E_150[i]->SetPointError(j, eE[j],(rms_clustersizes.at(i)).at(j));
      }
    }

    TMultiGraph *mg150 = new TMultiGraph();
    mg150->SetName("Clustersizes150");

    for(int i = 0; i < 10; i++){
      h_clustersize_E_150[i]->SetLineColor(1);
      h_clustersize_E_150[i]->SetMarkerStyle(20);
      h_clustersize_E_150[i]->SetMarkerSize(1.5);
      h_clustersize_E_150[i]->SetMarkerColor(color[i]);
      mg150->Add(h_clustersize_E_150[i],"p");
    }

    cmg150->cd();
    mg150->Draw("ap");
    mg150->SetMinimum(.0);
    //gPad->SetLogx();
    mg150->GetXaxis()->SetTitle("Energy [GeV]");
    mg150->GetYaxis()->SetTitle("<CS>");
    mg150->SetTitle("Mean clustersize vs energy");

    TLegend *leg150 = new TLegend(0.1,0.7,0.48,0.9);
    leg150->SetName("my_leg");
    leg150->SetFillColor(0);
    leg150->SetFillStyle(0);
    leg150->SetHeader("Thresholds");
    for(int i = 0; i < 10; i++)
      leg150->AddEntry(h_clustersize_E_150[i], std::to_string((int)thres.at(i)).c_str(),"lep");   

    leg150->Draw();

    std::cout << "Write to File" << std::endl;

    cmg150->Modified();
    cmg150->Update();

    mg150->Write();
    cmg150->Write();

    // now a.f.o threshold
    TCanvas *cmg150_thr = new TCanvas("cmg150_thr","cmg150_thr",10,10,800,600);
    cmg150_thr->cd();

    TGraphErrors *h_clustersize_E_150_thr[5];
    for(int i = 0; i < 5; i++)
      h_clustersize_E_150_thr[i]= new TGraphErrors(10);

    for(int i = 0; i < 5; i++){
      for(int j = 0; j < 10; j++){
	h_clustersize_E_150_thr[i]->SetPoint(j, j+3, (clustersizes.at(j)).at(i));
	h_clustersize_E_150_thr[i]->SetPointError(j, 0,(rms_clustersizes.at(j)).at(i));
      }
    }

    TMultiGraph *mg150_thr = new TMultiGraph();
    mg150_thr->SetName("Clustersizes150_thr");

    for(int i = 0; i < 5; i++){
      h_clustersize_E_150_thr[i]->SetLineColor(1);
      h_clustersize_E_150_thr[i]->SetMarkerStyle(20);
      h_clustersize_E_150_thr[i]->SetMarkerSize(1.5);
      h_clustersize_E_150_thr[i]->SetMarkerColor(color[i]);
      mg150_thr->Add(h_clustersize_E_150_thr[i],"p");
    }

    cmg150_thr->cd();
    mg150_thr->Draw("ap");
    mg150_thr->SetMinimum(.0);
    //gPad->SetLogx();
    mg150_thr->GetXaxis()->SetTitle("M26 threshold #xi_{n}");
    mg150_thr->GetYaxis()->SetTitle("<CS>");
    mg150_thr->SetTitle("Mean clustersize vs threshold");

    TLegend *leg150_thr = new TLegend(0.5,0.7,0.88,0.9);
    leg150_thr->SetName("my_leg");
    leg150_thr->SetFillColor(0);
    leg150_thr->SetFillStyle(0);
    leg150_thr->SetHeader("Beam momenta");
    for(int i = 0; i < 5; i++)
      leg150_thr->AddEntry(h_clustersize_E_150_thr[i], std::to_string((int)E[i]).c_str(),"lep");   

    leg150_thr->Draw();

    std::cout << "Write to File" << std::endl;

    cmg150_thr->Modified();
    cmg150_thr->Update();

    mg150_thr->Write();
    cmg150_thr->Write();

    // now for 20 mm

    TCanvas *cmg20 = new TCanvas("cmg20","cmg20",10,10,800,600);
    cmg20->cd();

    TGraphErrors *h_clustersize_E_20[10];
    for(int i = 0; i < 10; i++)
      h_clustersize_E_20[i]= new TGraphErrors(5);

    for(int i = 0; i < 10; i++){
      for(int j = 0; j < 5; j++){
	h_clustersize_E_20[i]->SetPoint(j, E[j], (clustersizes.at(i)).at(j+5));
	h_clustersize_E_20[i]->SetPointError(j, eE[j], (rms_clustersizes.at(i)).at(j+5));
      }
    }

    TMultiGraph *mg20 = new TMultiGraph();
    mg20->SetName("Clustersizes20");

    for(int i = 0; i < 10; i++){
      h_clustersize_E_20[i]->SetLineColor(1);
      h_clustersize_E_20[i]->SetMarkerStyle(25);
      h_clustersize_E_20[i]->SetMarkerSize(1.5);
      h_clustersize_E_20[i]->SetMarkerColor(color[i]);
      mg20->Add(h_clustersize_E_20[i],"p");
    }

    cmg20->cd();
    mg20->Draw("ap");
    mg20->SetMinimum(.0);
    //gPad->SetLogx();
    mg20->GetXaxis()->SetTitle("Energy [GeV]");
    mg20->GetYaxis()->SetTitle("<CS>");
    mg20->SetTitle("Mean clustersize vs energy");

    TLegend *leg20 = new TLegend(0.1,0.7,0.48,0.9);
    leg20->SetName("my_leg");
    leg20->SetFillColor(0);
    leg20->SetFillStyle(0);
    leg20->SetHeader("Thresholds");
    leg20->SetFillColor(kWhite);
    for(int i = 0; i < 10; i++)
      leg20->AddEntry(h_clustersize_E_20[i], std::to_string((int)thres.at(i)).c_str(),"lep");   

    leg20->Draw();

    std::cout << "Write to File" << std::endl;

    cmg20->Modified();
    cmg20->Update();

    mg20->Write();
    cmg20->Write();


    // now a.f.o threshold
    TCanvas *cmg20_thr = new TCanvas("cmg20_thr","cmg20_thr",10,10,800,600);
    cmg20_thr->cd();

    TGraphErrors *h_clustersize_E_20_thr[5];
    for(int i = 0; i < 5; i++)
      h_clustersize_E_20_thr[i]= new TGraphErrors(10);

    for(int i = 0; i < 5; i++){
      for(int j = 0; j < 10; j++){
	h_clustersize_E_20_thr[i]->SetPoint(j, j+3, (clustersizes.at(j)).at(i+5));
	h_clustersize_E_20_thr[i]->SetPointError(j, 0,(rms_clustersizes.at(j)).at(i+5));
      }
    }

    TMultiGraph *mg20_thr = new TMultiGraph();
    mg20_thr->SetName("Clustersizes20_thr");

    for(int i = 0; i < 5; i++){
      h_clustersize_E_20_thr[i]->SetLineColor(1);
      h_clustersize_E_20_thr[i]->SetMarkerStyle(20);
      h_clustersize_E_20_thr[i]->SetMarkerSize(1.5);
      h_clustersize_E_20_thr[i]->SetMarkerColor(color[i]);
      mg20_thr->Add(h_clustersize_E_20_thr[i],"p");
    }

    cmg20_thr->cd();
    mg20_thr->Draw("ap");
    mg20_thr->SetMinimum(.0);
    //gPad->SetLogx();
    mg20_thr->GetXaxis()->SetTitle("M26 threshold #xi_{n}");
    mg20_thr->GetYaxis()->SetTitle("<CS>");
    mg20_thr->SetTitle("Mean clustersize vs threshold");

    TLegend *leg20_thr = new TLegend(0.5,0.7,0.88,0.9);
    leg20_thr->SetName("my_leg");
    leg20_thr->SetFillColor(0);
    leg20_thr->SetFillStyle(0);
    leg20_thr->SetHeader("Beam momenta");
    for(int i = 0; i < 5; i++)
      leg20_thr->AddEntry(h_clustersize_E_20_thr[i], std::to_string((int)E[i]).c_str(),"lep");   

    leg20_thr->Draw();

    std::cout << "Write to File" << std::endl;

    cmg20_thr->Modified();
    cmg20_thr->Update();

    mg20_thr->Write();
    cmg20_thr->Write();



    TCanvas *cmgCombo = new TCanvas("cmgCombo","cmgcombo",10,10,800,600);
    cmgCombo->cd();

    TMultiGraph *mgCombo = new TMultiGraph();
    mgCombo->SetName("ClustersizesCombo");

    mgCombo->Add(mg20,"");
    mgCombo->Add(mg150,"");

    mgCombo->Draw("ap");
    mgCombo->SetMinimum(.0);
    mgCombo->GetXaxis()->SetTitle("Energy [GeV]");
    mgCombo->GetYaxis()->SetTitle("<CS>");

    TLegend *legCombo = new TLegend(0.5,0.5,0.88,0.88);
    legCombo->SetHeader("Thresholds");
    legCombo->SetFillColor(kWhite);
    for(int i = 0; i < 10; i++)
      legCombo->AddEntry(h_clustersize_E_150[i], std::to_string((int)thres.at(i)).c_str(),"lep");   

    legCombo->AddEntry((TObject*)0, " FULL symbols 150, open symbols 20 mm","");   
    legCombo->Draw();

    cmgCombo->Modified();
    cmgCombo->Update();

    mgCombo->Write();
    cmgCombo->Write();

    _outputFile->Close();
  }

  if (runmode == 111)
  {

    std::cout << " read pulls and write to file" << std::endl;

    whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT3";
    //whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT3_Tp35";
    IsBiased = true;
    if(!IsBiased) whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT3_unb";
    //if(!IsBiased) whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT3_unb_Ep5";

    //pulls20.open("pulls20-xx.txt");
    //pulls150.open("pulls150-xx.txt");
    if (IsBiased){
      pulls20.open("pulls20-09.txt");
      pulls150.open("pulls150-09.txt");
    } else {
      pulls20.open("pulls20-26_unb.txt");
      pulls150.open("pulls150-26_unb.txt");
    }

    telescopebuild = "AnaTel_narrow.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    fitter(117,6.0);
    //fitter(127,5.0);
    //fitter(137,4.0);
    //fitter(153,3.0);
    //fitter(165,2.0);
    //for( int i : run_6_20) fitter(i,6.0);

    telescopebuild = "AnaTel_wide.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    fitter(63,6.0);
    //fitter(73,5.0);
    //fitter(84,4.0);
    //fitter(96,3.0);
    fitter(106,2.0);


    //std::cout << "fill tgraphs" << std::endl;
    //FillGraph();

    pulls20.close();
    pulls150.close();

  }
  if (runmode == 11)
  {

    std::cout << " read pulls and write to file" << std::endl;
    whichfitter = "GBLfitter-CS";

    //pulls20.open("pulls20-xx.txt");
    //pulls20.open("pulls20-12.txt");

    telescopebuild = "AnaTel_narrow.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    for( int i : run_5_20) fitter(i,5.0);
    for( int i : run_4_20) fitter(i,4.0);
    for( int i : run_3_20) fitter(i,3.0);
    for( int i : run_2_20) fitter(i,2.0);
    for( int i : run_1_20) fitter(i,1.0);


 

    //std::cout << "fill tgraphs" << std::endl;
    FillGraph();

    //pulls20.close();

  }

  // test mode
  if (runmode == 9)
  {
    global_thresh=6;



    telescopebuild = "AnaTel_narrow.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    fitter(117,6.);
    //fitter(127,5.);
    //fitter(137,4.);
    //fitter(153,3.);
    //fitter(165,2.);

    std::cout << "fill tgraphs" << std::endl;
    FillGraph();
  }





  // Simply run over everything
  if(runmode == 0)
  {

    std::string whichfitter = "GBLKinkEstimator_kappa100";

    std::cout << " " << std::endl;
    std::cout << "Thin Geometry" << std::endl;
    std::cout << " " << std::endl;

    telescopebuild = "AnaTel_narrow.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = 20.0*j;

    for(int i=0;i<run_5_20.size();i++)
    {
      fitter( run_5_20[i], 5.0 );
      fitter( run_4_20[i], 4.0 );
      fitter( run_3_20[i], 3.0 );
      fitter( run_2_20[i], 2.0 );
      fitter( run_1_20[i], 1.0 );
    }


    std::cout << "fill tgraphs" << std::endl;
    FillGraph();

  }


  // Run over alu thicknesses
  if(runmode == 2)
  {
    std::cout << "Thickness mode" << std::endl;
    Double_t threshresult_1_20[alucount];
    Double_t threshresult_2_20[alucount];
    Double_t threshresult_3_20[alucount];
    Double_t threshresult_4_20[alucount];
    Double_t threshresult_5_20[alucount];
    Double_t threshresult120[alucount] = {0.};


    Double_t thresherror_1_20[alucount];
    Double_t thresherror_2_20[alucount];
    Double_t thresherror_3_20[alucount];
    Double_t thresherror_4_20[alucount];
    Double_t thresherror_5_20[alucount];


    Double_t x[alucount];
    Double_t xerrorthresh[alucount];
    for (Int_t j=0; j<alucount; j++)
    {
      std::cout << " j = " << j << std::endl;

 
      telescopebuild = "AnaTel_narrow.geom";
      planedistance = 20;
      for(int jj=0;jj<nplanes;jj++)
	posx[jj] = planedistance*jj;


      fitter( run_1_20[j], 1.0 );
      threshresult_1_20[j] = m26_resolution*1000.0;
      thresherror_1_20[j] = global_plot_error;

      fitter( run_2_20[j], 2.0 );
      threshresult_2_20[j] = m26_resolution*1000.0;
      thresherror_2_20[j] = global_plot_error;

      fitter( run_3_20[j], 3.0 );
      threshresult_3_20[j] = m26_resolution*1000.0;
      thresherror_3_20[j] = global_plot_error;

      fitter( run_4_20[j], 4.0 );
      threshresult_4_20[j] = m26_resolution*1000.0;
      thresherror_4_20[j] = global_plot_error;

      fitter( run_5_20[j], 5.0 );
      threshresult_5_20[j] = m26_resolution*1000.0;
      thresherror_5_20[j] = global_plot_error;


    }

    std::cout << "fill tgraphs" << std::endl;
    FillGraph();

    TGraphErrors *gr1_20  = new TGraphErrors((alucount),x,threshresult_1_20,xerrorthresh,thresherror_1_20);
    TGraphErrors *gr2_20  = new TGraphErrors((alucount),x,threshresult_2_20,xerrorthresh,thresherror_2_20);
    TGraphErrors *gr3_20  = new TGraphErrors((alucount),x,threshresult_3_20,xerrorthresh,thresherror_3_20);
    TGraphErrors *gr4_20  = new TGraphErrors((alucount),x,threshresult_4_20,xerrorthresh,thresherror_4_20);
    TGraphErrors *gr5_20  = new TGraphErrors((alucount),x,threshresult_5_20,xerrorthresh,thresherror_5_20);

    TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
    threshold->cd();


    //TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
    //gStyle->SetPadBorderMode(0);
    //gStyle->SetOptStat(0);
    threshold->SetFillColor(0);
    //gPad->SetGridx();
    //gPad->SetGridy();

    gr1_20->SetMarkerStyle(28);
    gr1_20->SetMarkerColor(kOrange);
    gr1_20->SetMarkerSize(3);
    gr2_20->SetMarkerStyle(24);
    gr2_20->SetMarkerColor(kBlack);
    gr2_20->SetMarkerSize(3);
    gr3_20->SetMarkerStyle(25);
    gr3_20->SetMarkerColor(kGreen);
    gr3_20->SetMarkerSize(3);
    gr4_20->SetMarkerStyle(26);
    gr4_20->SetMarkerColor(kGreen);
    gr4_20->SetMarkerSize(3);
    gr5_20->SetMarkerStyle(27);
    gr5_20->SetMarkerColor(kBlue);
    gr5_20->SetMarkerSize(3);

    TMultiGraph* mg = new TMultiGraph();

    mg->Add(gr1_20,"p");
    mg->Add(gr2_20,"p");
    mg->Add(gr3_20,"p");
    mg->Add(gr4_20,"p");
    mg->Add(gr5_20,"p");

    mg->Draw("apl");
    mg->SetMinimum(0);
    mg->GetXaxis()->SetLimits(2,15);
    mg->GetXaxis()->SetTitle("Threshold #xi");
    mg->GetYaxis()->SetTitle("#sigma_{M26} [#mum]");


    TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Energy:");

    leg->AddEntry(gr1_20,"p = 1 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr2_20,"p = 2 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr3_20,"p = 3 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr4_20,"p = 4 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr5_20,"p = 5 GeV, #Delta_{z} = 20 mm","p");

    leg->Draw();

    // Output
    threshold->Print("pics/threshold_all.eps");
    _outputFile->cd();
    threshold->Write();
    threshold->Close();
  }

  _outputFile->Close();

  // And we're done
  std::cout << "\nDone.\n" << std::endl;

}
