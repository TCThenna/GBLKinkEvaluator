/******************************************************************************************
 * Original code by Joerg Behr & Igor Rubinskiy                                           *
 * 2012                                                                                   *
 *           - functions to plot 6 Mimosa26 planes residuals                              *
 *             and extract plane intrinsic resolution                                     *
 *             along with telescope pointing resolution                                   *
 *                                                                                        *
 *                                                                                        *
 ******************************************************************************************/

// rewrite as of 03.07.2013: Thomas Eichhorn
// rewrite as of 2015: Hendrik Jansen Email: <hendrik.jansen@desy.de>
// The interface to AnaTel.c is the same, but AnaTel now uses GBL for pointing res prediction
// 	- usually unbiased residuals were used
// 	- Oct15: starting with biased residuals, as suggsted by Claus Kleinwort
// 	- Dez15: Also getting pulls, which makes minimisation for calculation of sigma int obsolete
// 		Rather iterate the tracking varying the input sigma int until width<pull> = 1
// 		Seems that Highland needs correction of ~ -25% to bring 20 mm and 150mm results to an agreement
// 		correction seems to be stable over energy for biased track fits


//#include "configure_histos.h"
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

std::string inputDir20;
std::string inputDir150;
std::string inputFile;
TFile* _outputFile;
ofstream pulls20;
ofstream pulls150;

//std::string whichfitter = "fitter";
//std::string whichfitter = "daffitter";
//std::string whichfitter = "daffitter_fixedAlign";
//std::string whichfitter = "daffitter-alignGBL";
//std::string whichfitter = "daffitter_fixedAlign-alignGBL";
//std::string whichfitter = "GBLfitter-alignGBL";
//std::string whichfitter = "GBLfitter_fixedAlign-alignGBL";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir";
//std::string whichfitter = "GBLfitter-alignGBL_wAir";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSens50";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSens40";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSens";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSens70";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT0";
//std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT1";
std::string whichfitter = "GBLfitter-fixedAlign-alignGBL_wAir_thickSensIT3";

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


const Int_t threshcount = 10;
const Int_t planescount = 12;

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
Bool_t verbose = true;

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

// Cut in clustersize, set to 0 for all
#define clusterlimit 0

// What subsensor information do we want? Format: A-D for submatrix, 1-7 for clustersize, empty for all.
// e.g.: submask = "A1" or submask = "3" or submask = "C"
TString submask="";

// 2 Dimensional chi2 (X, Y) with uncertainty in the fit for beam, spread and sensor thickness
class  MyFunctionObject2D_beam {
  public:
    MyFunctionObject2D_beam(AnaTel *t, Double_t *arr1, Double_t *err1, Double_t *arr2, Double_t *err2)
    {
      tl = t;
      measured1 = arr1;
      error1    = err1;
      measured2 = arr2;
      error2    = err2;
    }
    double getchi2_2D(Double_t *par, Double_t *p)
    {
      // std::cout << "chi2call" << std::endl;
      //std::cout << " Current intrinsic resolution set from gMinuit: " << par[0] << std::endl;
      tl->SetResolution(par[0]);
      //tl->SetBeam(par[1], par[2]);
      tl->SetBeam(par[1], 0.0);
      tl->SetThickness(par[3]);
      Double_t chi2 =0.0;
      Double_t chi2GBL =0.0;
      float w = 0.;
      float wGBL = 0.;
      for(Int_t j = 0; j < nplanes; j++)
      {
	if(!IsBiased)
	{
	  w    = tl->GetWidth(j,0);
	  wGBL = tl->GetWidthGBL(j,0);
	}else{
	  w    = tl->GetWidth(j,1);
	  wGBL = tl->GetWidthGBL(j,1);
	}


	chi2    += pow( (w    - measured1[j])/(error1[j]) ,2) + pow( (w    - measured2[j])/(error2[j]) ,2) ; 
	chi2GBL += pow( (wGBL - measured1[j])/(error1[j]) ,2) + pow( (wGBL - measured2[j])/(error2[j]) ,2) ; 

      }
      //return chi2;
      return chi2GBL;

    }

    double getchi2_2D_nomean(Double_t *par, Double_t *p)
    {
      //cout << "chi2call" << std::endl;
      std::cout << " Current intrinsic resolutions set from gMinuit: " << par[0] << " -- " << par[1] << " -- " << par[2] << " -- " << par[3] << " -- " << par[4] << " -- " << par[5] << std::endl;
      double res[6] = {par[0], par[1], par[2], par[3], par[4], par[5]};
      if(std::isnan(par[0])){
	for(unsigned int i = 0; i<6; i++) res[i] = 4.e-3;
      }
      tl->SetResolution(res);
      //tl->SetBeam(par[6], par[7]);
      tl->SetBeam(par[6], 0.0);
      tl->SetThickness(par[8]);
      Double_t chi2 =0.0;
      Double_t chi2GBL =0.0;
      float w = 0.;
      float wGBL = 0.;
      for(Int_t j = 0; j < nplanes; j++)
      {
	if(!IsBiased)
	{
	  w    = tl->GetWidth(j,0);
	  wGBL = tl->GetWidthGBL(j,0);
	}else{
	  w    = tl->GetWidth(j,1);
	  wGBL = tl->GetWidthGBL(j,1);
	}


	chi2    += pow( (w    - measured1[j])/(error1[j]) ,2) + pow( (w    - measured2[j])/(error2[j]) ,2) ; 
	chi2GBL += pow( (wGBL - measured1[j])/(error1[j]) ,2) + pow( (wGBL - measured2[j])/(error2[j]) ,2) ; 

      }
      //return chi2;
      return chi2GBL;

    }

    double getchi2_1D(Double_t *par, Double_t *p)
    {
      tl->SetResolution(par[0]);
      //tl->SetBeam(par[1], par[2]);
      tl->SetBeam(par[1], 0.0);
      tl->SetThickness(par[3]);
      Double_t chi2 =0.0;
      for(Int_t j = 0; j < nplanes; j++ )
      {
	chi2 +=  pow( (tl->GetWidth(j,0)*1000.0 - measured1[j])/(error1[j]) ,2) ; 
      }
      return chi2;
    }

    AnaTel *tl;
    Double_t *measured1;
    Double_t *error1   ;
    Double_t *measured2;
    Double_t *error2   ;
};

/////////////////////////////////////////////////////////////////////////////////////////////
MyFunctionObject2D_beam* fcn_beam=0;

// use MINUIT: 
void fcn_wrapper(int &npar, double *gin, double &f, double *par, int iflag)
{

  double p[4];
  p[0] = 0.0;
  p[1] = 0.0;
  p[2] = 0.0;
  p[3] = 0.0;
  f = fcn_beam->getchi2_2D(par,p);
}
//////////////////////////////////////////////////////////////////////////////////////////////

// use MINUIT: 
void fcn_wrapper_nomean(int &npar, double *gin, double &f, double *par, int iflag)
{

  double p[9];
  p[0] = 0.0;
  p[1] = 0.0;
  p[2] = 0.0;
  p[3] = 0.0;
  p[4] = 0.0;
  p[5] = 0.0;
  p[6] = 0.0;
  p[7] = 0.0;
  p[8] = 0.0;
  f = fcn_beam->getchi2_2D_nomean(par,p);
}
//////////////////////////////////////////////////////////////////////////////////////////////


// 1 Dimensional chi2 (X or Y)
class  MyFunctionObject1D {
  public:
    MyFunctionObject1D(AnaTel *t, Double_t *arr1, Double_t *err1)
    {
      tl = t;
      measured1 = arr1;
      error1    = err1;
    }
    double getchi2_1D(double resolution)
    {
      tl->SetResolution(resolution);
      Double_t chi2 = 0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 += pow( (tl->GetWidth(j,0)*1000.0 - measured1[j])/(error1[j]) ,2);
      }
      return chi2;
    }
    double operator() (double *x, double *p)
    {
      Double_t chi2 = 0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 = getchi2_1D(x[0]) + p[0];
      }
      return chi2;
    }
    AnaTel *tl;
    Double_t *measured1;
    Double_t *error1;
};


// 2 Dimensional chi2 (X,Y)
class  MyFunctionObject2D {
  public:
    MyFunctionObject2D(AnaTel *t, Double_t *arr1, Double_t *arr2, Double_t *err1, Double_t *err2)
    {
      tl = t;
      measured1 = arr1;
      measured2 = arr2;
      error1 = err1;
      error2 = err2;
    }
    double getchi2_2D(double resolution)
    {



      tl->SetResolution(resolution);
      Double_t chi2 =0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 += pow((tl->GetWidth(j,0)*1000.0 - measured1[j])/(error1[j]) ,2) + pow( (tl->GetWidth(j,0)*1000.0 - measured2[j])/(error2[j]) ,2);

      }
      return chi2;
    }
    double operator() (double *x, double *p)
    {
      Double_t chi2 = 0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 = getchi2_2D(x[0]) + p[0];
      }
      return chi2;
    }
    AnaTel *tl;
    Double_t *measured1;
    Double_t *measured2;
    Double_t *error1;
    Double_t *error2;
};


// Get the telescope resolution from the measured input via k-factor
Double_t resol_estimate(const Double_t measured, const Int_t pos)
{
  std::vector<Double_t> vec_z;
  for(int i=0; i < nplanes; i++)
  {
    vec_z.push_back(posx[i]);
  }

  const Double_t offset = vec_z[pos];
  const Double_t n = (Double_t)vec_z.size();
  Double_t s2sum = 0.0;
  Double_t sum = 0.0;
  for(Int_t i =0; i < nplanes;i++)
  {
    if(i != pos)
    {
      s2sum += pow((vec_z[i] - offset), 2);
      sum += (vec_z[i] - offset);
    }
  }
  const Double_t k = s2sum / (n*s2sum-pow(sum,2));

  Double_t tel_res = k /(1+k) * measured * measured;
  Double_t single_point = TMath::Sqrt(measured * measured - k * tel_res);
  return TMath::Sqrt(tel_res);
}

// Get the telescope resolution from the measured input via k-factor
Double_t get_k(const Int_t pos)
{
  std::vector<Double_t> vec_z;
  for(int i=0; i < nplanes; i++)
  {
    vec_z.push_back(posx[i]);
  }

  const Double_t offset = vec_z[pos];
  const Double_t n = (Double_t)vec_z.size();
  Double_t s2sum = 0.0;
  Double_t sum = 0.0;
  for(Int_t i =0; i < nplanes;i++)
  {
    if(i != pos)
    {
      s2sum += pow((vec_z[i] - offset), 2);
      sum += (vec_z[i] - offset);
    }
  }
  const Double_t k = s2sum / ((n-1)*s2sum-pow(sum,2));

  return k;
}


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
/*double RMS_help(std::vector<double> vec, const int length)
  {
  double array[length];
  for(int i = 0; i < length; i++) array[i] = vec.at(i);
  return TMath::RMS(length, array);
  }

  double RMS(std::vector<double> vec)
  {
  return RMS_help(vec, vec.size()); 
  }*/

void run_res_offset( Double_t ebeam, Double_t *obsresol_x, Double_t* obsresol_error_x, Double_t* obsresol_y, Double_t* obsresol_error_y )
{
  v_residual_offset.clear();

  std::cout << "Run offset" << std::endl;
  // Create telescope
  AnaTel *tscope = new AnaTel(telescopebuild.c_str(), ebeam);
  std::cout << "Tscope created" << std::endl;
  // Give original beam energy with spread // assume spread to be negligible
  tscope->SetBeam( ebeam, 0.0 );

  // calculate residual offset to achieve/force sigma_int = 3.42e-3
  double sigma_int_est = 3.42e-3;
  tscope->SetResolution( sigma_int_est);

  for(Int_t i = 0; i < nplanes; i++ )
  {
    double residual_offsetx = -( pow(sigma_int_est,2) + pow(tscope->GetPointingResGBL(i,0),2) -pow(obsresol_x[i] ,2)) / (2.0*obsresol_x[i]);
    double residual_offsety = -( pow(sigma_int_est,2) + pow(tscope->GetPointingResGBL(i,0),2) -pow(obsresol_y[i] ,2)) / (2.0*obsresol_y[i]);
    v_residual_offset.push_back(residual_offsetx);
    v_residual_offset.push_back(residual_offsety);
    std::cout << "Calculated needed offset in X = " << residual_offsetx << "   in y = " << residual_offsety << std::endl;
  }

  double mean_res_offset = 0.0;
  double mean_res_offset_inner = 0.0;
  for(Int_t i = 0; i < 2*nplanes; i++ )
  {
    mean_res_offset += v_residual_offset.at(i)/(2*nplanes);
  }
  for(Int_t i = 2; i < 2*nplanes - 2; i++ )
  {
    mean_res_offset_inner += v_residual_offset.at(i)/(2*nplanes - 4);
  }

  std::cout << " ---- mean = " << mean_res_offset << std::endl;
  std::cout << " ---- mean inner = " << mean_res_offset_inner << std::endl;
  v_mean_res_offset.push_back(mean_res_offset);
  v_mean_res_offset_inner.push_back(mean_res_offset_inner);
}

// Run minuit
void run_global( Double_t ebeam, Double_t *obsresol_x, Double_t* obsresol_error_x, Double_t* obsresol_y, Double_t* obsresol_error_y )
{
  std::cout << "Run global" << std::endl;
  // Create telescope
  AnaTel *tlb = new AnaTel(telescopebuild.c_str(), ebeam);
  std::cout << "Tscope created" << std::endl;
  // Give original beam energy with spread // assume spread to be negligible
  tlb->SetBeam( ebeam, 0.0 );
  //tlb->SetResolution(3.42e-3);
  //cout << " E = " << ebeam << std::endl;
  //cout << " Pointing reso estimation at plane 0 using given initial reso (" << *(tlb->GetResolution()) << ") = " << tlb->GetPointingRes(0,0) << "   Residual estimate = " << tlb->GetWidth(0,0) << std::endl;
  //cout << " Pointing reso estimation at plane 3 using given initial reso (" << *(tlb->GetResolution()) << ") = " << tlb->GetPointingRes(3,0) << "   Residual estimate = " << tlb->GetWidth(3,0) << std::endl;
  //cout << " GBL Pointing reso estimation at plane 3 using given initial reso (" << *(tlb->GetResolution()) << ") = " << tlb->GetPointingResGBL(3,0) << "   Residual estimate = " << tlb->GetWidthGBL(3,0) << std::endl;
  //cout << " GBL Pointing reso estimation at plane 0 using given initial reso (" << *(tlb->GetResolution()) << ") = " << tlb->GetPointingResGBL(0,0) << "   Residual estimate = " << tlb->GetWidthGBL(0,0) << std::endl;

  fcn_beam = new MyFunctionObject2D_beam(tlb, obsresol_x, obsresol_error_x, obsresol_y, obsresol_error_y );
  /*
     printf( "Plane 0  %8.5f %8.5f \n", tlb->GetWidth(0,0), tlb->GetWidth(0,1) );
     printf( "Plane 1  %8.5f %8.5f \n", tlb->GetWidth(1,0), tlb->GetWidth(1,1) );
     printf( "Plane 2  %8.5f %8.5f \n", tlb->GetWidth(2,0), tlb->GetWidth(2,1) );
     printf( "Plane 3  %8.5f %8.5f \n", tlb->GetWidth(3,0), tlb->GetWidth(3,1) );
     printf( "Plane 4  %8.5f %8.5f \n", tlb->GetWidth(4,0), tlb->GetWidth(4,1) );
     printf( "Plane 5  %8.5f %8.5f \n", tlb->GetWidth(5,0), tlb->GetWidth(5,1) );
     */

  /*
     bool firstminuitcall = true;
     if(firstminuitcall)
     {
     gSystem->Load("libMinuit");//is this really needed?
     firstminuitcall = false;
     }*/

  //initialize TMinuit with a maximum of 4 params
  TMinuit *gMinuit = new TMinuit(4);

  // set print level (-1 = quiet, 0 = normal, 1 = verbose)
  gMinuit->SetPrintLevel(0);

  // give the function
  bool AverageOverMimosas = true;
  if(AverageOverMimosas) {
    gMinuit->SetFCN(fcn_wrapper);
    double arglist[10];
    int ierflg = 0;

    // minimization strategy (1 = standard, 2 = slower)
    arglist[0] = 2;
    gMinuit->mnexcm("SET STR",arglist,2,ierflg);

    // set error definition (1 = for chi square)
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

    // Set starting values and step sizes for parameters
    Double_t vstart[4] = {0.0040, ebeam,   0.0, 0.05 };
    Double_t step[4]   = {0.0001,   0.1, 0.001, 0.001};

    /* orig:
       gMinuit->mnparm(0, "m26resol" , vstart[0], step[0],    0.0001,      0.020, ierflg);
       gMinuit->mnparm(1, "pbeam"    , vstart[1], step[1], 0.1*ebeam,    2*ebeam, ierflg);
       gMinuit->mnparm(2, "spread"   , vstart[2], step[2],      0.00,        0.1, ierflg);
       gMinuit->mnparm(3, "thickness", vstart[3], step[3],      0.05,      0.100, ierflg);

*/
    gMinuit->mnparm(0, "m26resol" , vstart[0], step[0],    0.001,       0.020, ierflg);
    gMinuit->mnparm(1, "pbeam"    , vstart[1], step[1], 0.1*ebeam,  10.*ebeam, ierflg);
    gMinuit->mnparm(2, "spread"   , vstart[2], step[2],     -0.1,         0.1, ierflg);
    gMinuit->mnparm(3, "thickness", vstart[3], step[3],     0.04,        0.20, ierflg);

    //gMinuit->FixParameter(0);
    gMinuit->FixParameter(1);
    gMinuit->FixParameter(2);
    gMinuit->FixParameter(3);

    // Now ready for minimization step
    arglist[0] = 1000 ;
    arglist[1] = 0.001 ;
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
    //gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
    //gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);

    //gMinuit->Release(2);
    //gMinuit->Release(3);

    ////  Now ready for minimization step
    //arglist[0] = 8000;
    //arglist[1] = 1.0;
    //gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);


    //           calculate errors using MINOS. do we need this?
    //           arglist[0] = 2000;
    //           arglist[1] = 0.1;
    //           gMinuit->mnexcm("MINOS",arglist,1,ierflg);

    //   get results from migrad
    double par[4];
    double perr[4];
    gMinuit->GetParameter(0, par[0], perr[0]);
    gMinuit->GetParameter(1, par[1], perr[1]);
    gMinuit->GetParameter(2, par[2], perr[2]);
    gMinuit->GetParameter(3, par[3], perr[3]);



    /*          printf("%s 0 %8.4f %8.4f\n", par[0], perr[0] );
		printf("%s 1 %8.4f %8.4f\n", par[1], perr[1] );
		printf("%s 2 %8.4f %8.4f\n", par[2], perr[2] );
		printf("%s 3 %8.4f %8.4f\n", par[3], perr[3] );

		printf("%s chi2 %8.4f \n", fcn_beam->getchi2_2D(par, perr) );
		printf("%s tel_resolution at 0 %8.4f \n", resol_estimate(par[0],0));
		printf("%s tel_resolution at 1 %8.4f \n", resol_estimate(par[0],1));
		printf("%s tel_resolution at 2 %8.4f \n", resol_estimate(par[0],2));
		printf("%s tel_resolution at 3 %8.4f \n", resol_estimate(par[0],3));
		printf("%s tel_resolution at 4 %8.4f \n", resol_estimate(par[0],4));
		printf("%s tel_resolution at 5 %8.4f \n", resol_estimate(par[0],5));
		*/

    std::cout << "Resolution was " << m26_resolution*1000.0 << " mu m" << std::endl;
    m26_resolution = par[0];
    std::cout << "Resolution is " << m26_resolution*1000.0 << " mu m" << std::endl;
    std::cout << "Resolution error was " << m26_res_error*1000.0 << " mu m" << std::endl;
    m26_res_error  = perr[0];
    std::cout << "Resolution error is " << m26_res_error*1000.0 << " mu m" << std::endl;
    std::cout << "Beam was " << global_beam << std::endl;
    global_beam      = par[1];
    std::cout << "global beam now " << global_beam << std::endl;
    std::cout << "Spread was " << global_spread << std::endl;
    global_spread    = par[2];
    std::cout << "Spread is " << global_spread << std::endl;
    std::cout << "Thickness was " << global_thickness << std::endl;
    global_thickness = par[3];
    std::cout << "Thickness is " << global_thickness << std::endl;

    Double_t tempdist = posx[1] - posx[0];
    Double_t scatterer= 0.0136/global_beam * sqrt(global_thickness/93.66 + tempdist/304200.0 + 0.1/286.0) * (1.+0.038*std::log(global_thickness/93.66 + tempdist/304200.0 + 0.05/286.0)) ;

    std::cout << "Scattering is " << scatterer << std::endl;
  } else {
    gMinuit->SetFCN(fcn_wrapper_nomean);
    double arglist[10];
    int ierflg = 0;

    // minimization strategy (1 = standard, 2 = slower)
    arglist[0] = 2;
    gMinuit->mnexcm("SET STR",arglist,2,ierflg);

    // set error definition (1 = for chi square)
    arglist[0] = 1;
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

    // Set starting values and step sizes for parameters
    Double_t vstart[9] = {0.005,0.005,0.005,0.005,0.005,0.005,      ebeam,   0.0, 0.05 };
    Double_t step[9]   = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,  0.1, 0.001, 0.001};

    gMinuit->mnparm(0, "m26resol0" , vstart[0], step[0],    0.001,       0.020, ierflg);
    gMinuit->mnparm(1, "m26resol1" , vstart[1], step[1],    0.001,       0.020, ierflg);
    gMinuit->mnparm(2, "m26resol2" , vstart[2], step[2],    0.001,       0.020, ierflg);
    gMinuit->mnparm(3, "m26resol3" , vstart[3], step[3],    0.001,       0.020, ierflg);
    gMinuit->mnparm(4, "m26resol4" , vstart[4], step[4],    0.001,       0.020, ierflg);
    gMinuit->mnparm(5, "m26resol5" , vstart[5], step[5],    0.001,       0.020, ierflg);
    gMinuit->mnparm(6, "pbeam"    ,  vstart[6], step[6], 0.1*ebeam,  10.*ebeam, ierflg);
    gMinuit->mnparm(7, "spread"   ,  vstart[7], step[7],    -0.01,         0.1, ierflg);
    gMinuit->mnparm(8, "thickness",  vstart[8], step[8],     0.04,        0.06, ierflg);

    //gMinuit->FixParameter(0);
    gMinuit->FixParameter(6);
    gMinuit->FixParameter(7);
    gMinuit->FixParameter(8);

    // Now ready for minimization step
    arglist[0] = 1000 ;
    arglist[1] = 0.001 ;
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
    //gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
    //gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);

    //gMinuit->Release(2);
    //gMinuit->Release(3);

    ////  Now ready for minimization step
    //arglist[0] = 8000;
    //arglist[1] = 1.0;
    //gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);


    //           calculate errors using MINOS. do we need this?
    //           arglist[0] = 2000;
    //           arglist[1] = 0.1;
    //           gMinuit->mnexcm("MINOS",arglist,1,ierflg);

    //   get results from migrad
    double par[9];
    double perr[9];
    gMinuit->GetParameter(0, par[0], perr[0]);
    gMinuit->GetParameter(1, par[1], perr[1]);
    gMinuit->GetParameter(2, par[2], perr[2]);
    gMinuit->GetParameter(3, par[3], perr[3]);
    gMinuit->GetParameter(4, par[4], perr[4]);
    gMinuit->GetParameter(5, par[5], perr[5]);
    gMinuit->GetParameter(6, par[6], perr[6]);
    gMinuit->GetParameter(7, par[7], perr[7]);
    gMinuit->GetParameter(8, par[8], perr[8]);



    /*          printf("%s 0 %8.4f %8.4f\n", par[0], perr[0] );
		printf("%s 1 %8.4f %8.4f\n", par[1], perr[1] );
		printf("%s 2 %8.4f %8.4f\n", par[2], perr[2] );
		printf("%s 3 %8.4f %8.4f\n", par[3], perr[3] );

		printf("%s chi2 %8.4f \n", fcn_beam->getchi2_2D(par, perr) );
		printf("%s tel_resolution at 0 %8.4f \n", resol_estimate(par[0],0));
		printf("%s tel_resolution at 1 %8.4f \n", resol_estimate(par[0],1));
		printf("%s tel_resolution at 2 %8.4f \n", resol_estimate(par[0],2));
		printf("%s tel_resolution at 3 %8.4f \n", resol_estimate(par[0],3));
		printf("%s tel_resolution at 4 %8.4f \n", resol_estimate(par[0],4));
		printf("%s tel_resolution at 5 %8.4f \n", resol_estimate(par[0],5));
		*/

    std::cout << "Resolution 0 is " << par[0]*1000.0 << " mu m" << std::endl;
    std::cout << "Resolution 1 is " << par[1]*1000.0 << " mu m" << std::endl;
    std::cout << "Resolution 2 is " << par[2]*1000.0 << " mu m" << std::endl;
    std::cout << "Resolution 3 is " << par[3]*1000.0 << " mu m" << std::endl;
    std::cout << "Resolution 4 is " << par[4]*1000.0 << " mu m" << std::endl;
    std::cout << "Resolution 5 is " << par[5]*1000.0 << " mu m" << std::endl;
    m26_resolution = (par[0] + par[1] + par[2] + par[3] + par[4] + par[5])/6.;
    std::cout << "Avg. resolution is " << m26_resolution*1000.0 << " mu m" << std::endl;

  }

}

void filehelper(std::string &help, Int_t RN)
{
  if(RN < 113 || RN > 700) help = inputDir150 + inputFile;
  else help = inputDir20 + inputFile;
  //std::cout << "help: " << help << ", runnumber = " << RN << std::endl;
  //std::cout << "run number = " << RN << std::endl;
  //TString filename(help.c_str());
  if (RN <= 999)
    help+="0";
  if (RN <= 99)
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

  /*h_clustersize[0] = (TH1D*) (tfile->Get("FilterHisto/detector_0/clusterSignal_d0"))->Clone();
    h_clustersize[1] = (TH1D*) (tfile->Get("FilterHisto/detector_1/clusterSignal_d1"))->Clone();
    h_clustersize[2] = (TH1D*) (tfile->Get("FilterHisto/detector_2/clusterSignal_d2"))->Clone();
    h_clustersize[3] = (TH1D*) (tfile->Get("FilterHisto/detector_3/clusterSignal_d3"))->Clone();
    h_clustersize[4] = (TH1D*) (tfile->Get("FilterHisto/detector_4/clusterSignal_d4"))->Clone();
    h_clustersize[5] = (TH1D*) (tfile->Get("FilterHisto/detector_5/clusterSignal_d5"))->Clone();*/

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

  TCanvas *canv_unb_pulls;
  canv_name = "m26fitter_unb_pulls_";
  canv_name += std::to_string(runnumber);
  canv_unb_pulls = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canv_unb_pulls);
  canv_unb_pulls->Divide(2,6);

  TCanvas *canv_kink_pulls;
  canv_name = "m26fitter_kink_pulls_";
  canv_name += std::to_string(runnumber);
  canv_kink_pulls = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canv_kink_pulls);
  canv_kink_pulls->Divide(2,4);

  TCanvas *canv_kink_unb_pulls;
  canv_name = "m26fitter_kink_unb_pulls_";
  canv_name += std::to_string(runnumber);
  canv_kink_unb_pulls = new TCanvas(canv_name.c_str(),canv_name.c_str(),900,10,600,800); 
  CanvasSetter(canv_kink_unb_pulls);
  canv_kink_unb_pulls->Divide(2,4);





  // histos
  TH1D *h_m26[planescount];
  TH1D *h_m26_biased[planescount];
  TH1D *h_pull_biased[planescount];
  TH1D *h_pull_unbiased[planescount];
  TH1D *h_kink_pull[8];
  TH1D *h_kink_unb_pull[8];
  tfile->cd();

  // If the clustersize is limited, put clustersize 1 into the original histogram
  if(clusterlimit >0) submask="1";

  std::cout << "clone histos";
  // Load file
  //
  /*
     h_m26[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftX"+submask))->Clone();
     h_m26[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftY"+submask))->Clone();
     h_m26[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftX"+submask))->Clone();
     h_m26[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftY"+submask))->Clone();
     h_m26[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftX"+submask))->Clone();
     h_m26[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftY"+submask))->Clone();
     h_m26[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX"+submask))->Clone();
     h_m26[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY"+submask))->Clone();
     h_m26[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftX"+submask))->Clone();
     h_m26[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftY"+submask))->Clone();
     h_m26[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftX"+submask))->Clone();
     h_m26[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftY"+submask))->Clone();
     */
  /*h_m26[0] = (TH1D*) (tfile->Get("Fitter00/pl0_residualX"))->Clone();
    h_m26[1] = (TH1D*) (tfile->Get("Fitter00/pl0_residualY"))->Clone();
    h_m26[2] = (TH1D*) (tfile->Get("Fitter01/pl1_residualX"))->Clone();
    h_m26[3] = (TH1D*) (tfile->Get("Fitter01/pl1_residualY"))->Clone();
    h_m26[4] = (TH1D*) (tfile->Get("Fitter02/pl2_residualX"))->Clone();
    h_m26[5] = (TH1D*) (tfile->Get("Fitter02/pl2_residualY"))->Clone();
    h_m26[6] = (TH1D*) (tfile->Get("Fitter03/pl3_residualX"))->Clone();
    h_m26[7] = (TH1D*) (tfile->Get("Fitter03/pl3_residualY"))->Clone();
    h_m26[8] = (TH1D*) (tfile->Get("Fitter04/pl4_residualX"))->Clone();
    h_m26[9] = (TH1D*) (tfile->Get("Fitter04/pl4_residualY"))->Clone();
    h_m26[10] = (TH1D*)(tfile->Get("Fitter05/pl5_residualX"))->Clone();
    h_m26[11] = (TH1D*)(tfile->Get("Fitter05/pl5_residualY"))->Clone();
    */

  h_m26[0] = (TH1D*) (tfile->Get("Fitter00/GBL/gblrx0"))->Clone();
  h_m26[1] = (TH1D*) (tfile->Get("Fitter00/GBL/gblry0"))->Clone();
  h_m26[2] = (TH1D*) (tfile->Get("Fitter01/GBL/gblrx1"))->Clone();
  h_m26[3] = (TH1D*) (tfile->Get("Fitter01/GBL/gblry1"))->Clone();
  h_m26[4] = (TH1D*) (tfile->Get("Fitter02/GBL/gblrx2"))->Clone();
  h_m26[5] = (TH1D*) (tfile->Get("Fitter02/GBL/gblry2"))->Clone();
  h_m26[6] = (TH1D*) (tfile->Get("Fitter03/GBL/gblrx3"))->Clone();
  h_m26[7] = (TH1D*) (tfile->Get("Fitter03/GBL/gblry3"))->Clone();
  h_m26[8] = (TH1D*) (tfile->Get("Fitter04/GBL/gblrx4"))->Clone();
  h_m26[9] = (TH1D*) (tfile->Get("Fitter04/GBL/gblry4"))->Clone();
  h_m26[10] = (TH1D*)(tfile->Get("Fitter05/GBL/gblrx5"))->Clone();
  h_m26[11] = (TH1D*)(tfile->Get("Fitter05/GBL/gblry5"))->Clone();

  h_m26_biased[0] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx0"))->Clone();
  h_m26_biased[1] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry0"))->Clone();
  h_m26_biased[2] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx1"))->Clone();
  h_m26_biased[3] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry1"))->Clone();
  h_m26_biased[4] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx2"))->Clone();
  h_m26_biased[5] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry2"))->Clone();
  h_m26_biased[6] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx3"))->Clone();
  h_m26_biased[7] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry3"))->Clone();
  h_m26_biased[8] = (TH1D*) (tfile->Get("Fitter06/GBL/gblrx4"))->Clone();
  h_m26_biased[9] = (TH1D*) (tfile->Get("Fitter06/GBL/gblry4"))->Clone();
  h_m26_biased[10] = (TH1D*)(tfile->Get("Fitter06/GBL/gblrx5"))->Clone();
  h_m26_biased[11] = (TH1D*)(tfile->Get("Fitter06/GBL/gblry5"))->Clone();

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

  h_pull_unbiased[0] = (TH1D*) (tfile->Get("Fitter00/GBL/gblpx0_unb"))->Clone();
  h_pull_unbiased[1] = (TH1D*) (tfile->Get("Fitter00/GBL/gblpy0_unb"))->Clone();
  h_pull_unbiased[2] = (TH1D*) (tfile->Get("Fitter01/GBL/gblpx1_unb"))->Clone();
  h_pull_unbiased[3] = (TH1D*) (tfile->Get("Fitter01/GBL/gblpy1_unb"))->Clone();
  h_pull_unbiased[4] = (TH1D*) (tfile->Get("Fitter02/GBL/gblpx2_unb"))->Clone();
  h_pull_unbiased[5] = (TH1D*) (tfile->Get("Fitter02/GBL/gblpy2_unb"))->Clone();
  h_pull_unbiased[6] = (TH1D*) (tfile->Get("Fitter03/GBL/gblpx3_unb"))->Clone();
  h_pull_unbiased[7] = (TH1D*) (tfile->Get("Fitter03/GBL/gblpy3_unb"))->Clone();
  h_pull_unbiased[8] = (TH1D*) (tfile->Get("Fitter04/GBL/gblpx4_unb"))->Clone();
  h_pull_unbiased[9] = (TH1D*) (tfile->Get("Fitter04/GBL/gblpy4_unb"))->Clone();
  h_pull_unbiased[10] = (TH1D*)(tfile->Get("Fitter05/GBL/gblpx5_unb"))->Clone();
  h_pull_unbiased[11] = (TH1D*)(tfile->Get("Fitter05/GBL/gblpy5_unb"))->Clone();

  h_kink_pull[0] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx1"))->Clone();
  //h_kink_pull[1] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty1"))->Clone();
  h_kink_pull[2] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx2"))->Clone();
  //h_kink_pull[3] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty2"))->Clone();
  h_kink_pull[4] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx3"))->Clone();
  //h_kink_pull[5] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty3"))->Clone();
  h_kink_pull[6] = (TH1D*) (tfile->Get("Fitter06/GBL/gbltx4"))->Clone();
  //h_kink_pull[7] = (TH1D*) (tfile->Get("Fitter06/GBL/gblty4"))->Clone();

  h_kink_unb_pull[0] = (TH1D*) (tfile->Get("Fitter01/GBL/gbltx1"))->Clone();
  //h_kink_unb_pull[1] = (TH1D*) (tfile->Get("Fitter01/GBL/gblty1"))->Clone();
  h_kink_unb_pull[2] = (TH1D*) (tfile->Get("Fitter02/GBL/gbltx2"))->Clone();
  //h_kink_unb_pull[3] = (TH1D*) (tfile->Get("Fitter02/GBL/gblty2"))->Clone();
  h_kink_unb_pull[4] = (TH1D*) (tfile->Get("Fitter03/GBL/gbltx3"))->Clone();
  //h_kink_unb_pull[5] = (TH1D*) (tfile->Get("Fitter03/GBL/gblty3"))->Clone();
  h_kink_unb_pull[6] = (TH1D*) (tfile->Get("Fitter04/GBL/gbltx4"))->Clone();
  //h_kink_unb_pull[7] = (TH1D*) (tfile->Get("Fitter04/GBL/gblty4"))->Clone();

  std::cout << "  - done" << std::endl;
  // Add the clustersizes 2-4 into the histos
  //
  // --- No clustersize dependence possible at the moment with GBL tracks
  if(clusterlimit >0)
  {
    std::cout << "Clustersize is limited. This is only implemented up to 4!" << std::endl;
    TH1D *h_m26_2[planescount];
    TH1D *h_m26_3[planescount];
    TH1D *h_m26_4[planescount];
    h_m26_2[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftX2"))->Clone();
    h_m26_2[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftY2"))->Clone();
    h_m26_2[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftX2"))->Clone();
    h_m26_2[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftY2"))->Clone();
    h_m26_2[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftX2"))->Clone();
    h_m26_2[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftY2"))->Clone();
    h_m26_2[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX2"))->Clone();
    h_m26_2[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY2"))->Clone();
    h_m26_2[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftX2"))->Clone();
    h_m26_2[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftY2"))->Clone();
    h_m26_2[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftX2"))->Clone();
    h_m26_2[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftY2"))->Clone();

    h_m26_3[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftX3"))->Clone();
    h_m26_3[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftY3"))->Clone();
    h_m26_3[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftX3"))->Clone();
    h_m26_3[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftY3"))->Clone();
    h_m26_3[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftX3"))->Clone();
    h_m26_3[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftY3"))->Clone();
    h_m26_3[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX3"))->Clone();
    h_m26_3[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY3"))->Clone();
    h_m26_3[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftX3"))->Clone();
    h_m26_3[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftY3"))->Clone();
    h_m26_3[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftX3"))->Clone();
    h_m26_3[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftY3"))->Clone();

    h_m26_4[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftX4"))->Clone();
    h_m26_4[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftY4"))->Clone();
    h_m26_4[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftX4"))->Clone();
    h_m26_4[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftY4"))->Clone();
    h_m26_4[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftX4"))->Clone();
    h_m26_4[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftY4"))->Clone();
    h_m26_4[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX4"))->Clone();
    h_m26_4[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY4"))->Clone();
    h_m26_4[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftX4"))->Clone();
    h_m26_4[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftY4"))->Clone();
    h_m26_4[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftX4"))->Clone();
    h_m26_4[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftY4"))->Clone();

    // Add the clustersize 2-4 histos into the "original" one and continue
    for(int i=0;i<planescount;i++)
    {
      h_m26_3[i]->Add(h_m26_4[i]);
      h_m26_2[i]->Add(h_m26_3[i]);
      h_m26[i]->Add(h_m26_2[i]);
    }
  }

  //  range
  for(int i=0;i<planescount;i++){
    h_m26[i]->GetXaxis()->SetRangeUser(-100.,100);
    h_m26_biased[i]->GetXaxis()->SetRangeUser(-100.,100);
  }

  // TCanvas *canv2;
  // canv2 = new TCanvas("m26fitter2","m21fitter(",900,10,600,800); 
  // canv2->SetFillColor(0);
  // canv2->Divide(2,3);

  if (verbose1)
  {
    std::cout << "\n  Reading residuals from file...\n" << std::endl;
  }

  std::vector<double> vec_res_sigma_unbiasedX;
  std::vector<double> vec_res_sigma_unbiasedY;
  std::vector<double> vec_res_sigma_biasedX;
  std::vector<double> vec_res_sigma_biasedY;
  std::vector<double> vec_pull_meanX;
  std::vector<double> vec_pull_meanY;
  std::vector<double> vec_pull_sigmaX;
  std::vector<double> vec_pull_sigmaY;
  std::vector<double> vec_unb_pull_meanX;
  std::vector<double> vec_unb_pull_meanY;
  std::vector<double> vec_unb_pull_sigmaX;
  std::vector<double> vec_unb_pull_sigmaY;
  std::vector<double> vec_kink_pull_meanX;
  std::vector<double> vec_kink_pull_sigmaX;
  std::vector<double> vec_kink_unb_pull_meanX;
  std::vector<double> vec_kink_unb_pull_sigmaX;

  // Sort into X or Y
  for(Int_t i = 0; i < nplanes*2; i++ )
  {
    canv->cd(i+1);
    h_m26_biased[i]->Draw();
    gPad->Update();
    TPaveStats * ps = (TPaveStats *)h_m26_biased[i]->FindObject("stats");
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

    if(plot_residuals)
    {
      TString xtitle;
      if(i == 0 || i == 2 || i == 4 ||  i == 6 || i == 8 || i == 10)
	xtitle = "(x_{pred.} - x_{meas.}) / mm";
      else
	xtitle = "(y_{pred.} - y_{meas.}) / mm";

      //format_mc(h_m26[i], 1, kBlue);
      //histo_cfg(h_m26[i],xtitle , "tracks", st);
      h_m26[i]->GetYaxis()->SetNoExponent();
      h_m26[i]->GetYaxis()->SetTitleOffset(1.68);
      h_m26[i]->GetXaxis()->SetTitleOffset(0.97);
      h_m26[i]->GetYaxis()->SetTitleSize(0.035);
      h_m26[i]->GetXaxis()->SetTitleSize(0.035);
      h_m26[i]->GetYaxis()->SetLabelSize(0.035);
      h_m26[i]->GetXaxis()->SetLabelSize(0.035);
      h_m26[i]->GetXaxis()->SetLabelOffset(0.004);
      //      gPad->SetLeftMargin(0.2238535);
      h_m26[i]->SetMaximum(h_m26[i]->GetMaximum()*1.3);
      h_m26[i]->DrawCopy("hist");
    }

    // Fit the residuals
    TF1 *f1 = new TF1("f1","gaus");
    TF1 *f1_biased = new TF1("f1_biased","gaus");

    float mean1 = h_m26[i]->GetMean();
    float rms1 = h_m26[i]->GetRMS();

    float mean_biased = h_m26_biased[i]->GetMean();
    float rms_biased  = h_m26_biased[i]->GetRMS();

    //std::cout << "Mean = " << mean1 << "   RMS = " << rms1 << std::endl;

    f1->SetParameter(1,mean1);
    f1->SetParameter(2,rms1);
    f1->SetParLimits(2,rms1*0.5,rms1*2.);
    //f1->FixParameter(2,16.);
    f1_biased->SetParameter(1,mean_biased);
    f1_biased->SetParameter(2,rms_biased);

    h_m26[i]->Fit(f1,"QB");
    h_m26_biased[i]->Fit(f1_biased,"EMQI0");

    Double_t mean_1 = f1->GetParameter(1);
    Double_t sigma_1 = f1->GetParameter(2);
    Double_t sigma_error_1 = f1->GetParError(2);


    //std::cout << " --- simga from first fit = " << sigma_1 << std::endl;

    Double_t mean_1_biased = f1_biased->GetParameter(1);
    Double_t sigma_1_biased = f1_biased->GetParameter(2);
    Double_t sigma_error_1_biased = f1_biased->GetParError(2);

    // Repeat within 2 sigma

    float nsig = 2.;
    float low_lim = mean_1        - nsig*sigma_1;
    float high_lim = mean_1        + nsig*sigma_1;
    if(low_lim < -25.) low_lim = -24.9;
    if(high_lim > 25.) high_lim = 24.9;
    float low_lim_biased = mean_1_biased        - nsig*sigma_1_biased;
    float high_lim_biased = mean_1_biased        + nsig*sigma_1_biased;
    TF1 *f2 = new TF1("f2","gaus", low_lim, high_lim);
    TF1 *f2_biased = new TF1("f2_biased","gaus", low_lim_biased, high_lim_biased);
    f2->SetLineWidth(2);
    f2->SetLineStyle(1);
    f2->SetLineColor(kBlack);
    f2_biased->SetLineWidth(2);
    f2_biased->SetLineStyle(1);
    f2_biased->SetLineColor(kBlack);
    h_m26[i]->       Fit(f2,"EQMIR", "");
    h_m26_biased[i]->Fit(f2_biased,"EQMIR","");
    Double_t mean_2 = f2->GetParameter(1);
    Double_t emean_2 = f2->GetParError(1);
    Double_t mean_2_biased = f2_biased->GetParameter(1);
    Double_t emean_2_biased = f2_biased->GetParError(1);

    if (verbose0)
      std::cout << " Measured residual mean plane " << j << " is " << (mean_2 ) << " mu m" << std::endl; // removed 1e3

    Double_t sigma_2 = f2->GetParameter(2);
    Double_t esigma_2 = f2->GetParError(2);
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
    Double_t sigma_error_2 = 0.05 * sigma_2;
    Double_t sigma_error_2_biased = 0.05 * sigma_2;

    // Now get resolution estimate
    //Double_t tel_resol = resol_estimate(sigma_2,j);

    //if (verbose0)
    //cout << " Measured residual width plane " << j << " is " <<  tel_resol << " mu m" << std::endl;

    // Fill fitted resolution into array
    if(i == 0 || i == 2 || i == 4 ||  i == 6 || i == 8 || i == 10)
    {
      obsresol_x[j] = sigma_2*1.e-3;
      obsresol_error_x[j] = sigma_error_2*1.e-3;
      if(IsBiased)
      {
	obsresol_x[j] = sigma_2_biased*1.e-3;
	obsresol_error_x[j] = sigma_error_2_biased*1.e-3;
      }
    }
    else
    {
      obsresol_y[j] = sigma_2*1.e-3; 
      obsresol_error_y[j] = sigma_error_2*1.e-3;
      if(IsBiased)
      {
	obsresol_y[j] = sigma_2_biased*1.e-3;
	obsresol_error_y[j] = sigma_error_2_biased*1.e-3;
      }
    }

    //std::cout << " perform fits " << std::endl;
    double mean = 0.;
    double sigma = 0.;

    canv_pulls->cd(i+1);
    DoGausFit(h_pull_biased[i], mean, sigma);
    canv_pulls->Modified();
    canv_pulls->Update();

    if(i%2 == 0) {
      vec_res_sigma_unbiasedX.push_back(sigma_2);
      vec_res_sigma_biasedX.push_back(sigma_2_biased);
      vec_pull_meanX.push_back(mean);
      vec_pull_sigmaX.push_back(sigma);
    }
    if(i%2 == 1) {
      vec_res_sigma_unbiasedY.push_back(sigma_2);
      vec_res_sigma_biasedY.push_back(sigma_2_biased);
      vec_pull_meanY.push_back(mean);
      vec_pull_sigmaY.push_back(sigma);
    }

    mean = 0;
    sigma = 0;
    canv_unb_pulls->cd(i+1);
    DoGausFit(h_pull_unbiased[i], mean, sigma);
    canv_unb_pulls->Modified();
    canv_unb_pulls->Update();

    if(i%2 == 0) {
      vec_unb_pull_meanX.push_back(mean);
      vec_unb_pull_sigmaX.push_back(sigma);
    }
    if(i%2 == 1) {
      vec_unb_pull_meanY.push_back(mean);
      vec_unb_pull_sigmaY.push_back(sigma);

    }

    mean = 0;
    sigma = 0;
    if( i > 1 && i < 10){
      canv_kink_pulls->cd(i+1);
      if(i%2 == 0) {
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
    if( i > 1 && i < 10){
      canv_kink_unb_pulls->cd(i+1);
      if(i%2 == 0) {
	DoGausFit(h_kink_unb_pull[i-2], mean, sigma);
	vec_kink_unb_pull_meanX.push_back(mean);
	vec_kink_unb_pull_sigmaX.push_back(sigma);
      }
      canv_kink_unb_pulls->Modified();
      canv_kink_unb_pulls->Update();
    }

    //if(i%2 == 1) {
    //  vec_kinks_meanY.push_back(mean);
    //  vec_kinks_sigmaY.push_back(sigma);
    //}



    // Plot this
    if(plot_residuals)
    {
      TPaveText *label;
      label = new TPaveText(0.2588087,0.7992021,0.6363255,0.8823138, "brNDC");
      label->SetTextAlign(11);
      label->SetTextFont(22);
      label->SetTextSize(0.08);
      label->SetFillStyle(0);
      label->SetBorderSize(0);
      TString sigmatext = "";
      char tmpstring2[500];
      sprintf(tmpstring2, "#sigma = (%1.3f #pm %1.3f) #mum", sigma_2, sigma_error_2);
      sigmatext = tmpstring2;
      label->AddText(sigmatext);
      label->Draw();
    }
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
    if(!IsBiased)
    {
      v_mean.push_back(mean_2);
      v_emean.push_back(emean_2);
      v_sigma.push_back(sigma_2);
      v_esigma.push_back(esigma_2);
    } else {
      v_mean.push_back(mean_2_biased);
      v_emean.push_back(emean_2_biased);
      v_sigma.push_back(sigma_2_biased);
      v_esigma.push_back(esigma_2_biased);
    }

    delete f1;
    delete f1_biased;
    delete f2;
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
  std::cout  << "   avgX = " << avg_res_sigma;
  DoRMS(vec_res_sigma_biasedX,rms_res_sigma);
  std::cout  << "   rmsX = " << rms_res_sigma << std::endl;

  double avg_unb_pull_sigma = .0;
  double rms_unb_pull_sigma = .0;
  std::cout << " vec_unb_pull_sigmaX = ";
  for(unsigned int i = 0; i<6; i++){
    std::cout << vec_unb_pull_sigmaX.at(i) << " ";
    avg_unb_pull_sigma += vec_unb_pull_sigmaX.at(i)/6.;
  }
  std::cout  << " unb avgX = " << avg_unb_pull_sigma;
  DoRMS(vec_unb_pull_sigmaX,rms_unb_pull_sigma);
  std::cout  << "   rmsX = " << rms_unb_pull_sigma << std::endl;

  avg_unb_pull_sigma = .0;
  rms_unb_pull_sigma = .0;
  std::cout << " vec_unb_pull_sigmaY = ";
  for(unsigned int i = 0; i<6; i++){
    std::cout << vec_unb_pull_sigmaY.at(i) << " ";
    avg_unb_pull_sigma += vec_unb_pull_sigmaY.at(i)/6.;
  }
  std::cout  << " unb avgY = " << avg_unb_pull_sigma;
  DoRMS(vec_unb_pull_sigmaY,rms_unb_pull_sigma);
  std::cout  << "   rmsY = " << rms_unb_pull_sigma << std::endl;


  /*float avg_sig_intX = .0;
    std::cout << " sqrt(vec_res_sigma_unbiasedX * vec_res_sigma_biasedX) = ";
    for(unsigned int i = 0; i<6; i++){
    std::cout << "sqrt( " << vec_res_sigma_unbiasedX.at(i) << " * " << vec_res_sigma_biasedX.at(i) << ") = " 
    << sqrt(vec_res_sigma_unbiasedX.at(i)*vec_res_sigma_biasedX.at(i)) << " ";
    avg_sig_intX += sqrt(vec_res_sigma_unbiasedX.at(i)*vec_res_sigma_biasedX.at(i))/6.;
    }
    std::cout  << "   avgX = " << avg_sig_intX << std::endl;

    float avg_sig_intY = .0;
    std::cout << " sqrt(vec_res_sigma_unbiasedY * vec_res_sigma_biasedY) = ";
    for(unsigned int i = 0; i<6; i++){
    std::cout << sqrt(vec_res_sigma_unbiasedY.at(i)*vec_res_sigma_biasedY.at(i)) << " ";
    avg_sig_intY += sqrt(vec_res_sigma_unbiasedY.at(i)*vec_res_sigma_biasedY.at(i))/6.;
    }
    std::cout  << "   avgY = " << avg_sig_intY << std::endl;
    */

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

  double avg_kink_unb_pull_sigma = .0;
  double rms_kink_unb_pull_sigma = .0;
  std::cout << " vec_kink_unb_pull_sigmaX = ";
  for(unsigned int i = 0; i<4; i++){
    std::cout << vec_kink_unb_pull_sigmaX.at(i) << " ";
    avg_kink_unb_pull_sigma += vec_kink_unb_pull_sigmaX.at(i)/4.;
  }
  std::cout  << " avg kink unb_pull X = " << avg_kink_unb_pull_sigma;
  DoRMS(vec_kink_unb_pull_sigmaX,rms_kink_unb_pull_sigma);
  std::cout  << "   rmsX = " << rms_kink_unb_pull_sigma << std::endl;



  v_runnumber.push_back((double)runnumber);
  v_erunnumber.push_back(0.0);

  // use this for biased
  if(IsBiased){
    if(runnumber > 113 && runnumber < 700) pulls20  << runnumber << " " << std::fixed << std::setprecision(6) <<  (vec_pull_sigmaX.at(0)+ vec_pull_sigmaY.at(0))/2. << " " <<  (vec_pull_sigmaX.at(1)+vec_pull_sigmaY.at(1))/2. << " " <<  (vec_pull_sigmaX.at(2)+vec_pull_sigmaY.at(2))/2. << " " <<  (vec_pull_sigmaX.at(3)+vec_pull_sigmaY.at(3))/2. << " " <<  (vec_pull_sigmaX.at(4)+vec_pull_sigmaY.at(4))/2. << " " <<  (vec_pull_sigmaX.at(5)+vec_pull_sigmaY.at(5))/2. << "\n";
    if(runnumber < 113 || runnumber > 700) pulls150  << runnumber << " " << std::fixed << std::setprecision(6) <<  (vec_pull_sigmaX.at(0)+ vec_pull_sigmaY.at(0))/2. << " " <<  (vec_pull_sigmaX.at(1)+vec_pull_sigmaY.at(1))/2. << " " <<  (vec_pull_sigmaX.at(2)+vec_pull_sigmaY.at(2))/2. << " " <<  (vec_pull_sigmaX.at(3)+vec_pull_sigmaY.at(3))/2. << " " <<  (vec_pull_sigmaX.at(4)+vec_pull_sigmaY.at(4))/2. << " " <<  (vec_pull_sigmaX.at(5)+vec_pull_sigmaY.at(5))/2. << "\n"; }
  else {

    // use this for unbiased
    if(runnumber > 113 && runnumber < 700) pulls20  << runnumber << " " << std::fixed << std::setprecision(6) <<  (vec_unb_pull_sigmaX.at(0)+ vec_unb_pull_sigmaY.at(0))/2. << " " <<  (vec_unb_pull_sigmaX.at(1)+vec_unb_pull_sigmaY.at(1))/2. << " " <<  (vec_unb_pull_sigmaX.at(2)+vec_unb_pull_sigmaY.at(2))/2. << " " <<  (vec_unb_pull_sigmaX.at(3)+vec_unb_pull_sigmaY.at(3))/2. << " " <<  (vec_unb_pull_sigmaX.at(4)+vec_unb_pull_sigmaY.at(4))/2. << " " <<  (vec_unb_pull_sigmaX.at(5)+vec_unb_pull_sigmaY.at(5))/2. << "\n";
    if(runnumber < 113 || runnumber > 700) pulls150  << runnumber << " " << std::fixed << std::setprecision(6) <<  (vec_unb_pull_sigmaX.at(0)+ vec_unb_pull_sigmaY.at(0))/2. << " " <<  (vec_unb_pull_sigmaX.at(1)+vec_unb_pull_sigmaY.at(1))/2. << " " <<  (vec_unb_pull_sigmaX.at(2)+vec_unb_pull_sigmaY.at(2))/2. << " " <<  (vec_unb_pull_sigmaX.at(3)+vec_unb_pull_sigmaY.at(3))/2. << " " <<  (vec_unb_pull_sigmaX.at(4)+vec_unb_pull_sigmaY.at(4))/2. << " " <<  (vec_unb_pull_sigmaX.at(5)+vec_unb_pull_sigmaY.at(5))/2. << "\n";
  }

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
  canv_unb_pulls->Write();
  canv_unb_pulls->Close();
  canv_kink_pulls->Write();
  canv_kink_pulls->Close();
  canv_kink_unb_pulls->Write();
  canv_kink_unb_pulls->Close();

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
  g_obsresidualX->SetMarkerColor(runnumber%10);
  g_obsresidualX->SetMarkerSize(1.2);
  g_obsresidualY->SetMarkerStyle(runnumber%10 + 20);
  g_obsresidualY->SetMarkerColor(runnumber%10);
  g_obsresidualY->SetMarkerSize(1.2);


  mg_obsresX->Add(g_obsresidualX,"p");
  mg_obsresY->Add(g_obsresidualY,"p");

  delete tfile;

  delete canv;
  delete canv_pulls;
  delete canv_unb_pulls;
  delete canv_kink_pulls;
  delete canv_kink_unb_pulls;


  /*
   * 
   * 
   * 

  // Create 4 new "telescopes"
  AnaTel *tl = new AnaTel(telescopebuild);
  tl->SetBeam( ebeam );
  AnaTel *tl2 = new AnaTel(telescopebuild);
  tl2->SetBeam( ebeam );
  AnaTel *tl2x = new AnaTel(telescopebuild);
  tl2x->SetBeam( ebeam );
  AnaTel *tl2y = new AnaTel(telescopebuild);
  tl2y->SetBeam( ebeam );
  tl->SetThickness( global_thickness );
  tl2->SetThickness( global_thickness );
  tl2x->SetThickness( global_thickness );
  tl2y->SetThickness( global_thickness );
  tl->SetBeam( global_beam, global_spread);   
  tl2->SetBeam( global_beam, global_spread);  
  tl2x->SetBeam( global_beam, global_spread);
  tl2y->SetBeam( global_beam, global_spread);


  AnaTel *tsmile = new AnaTel(telescopebuild);
  tsmile->SetBeam( ebeam );
  tsmile->SetThickness( global_thickness );
  tsmile->SetBeam( global_beam, global_spread);
  tsmile->SetResolution(m26_resolution);

  // predicted resolution, fill pre with posx and prediction_y
  TGraph*  pre[ngraphs];
  Double_t prediction_y[nplanes];
  // = {
  //tl->GetWidth(0,0)*1000.0,
  //tl->GetWidth(1,0)*1000.0,
  //tl->GetWidth(2,0)*1000.0,
  //tl->GetWidth(3,0)*1000.0,
  //tl->GetWidth(4,0)*1000.0
  //};

  // if(verbose)
  //  std::cout << "Resolution estimate: " << tl->GetWidth(0,0) << std::endl;

  // The starting Mimosa26 singlepoint resolution
  Double_t singlepoint_resolution = 0.001;

  // Start predicting
  if (verbose1)
  {
  std::cout << "Calculating Resolution predictions" << std::endl;
  std::cout << " " << std::endl;
  }

  for(Int_t icurve = 0; icurve <  ngraphs; icurve++)
  {
  tl->SetResolution(singlepoint_resolution);
  if(verbose0)
  std::cout << " Resolution prediction on plane : " ;
  for(Int_t j = 0; j < nplanes; j++)
  {
  prediction_y[j] = tl->GetWidth(j,0)*1000.0;
  if(verbose0)
  printf(" %d: %8.3f", j, prediction_y[j] );
  }
  if(verbose0)
  std::cout << std::endl; 
  pre[icurve] = new TGraph(nplanes, posx, prediction_y);
  pre[icurve]->SetLineStyle(icurve+1);
  pre[icurve]->SetLineWidth(1.3);
  singlepoint_resolution += 0.0010; 
}

// hacks
AnaTel *teltest = new AnaTel(telescopebuild);
MyFunctionObject2D *fobjtest = new MyFunctionObject2D(teltest, obsresol_x, obsresol_y, obsresol_error_x, obsresol_error_y );

teltest->SetThickness( global_thickness );
teltest->SetBeam( global_beam, global_spread);  
teltest->SetResolution( m26_resolution);

Double_t tempchi2  = fobjtest->getchi2_2D(m26_resolution );

cout << "the chi2 is " << tempchi2 << std::endl;




cout << " " << std::endl;

// std::cout << " observed " << obsresol_x[0] << std::endl;
// std::cout << " observed error " << obsresol_error_x[0] << std::endl;
cout << " 26 fit " << m26_resolution << std::endl;
cout << " 26 er " << m26_res_error << std::endl;


double the_resi = m26_resolution*1000.0;
double the_error = 0.1;

//m26_res_error*1000.0;


double temptemp = 0;


for (int i = 0; i< 6;i++)
{
  //    temptemp+= (sqrt(obsresol_x[i]*obsresol_x[i]/(1+get_k(i))) - m26_resolution*1000.0)*(sqrt(obsresol_x[i]*obsresol_x[i]/(1+get_k(i))) - m26_resolution*1000.0)/ m26_res_error/1000.0/m26_res_error/1000.0 + (sqrt(obsresol_y[i]*obsresol_y[i]/(1+get_k(i))) - m26_resolution*1000.0)*(sqrt(obsresol_y[i]*obsresol_y[i]/(1+get_k(i))) - m26_resolution*1000.0)/ m26_res_error/1000.0/m26_res_error/1000.0;

  std::cout << " observed " << obsresol_x[i] << std::endl;
  std::cout << " observed error " << obsresol_error_x[i] << std::endl;

  std::cout << " k is " << get_k(i) << std::endl;
  std::cout << " sqrt ob^2/1+k is " << sqrt(obsresol_x[i]*obsresol_x[i]/(1+get_k(i))) << std::endl;

  // temptemp += pow(((sqrt((obsresol_x[i]*obsresol_x[i]-obsresol_error_x[i]*obsresol_error_x[i])/(1+get_k(i))) - the_resi)/the_error) , 2);
  // temptemp += pow(((sqrt((obsresol_y[i]*obsresol_y[i]-obsresol_error_y[i]*obsresol_error_y[i])/(1+get_k(i))) - the_resi)/the_error) , 2);

  temptemp += pow(((sqrt((obsresol_x[i]*obsresol_x[i]+obsresol_error_x[i]*obsresol_error_x[i])/(1+get_k(i))) - the_resi)/the_error) , 2);
  temptemp += pow(((sqrt((obsresol_y[i]*obsresol_y[i]+obsresol_error_y[i]*obsresol_error_y[i])/(1+get_k(i))) - the_resi)/the_error) , 2);

}

cout << std::endl;
cout << " CHI2/ndf " << temptemp/11.0 << std::endl;
cout << std::endl;





// Plot the chi2
TCanvas *canv_chi2 = new TCanvas("chi2","chi2",1400,10,500,500);
canv_chi2->Divide(1,1,0.02,0.02);
//  std::cout << "New resolution estimate: " << tl->GetWidth(0,0) << std::endl;
Double_t chi2_axis_Xmin = 0.0001;
Double_t chi2_axis_Xmax = 0.007;
Double_t epsilon = 0.1;
Int_t    niter = 10;

MyFunctionObject2D *fobj = new MyFunctionObject2D(tl2, obsresol_x, obsresol_y, obsresol_error_x, obsresol_error_y );
MyFunctionObject1D *fobjX = new MyFunctionObject1D(tl2x, obsresol_x, obsresol_error_x );
MyFunctionObject1D *fobjY = new MyFunctionObject1D(tl2y, obsresol_y, obsresol_error_y );

TF1* fchi2      = new TF1("fchi2",      fobj,  chi2_axis_Xmin, chi2_axis_Xmax, 1, "MyFunctionObject2D");
TF1* fchi2X     = new TF1("fchi2X",     fobjX, chi2_axis_Xmin, chi2_axis_Xmax, 1, "MyFunctionObject1D");
TF1* fchi2Y     = new TF1("fchi2Y",     fobjY, chi2_axis_Xmin, chi2_axis_Xmax, 1, "MyFunctionObject1D");


fchi2->SetLineColor(45);
fchi2->SetNpx(100);
fchi2->Draw("L");

fchi2X->SetLineColor(25);
fchi2X->SetNpx(100);
fchi2X->Draw("LSAME");

fchi2Y->SetLineColor(65);
fchi2Y->SetNpx(100);
fchi2Y->Draw("LSAME");

Double_t resolmin2D      = fchi2->GetMinimumX(chi2_axis_Xmin, chi2_axis_Xmax, epsilon, niter);
Double_t resolmin1DX     = fchi2X->GetMinimumX(chi2_axis_Xmin, chi2_axis_Xmax, epsilon, niter);
Double_t resolmin1DY     = fchi2Y->GetMinimumX(chi2_axis_Xmin, chi2_axis_Xmax, epsilon, niter);

Double_t chimin2D  = fobj->getchi2_2D(resolmin2D );
Double_t chimin1D_X = fobjX->getchi2_1D(resolmin1DX);
Double_t chimin1D_Y = fobjY->getchi2_1D(resolmin1DY);

Double_t  ymin = 0.;
Double_t  ymax = chimin2D*2.;

fchi2->SetMinimum(ymin);
fchi2X->SetMinimum(ymin);
fchi2Y->SetMinimum(ymin);

fchi2->SetMaximum(ymax);
fchi2X->SetMaximum(ymax);
fchi2Y->SetMaximum(ymax);

if(verbose)
  std::cout << "Minimal 2D resolution: " << resolmin2D << " . Minimal 2D Chi2: " << chimin2D << std::endl;

  std::cout << "hackchi2 " << std::endl;


  // starting values only:
  Double_t resol_p =  resolmin2D;
  Double_t resol_n = -resolmin2D/3.;
  Double_t resolX_p =  resolmin1DX;
  Double_t resolX_n = -resolmin1DX/3.;
  Double_t resolY_p =  resolmin1DY;
  Double_t resolY_n = -resolmin1DY/3.;

  Double_t value_x_p = fchi2->GetX( chimin2D+1., resolmin2D, resolmin2D+5., epsilon, niter  );
  //   Double_t value_x_n = fobj->GetX( chimin2D+1.);
  //cout << value_x_p <<  std::endl;

  Double_t fraction_step = TMath::Abs( value_x_p-resolmin2D ) /10.;

  Int_t nstep     = 100; // step to find the plus and minus range of the resolmin error bars
  Double_t dstep  = 0.01  ; // delta = 1/nstep and dstep/nstep sets the "resolution" of the error bars scan.
if(fraction_step < dstep/nstep )
  fraction_step = dstep/nstep;

  for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  double valuep= fobj->getchi2_2D(resolmin2D + delta );
  if( valuep>= chimin2D+1.0)
  {
    resol_p = + delta ;
    if(verbose)
      printf("p %d %8.3f %8.3f \n", i, resol_p,  valuep);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  if(resolmin2D-delta<1e-7)
    break;
  double valuen= fobj->getchi2_2D(resolmin2D - delta );
  if( valuen>= chimin2D+1.0)
  {
    resol_n = -delta;
    if(verbose)
      printf("n %d %8.3f %8.3f \n", i, resol_n,  valuen);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  double valuep= fobjX->getchi2_1D(resolmin1DX + delta );
  if( valuep>= chimin1D_X+1.0)
  {
    resolX_p = + delta ;
    if(verbose)
      printf("p %d %8.3f %8.3f \n", i, resolX_p,  valuep);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  if(resolmin1DX-delta<1e-7)
    break;
  double valuen= fobjX->getchi2_1D(resolmin1DX-delta );
  if( valuen>= chimin1D_X+1.0)
  {
    resolX_n = -delta;
    if(verbose)
      printf("n %d %8.3f %8.3f \n", i, resolX_n,  valuen);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  double valuep= fobjY->getchi2_1D(resolmin1DY + delta );
  if( valuep>= chimin1D_Y+1.0)
  {
    resolY_p = + delta ;
    if(verbose)
      printf("p %d %8.3f %8.3f \n", i, resolY_p,  valuep);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  if(resolmin1DY-delta<1e-7) break;
  double valuen= fobjY->getchi2_1D(resolmin1DY-delta );
  if( valuen>= chimin1D_Y+1.0)
  {
    resolY_n = -delta;
    if(verbose)
      printf("n %d %8.3f %8.3f \n", i, resolY_n,  valuen);
    break;
  }
}


TGraphErrors *x_direction = new TGraphErrors(nplanes, posx, obsresol_x, posx_error, obsresol_error_x);
TGraphErrors *y_direction = new TGraphErrors(nplanes, posx, obsresol_y, posx_error, obsresol_error_y);

resolmin2D=m26_resolution;

tl->SetResolution(resolmin2D);

// Get errors for t1
Double_t prediction_error[6];

// for the plot
float average_error = 0.0;


for(Int_t j = 0; j < nplanes; j++)
{
  tl->SetResolution(resolmin2D);
  prediction_y[j] = tl->GetWidth(j,0)*1000.0;
  tl->SetResolution(resolmin2D-0.0001);
  const Double_t up = TMath::Abs(tl->GetWidth(j,0)*1000.0 - prediction_y[j]);
  tl->SetResolution(resolmin2D+0.0001);
  const Double_t down = TMath::Abs(tl->GetWidth(j,0)*1000.0 - prediction_y[j]);
  prediction_error[j] = 0.5*(up+down);
  average_error += prediction_error[j];

}

for (int j = 0; j<6;j++)
cout << "asdfasdfasdfasdf " << get_k(j) << std::endl;


// fail!
// Errorband for the smilie plot
cout << "asdfasdf " << tsmile->GetWidth(0,0)*1000.0 << std::endl;
cout << "asdfasdf " << tsmile->GetWidth(1,0)*1000.0 << std::endl;
cout << "asdfasdf " << tsmile->GetWidth(2,0)*1000.0 << std::endl;
cout << "asdfasdf " << tsmile->GetWidth(3,0)*1000.0 << std::endl;


// fix error
for (int j=0;j<6;j++)
{
  //prediction_error[j] = prediction_error[j] / (sqrt(12.0));



  std::cout << std::endl;
  std::cout << " !!! " << std::endl;
  std::cout << std::endl;
  std::cout << "Prediction error on plane " << j << " is " << prediction_error[j] << std::endl;
  std::cout << std::endl;

  float tempup = 0.0;
  float tempdn = 0.0;

  tsmile->SetResolution(m26_resolution);

  float deltam26 = 0.1;


  double localscatter = sqrt(tsmile->GetWidth(j,0)*1000.0*tsmile->GetWidth(j,0)*1000.0 - m26_resolution*1000.0*m26_resolution*1000.0*(1 + get_k(j)));
  float deltams = 0.1*localscatter;
  std::cout << "scatter here " << localscatter << std::endl;

  float thedelta = sqrt(pow(((1+get_k(j))*m26_resolution*1000.0*deltam26/(tsmile->GetWidth(j,0)*1000.0)),2) + pow((localscatter*deltams/(tsmile->GetWidth(j,0)*1000.0)),2));

  std::cout << " THE DELTA " << thedelta << std::endl;
  std::cout << "error times k upper " << sqrt((m26_resolution*1000.0+prediction_error[j])*(m26_resolution*1000.0+prediction_error[j])*(1 + get_k(j)) + localscatter*localscatter) << std::endl;
  std::cout << "error times k lower " << sqrt((m26_resolution*1000.0-prediction_error[j])*(m26_resolution*1000.0-prediction_error[j])*(1 + get_k(j)) + localscatter*localscatter) << std::endl;



  tsmile->SetResolution(m26_resolution+prediction_error[j]/1000.0);
  tempup = tsmile->GetWidth(j,0)*1000.0;
  std::cout << " max " << m26_resolution*1000.0+prediction_error[j] << "    " << tsmile->GetWidth(j,0)*1000.0 << std::endl;
  tsmile->SetResolution(m26_resolution-prediction_error[j]/1000.0);
  std::cout << " max " << m26_resolution*1000.0-prediction_error[j] << "    " << tsmile->GetWidth(j,0)*1000.0 << std::endl;
  tempdn = tsmile->GetWidth(j,0)*1000.0;


  prediction_error[j] = (tempup - tempdn)/2.0;
  std::cout << "Prediction error on plane " << j << " is " << prediction_error[j] <<  " in sig meas! " << std::endl;
  std::cout << std::endl;


  prediction_error[j] = thedelta;

}




average_error = average_error / nplanes;

global_plot_error = average_error;


cout << "global_plot_error " << global_plot_error << std::endl;

global_plot_error = global_plot_error/ (sqrt(12.0));

cout << "global_plot_error / sqrt 12 " << global_plot_error << std::endl;


// Get fiterrors on t1 for plane 2
tl->SetResolution(resolmin2D);
Double_t fiterror_central =  tl->GetPointingRes(2,0);
tl->SetResolution(resolmin2D+0.0001);
Double_t fiterror_up = TMath::Abs(tl->GetPointingRes(2,0) - fiterror_central);
if(verbose)
  std::cout << "fiterror_up " << fiterror_up  << std::endl;
  tl->SetResolution(resolmin2D-0.0001);
  Double_t fiterror_down = TMath::Abs(tl->GetPointingRes(2,0) - fiterror_central);
if(verbose)
  std::cout << "fiterror_down " << fiterror_down  << std::endl;
  Double_t fiterror_error = 0.5*(fiterror_up + fiterror_down);

  canv_chi2->Print("pics/"+outputname+"_chi2.eps");
  _outputFile->cd();
  canv_chi2->Write();
  canv_chi2->Close();


  //TGraphErrors *sys_prediction = new TGraphErrors(nplanes, posx, prediction_y, posx_error, prediction_error);
  TGraphErrors *sys_prediction = new TGraphErrors(nplanes, posx, prediction_y);

  sys_prediction->SetLineWidth(3);

  TGraphErrors *sys_prediction2 = new TGraphErrors(nplanes, posx, prediction_y, posx_error, prediction_error);
  sys_prediction2->SetFillColor(kOrange);

  //


  TH1D *h_axis;

if (planedistance == 150)
{
  h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 150 mm",1, -75, 825.0);
  //h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 150 mm",100, -75, 825.0*2);
  std::cout << "planedistance 150" << std::endl;
} else {
  h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 20 mm",1, -10, 110.0);
  //h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 20 mm",100, -10, 110.0*2);
  std::cout << "planedistance 20" << std::endl;
}
TCanvas *smilie = new TCanvas("smilie","smilie",10,10,800,600);
smilie->SetRightMargin(0.3);
// Do smilie plot, x/y_direction has the measured residuals, pre the predictions and sys_prediction the errorband around the fit
//  gStyle->SetPadBorderMode(0);
//  gStyle->SetOptStat(0);
smilie->SetFillColor(0);
//  smilie->Divide(1,1);
//  smilie->cd(1);
gStyle->SetErrorX(0);

gPad->SetLogx(0);
gPad->SetLogy(0);
x_direction->SetMarkerStyle(21);
x_direction->SetMarkerColor(kRed);
x_direction->SetMarkerSize(2.5);
y_direction->SetMarkerStyle(22);
y_direction->SetMarkerColor(kBlue);
y_direction->SetMarkerSize(2.5);
h_axis->GetXaxis()->SetLabelFont(42);
h_axis->GetXaxis()->SetLabelSize(0.035);
h_axis->GetXaxis()->SetTitleSize(0.035);
h_axis->GetXaxis()->SetTitleFont(42);
h_axis->GetXaxis()->SetTitle("z in mm");
h_axis->GetYaxis()->SetLabelFont(42);
h_axis->GetYaxis()->SetLabelSize(0.035);
h_axis->GetYaxis()->SetTitleSize(0.035);
h_axis->GetYaxis()->SetTitleFont(42);
h_axis->GetYaxis()->SetTitle("#sigma_{meas} in #mum");
//histo_cfg(h_axis, "z in mm","#sigma_{meas.} in #mum","");
h_axis->SetMinimum(0.0);

// The y-axis can be adjusted to geometry too
if(telescopebuild == "AnaTel_wide.geom")
h_axis->SetMaximum(30.);
if(telescopebuild == "AnaTel_narrow.geom")

if (planedistance == 20)
  h_axis->SetMaximum(15.);
if (planedistance == 150)
  h_axis->SetMaximum(30.);

  h_axis->Draw("hist");

  sys_prediction2->Draw("LE3 same");
  sys_prediction->Draw("L same");
  x_direction->Draw("p same");
  y_direction->Draw("p same");

  // Draw predictions
  //for(Int_t i = 0; i < ngraphs; i++)
  for(Int_t i = 1; i < 5; i++)
  pre[i]->Draw("L same");

  // Legend for x and y
  TLegend *leg = new TLegend(0.19,0.70,0.50,0.80);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(x_direction,"x direction","p");
  leg->AddEntry(y_direction,"y direction","p");
  leg->SetTextSize(0.03);
  leg->Draw();

  // Create second legend with the predictions
  TLegend *leg2 = new TLegend(0.7,0.6,0.95,0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  //leg2->SetHeader("Prediction for intrinsic resolution:");
  char clabel[100]="";
  sprintf(clabel,  "#sigma_{M26} = (%4.2f #pm %4.2f) #mum", resolmin2D*1000., global_plot_error);
  //sprintf(clabel,  "#sigma_{DUT} = (%4.2f #pm %4.2f) #mum", resolmin2D*1000., m26_res_error*1000.);
  leg2->AddEntry(sys_prediction, clabel ,"l");
  //leg2->AddEntry(pre[0],"#sigma_{M26} = 1.0 #mum","l");
  leg2->AddEntry("","Prediction for intrinsic resolution:","h");
  leg2->AddEntry(pre[1],"#sigma_{M26} = 2.0 #mum","l");
  leg2->AddEntry(pre[2],"#sigma_{M26} = 3.0 #mum","l");
  leg2->AddEntry(pre[3],"#sigma_{M26} = 4.0 #mum","l");
  leg2->AddEntry(pre[4],"#sigma_{M26} = 5.0 #mum","l");
  //leg2->AddEntry(pre[5],"#sigma_{M26} = 6.0 #mum","l");
  leg2->SetTextSize(0.03);
  leg2->Draw();

  // Create third legend with some additional information for the plot
  TLegend *leg3 = new TLegend(0.7,0.10,0.95,0.30);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->SetFillStyle(0);
  char beambuffer[100]="";
  sprintf(beambuffer, "DESY %d GeV e^{-}", nominalbeam);
  leg3->AddEntry((TObject*)0, beambuffer, "");
  TString threshlabel = "M26 threshold: ";
  threshlabel += global_thresh;
  leg3->AddEntry((TObject*)0, threshlabel, "");
  char chi2buffer[100]="";
  sprintf(chi2buffer, "The #Chi^{2} of the fit is %4.2f", (chimin2D/12.0));
  //leg3->AddEntry((TObject*)0, chi2buffer, "");
  leg3->AddEntry(x_direction,"x direction","p");
  leg3->AddEntry(y_direction,"y direction","p");
  leg3->SetTextSize(0.030);
  leg3->Draw();

  // Output
  smilie->Print("pics/"+outputname+"_smilie.eps");
  _outputFile->cd();
  smilie->Write();
  smilie->Close();
  */

  }

// Fitting of each file -> noise and efficiency
void noise(Int_t runnumber)
{
  std::string help;
  if(runnumber < 113 || runnumber > 700) help = inputDir150 + inputFile;
  else help = inputDir20 + inputFile;
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




  /*
     Double_t noisevalue = 0;
     Double_t noisevalue_error = 0;
     Double_t efficiency = 0;
     Double_t efficiency_error = 0;

     noisevalue = h_m26[4]->GetMean(2);
     noisevalue_error = h_m26[4]->GetMeanError(2);
     efficiency = h_m26_eff[4]->GetMean(2);
     efficiency_error = h_m26_eff[4]->GetMeanError(2);

     avgnoise = noisevalue;
     avgnoise_error = noisevalue_error;
     avgeffi = efficiency;
     avgeffi_error = efficiency_error;
     */

  std::cout << "   - run done" << std::endl;

}


// FIXME update code to get clustersize from hits in GBL tracks
void getclusize(Int_t runnumber)
{
  std::string help;
  if(runnumber < 113 || runnumber > 700) help = inputDir150 + inputFile;
  else help = inputDir20 + inputFile;
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

  // Load file
  h_m26[0] = (TH1D*) (tfile->Get("DUTHisto00/clusterSizeX"))->Clone();
  h_m26[1] = (TH1D*) (tfile->Get("DUTHisto00/clusterSizeY"))->Clone();
  h_m26[2] = (TH1D*) (tfile->Get("DUTHisto01/clusterSizeX"))->Clone();
  h_m26[3] = (TH1D*) (tfile->Get("DUTHisto01/clusterSizeY"))->Clone();
  h_m26[4] = (TH1D*) (tfile->Get("DUTHisto02/clusterSizeX"))->Clone();
  h_m26[5] = (TH1D*) (tfile->Get("DUTHisto02/clusterSizeY"))->Clone();
  h_m26[6] = (TH1D*) (tfile->Get("DUTHisto03/clusterSizeX"))->Clone();
  h_m26[7] = (TH1D*) (tfile->Get("DUTHisto03/clusterSizeY"))->Clone();
  h_m26[8] = (TH1D*) (tfile->Get("DUTHisto04/clusterSizeX"))->Clone();
  h_m26[9] = (TH1D*) (tfile->Get("DUTHisto04/clusterSizeY"))->Clone();
  h_m26[10] = (TH1D*) (tfile->Get("DUTHisto05/clusterSizeX"))->Clone();
  h_m26[11] = (TH1D*) (tfile->Get("DUTHisto05/clusterSizeY"))->Clone();

  // Add the values for all 6 sensor planes, in x and y and divide by 12
  Double_t clustervalue = 0;
  Double_t clustervalue_error = 0;
  for(int i=0;i<planescount;i++)
  {
    clustervalue += h_m26[i]->GetMean(1);
    clustervalue_error += h_m26[i]->GetMeanError(1);
  }
  avgclustersize = clustervalue / (float)planescount;
  avgclustersize_error = clustervalue_error / (float)planescount;
}


void getpointing(Int_t runnumber, float sigm26)
{
  std::string help;
  if(runnumber < 113 || runnumber > 700) help = inputDir150 + inputFile;
  else help = inputDir20 + inputFile;
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

  std::cout << "File is " << runnumber << std::endl;

  // Load file
  h_m26[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftX"))->Clone();
  h_m26[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftY"))->Clone();
  h_m26[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftX"))->Clone();
  h_m26[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftY"))->Clone();
  h_m26[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftX"))->Clone();
  h_m26[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftY"))->Clone();
  h_m26[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX"))->Clone();
  h_m26[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY"))->Clone();
  h_m26[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftX"))->Clone();
  h_m26[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftY"))->Clone();
  h_m26[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftX"))->Clone();
  h_m26[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftY"))->Clone();

  Double_t sigmeas = 0;
  Double_t sigmeas_error = 0;
  for (int i = 4; i<6; i++)
  {
    TF1 *f1 = new TF1("f1","gaus");
    TF1 *f2 = new TF1("f2","gaus");
    f2->SetLineWidth(2);
    f2->SetLineStyle(1);
    f2->SetLineColor(kBlack);
    h_m26[i]->Fit(f1,"EMQI0","");
    Double_t mean_1 = f1->GetParameter(1);
    Double_t sigma_1 = f1->GetParameter(2);



    // Repeat within 2 sigma
    h_m26[i]->Fit(f2,"EQMI","", (mean_1 - 2.0*sigma_1), (mean_1 + 2.0*sigma_1));
    Double_t mean_2 = f2->GetParameter(1);
    Double_t sigma_2 = f2->GetParameter(2);
    Double_t sigma_2e = f2->GetParError(2);

    std::cout << "measured width: " << sigma_2*1000.0 << std::endl;

    sigmeas += sigma_2*1000.0;
    sigmeas_error += sigma_2e*1000.0;
  }


  //sigm26 = 3.42;
  //Double_t sigm26_e = 0.12;
  Double_t sigm26_e = 0.035;
  /*
     if (runnumber == 291)
     {
     sigm26 = 3.55;
     std::cout << "asdfasdf" << std::endl; 
     }

     if (runnumber == 755)
     {
     sigm26 = 3.18;
     std::cout << "asdf" << std::endl;  

     }*/

  // These lines for pointing at plane 3:

  /*
     avgmeas = sqrt( (sigmeas / 2.0)*(sigmeas / 2.0) - sigm26*sigm26);
  //avgmeas_error = sqrt( (sigmeas_error / 2.0)*(sigmeas_error / 2.0) + sigm26_e*sigm26_e);
  avgmeas_error = sqrt( (sigmeas_error / 2.0 * sigmeas / 2.0 / avgmeas)*(sigmeas_error / 2.0 * sigmeas / 2.0 / avgmeas) + (sigm26*sigm26_e/avgmeas)*(sigm26*sigm26_e/avgmeas)     );

  std::cout << "sigmeas is " << sigmeas/2.0 << " pm " << sigmeas_error/2.0 << std::endl;
  std::cout << "pointing res at plane 3 is " << avgmeas << " pm  " << avgmeas_error << std::endl;
  std::cout << std::endl;
  */

  // these lines for extrapolation to telescope center:


  float kfive = 0.2209302326;

  avgmeas = sqrt( (sigmeas / 2.0)*(sigmeas / 2.0) + sigm26*sigm26*((1.0/6.0) - kfive - 1.0) );
  avgmeas_error = sqrt((sigmeas/2.0*sigmeas_error/2.0/avgmeas)*(sigmeas/2.0*sigmeas_error/2.0/avgmeas) + (sigm26*sigm26_e/avgmeas*((1.0/6.0) - kfive - 1.0))*(sigm26*sigm26_e/avgmeas*((1.0/6.0) - kfive - 1.0)));

  std::cout << "sigmeas is " << sigmeas/2.0 << " pm " << sigmeas_error/2.0 << std::endl;
  std::cout << "pointing res at center telescope is " << avgmeas << " pm  " << avgmeas_error << std::endl;
  std::cout << std::endl;


}


void histoplot(Int_t runnumber)
{
  std::string help;
  if(runnumber < 113 || runnumber > 700) help = inputDir150 + inputFile;
  else help = inputDir20 + inputFile;
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

  std::cout << "File is " << runnumber << std::endl;

  // Load file
  h_m26[0] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftX"))->Clone();
  h_m26[1] = (TH1D*) (tfile->Get("DUTHisto00/DUTshiftY"))->Clone();
  h_m26[2] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftX"))->Clone();
  h_m26[3] = (TH1D*) (tfile->Get("DUTHisto01/DUTshiftY"))->Clone();
  h_m26[4] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftX"))->Clone();
  h_m26[5] = (TH1D*) (tfile->Get("DUTHisto02/DUTshiftY"))->Clone();
  h_m26[6] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX"))->Clone();
  h_m26[7] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY"))->Clone();
  h_m26[8] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftX"))->Clone();
  h_m26[9] = (TH1D*) (tfile->Get("DUTHisto04/DUTshiftY"))->Clone();
  h_m26[10] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftX"))->Clone();
  h_m26[11] = (TH1D*) (tfile->Get("DUTHisto05/DUTshiftY"))->Clone();

  double val[planescount];


  for (int i = 0; i<planescount; i++)
  {
    TF1 *f1 = new TF1("f1","gaus");
    TF1 *f2 = new TF1("f2","gaus");
    f2->SetLineWidth(2);
    f2->SetLineStyle(1);
    f2->SetLineColor(kBlack);
    h_m26[i]->Fit(f1,"EMQI0","");
    Double_t mean_1 = f1->GetParameter(1);
    Double_t sigma_1 = f1->GetParameter(2);



    // Repeat within 2 sigma
    h_m26[i]->Fit(f2,"EQMI","", (mean_1 - 2.0*sigma_1), (mean_1 + 2.0*sigma_1));
    Double_t mean_2 = f2->GetParameter(1);
    Double_t sigma_2 = f2->GetParameter(2);
    Double_t sigma_2e = f2->GetParError(2);

    std::cout << "measured width: " << sigma_2*1000.0 << std::endl;

    val[i] = sigma_2*1000.0;

    /*
       if (i==0 || i == 1 || i == 10 || i == 11)
       {
       delta0->Fill(sigma_2*1000.0);
       }
       if (i==2 || i == 3 || i == 8 || i == 9)
       {
       delta1->Fill(sigma_2*1000.0);
       }
       if (i==4 || i == 5 || i == 6 || i == 7)
       {
       delta2->Fill(sigma_2*1000.0);
       }*/

  }

  delta0->Fill(fabs(val[0]-val[1]));
  delta1->Fill(fabs(val[2]-val[3]));
  delta2->Fill(fabs(val[4]-val[5]));
  delta3->Fill(fabs(val[6]-val[7]));
  delta4->Fill(fabs(val[8]-val[9]));
  delta5->Fill(fabs(val[10]-val[11]));





  delta6->Fill(fabs(val[0]-val[10]));
  delta7->Fill(fabs(val[2]-val[8]));
  delta8->Fill(fabs(val[4]-val[6]));
  delta9->Fill(fabs(val[1]-val[11]));
  delta10->Fill(fabs(val[3]-val[9]));
  delta11->Fill(fabs(val[5]-val[7]));





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
    name = "Mean" + std::to_string(i);
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
    name = "Measured-Residual" + std::to_string(i);
    g_sigma[i]->SetNameTitle(name.c_str(), name.c_str());
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
    g_mean[i]->Draw("ap");
    g_mean[i]->Write();
  }


  for(int i = 0; i<2*nplanes; i++)
  {
    g_sigma[i]->Draw("ap");
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
int main()
{

  //gROOT->SetBatch();
  //gSystem->Load("libMinuit");
  //gSystem->Load("lib/libGBL.so");
  //gSystem->AddIncludePath("include/");
  gROOT->SetStyle("Plain");

  inputDir20  = "../../analysis-20mm/output/histograms/";
  inputDir150  = "../../analysis-150mm/output/histograms/";
  inputFile = "run00";

  _outputFile = new TFile("../bin/output.root", "RECREATE");
  _outputFile->cd();

  // Initial telescope constructor name, AnaTel_wide.geom for 150mm, AnaTel_narrow.geom for 20mm data
  //char telescopebuild[50];
  //int tempint = 0;

  //telescopebuild = "AnaTel_wide.geom";


  // Initialises the run vectors
  std::vector<int> run_6_150;
  std::vector<int> run_5_150;
  std::vector<int> run_4_150;
  std::vector<int> run_3_150;
  std::vector<int> run_2_150;
  std::vector<int> run_6_20;
  std::vector<int> run_5_20;
  std::vector<int> run_4_20;
  std::vector<int> run_3_20;
  std::vector<int> run_2_20;
  std::vector<int> run120;

  /*std::vector<int> testing;

    testing.push_back(291);
    testing.push_back(294);
    testing.push_back(295);
    testing.push_back(296);
    testing.push_back(297);
    testing.push_back(298);
    testing.push_back(299);
    testing.push_back(300);
    */

  // Fills the run vectors with runnumbers
  if(true)
  {
    // 150 mm
    // energy 6 GeV:
    run_6_150.push_back(60);	// thr 3
    run_6_150.push_back(61);	// thr 4
    run_6_150.push_back(62);	// thr 5
    run_6_150.push_back(63);	// thr 6
    run_6_150.push_back(64);	// thr 7
    run_6_150.push_back(65);	// thr 8
    run_6_150.push_back(66);	// thr 9
    run_6_150.push_back(67);	// thr 10
    run_6_150.push_back(68);	// thr 11
    run_6_150.push_back(69);	// thr 12

    // energy 5 GeV:
    run_5_150.push_back(70);	// thr 3
    run_5_150.push_back(71);	// thr 4
    run_5_150.push_back(72);	// thr 5
    run_5_150.push_back(73);	// thr 6
    run_5_150.push_back(74);	// thr 7
    run_5_150.push_back(75);	// thr 8
    run_5_150.push_back(76);	// thr 9
    run_5_150.push_back(77);	// thr 10
    run_5_150.push_back(78);	// thr 11
    run_5_150.push_back(79);	// thr 12

    // energy 4 GeV:
    run_4_150.push_back(81);	// thr 3
    run_4_150.push_back(82);	// thr 4
    run_4_150.push_back(83);	// thr 5
    run_4_150.push_back(84);	// thr 6
    run_4_150.push_back(85);	// thr 7
    run_4_150.push_back(86);	// thr 8
    run_4_150.push_back(89);	// thr 9
    run_4_150.push_back(90);	// thr 10
    run_4_150.push_back(91);	// thr 11
    run_4_150.push_back(92);	// thr 12

    // energy 3 GeV:
    run_3_150.push_back(93);	// thr 3
    run_3_150.push_back(94);	// thr 4
    run_3_150.push_back(95);	// thr 5
    run_3_150.push_back(96);	// thr 6
    run_3_150.push_back(97);	// thr 7
    run_3_150.push_back(98);	// thr 8
    run_3_150.push_back(99);	// thr 9
    run_3_150.push_back(100);	// thr 10
    run_3_150.push_back(101);	// thr 11
    run_3_150.push_back(102);	// thr 12

    // energy 2 GeV:
    run_2_150.push_back(103);	// thr 3
    run_2_150.push_back(104);	// thr 4
    run_2_150.push_back(105);	// thr 5
    run_2_150.push_back(106);	// thr 6
    run_2_150.push_back(107);	// thr 7
    run_2_150.push_back(108);	// thr 8
    run_2_150.push_back(109);	// thr 9
    run_2_150.push_back(110);	// thr 10
    run_2_150.push_back(111);	// thr 11
    run_2_150.push_back(112);	// thr 12

    // 20 mm
    // energy 6 GeV:
    run_6_20.push_back(114);	// thr 3
    run_6_20.push_back(115);	// thr 4
    run_6_20.push_back(116);	// thr 5
    run_6_20.push_back(117);	// thr 6
    run_6_20.push_back(118);	// thr 7
    run_6_20.push_back(119);	// thr 8
    run_6_20.push_back(120);	// thr 9
    run_6_20.push_back(121);	// thr 10
    run_6_20.push_back(122);	// thr 11
    run_6_20.push_back(123);	// thr 12

    // energy 5 GeV:
    run_5_20.push_back(124);	// thr 3
    run_5_20.push_back(125);	// thr 4
    run_5_20.push_back(126);	// thr 5
    run_5_20.push_back(127);	// thr 6
    run_5_20.push_back(128);	// thr 7
    run_5_20.push_back(129);	// thr 8
    run_5_20.push_back(130);	// thr 9
    run_5_20.push_back(131);	// thr 10
    run_5_20.push_back(132);	// thr 11
    run_5_20.push_back(133);	// thr 12

    // energy 4 GeV:
    run_4_20.push_back(134);	// thr 3
    run_4_20.push_back(135);	// thr 4
    run_4_20.push_back(136);	// thr 5
    run_4_20.push_back(137);	// thr 6
    run_4_20.push_back(140);	// thr 7
    run_4_20.push_back(141);	// thr 8
    run_4_20.push_back(142);	// thr 9
    run_4_20.push_back(143);	// thr 10
    run_4_20.push_back(144);	// thr 11
    run_4_20.push_back(145);	// thr 12

    // energy 3 GeV:
    run_3_20.push_back(146);	// thr 3
    run_3_20.push_back(147);	// thr 4
    run_3_20.push_back(148);	// thr 5
    run_3_20.push_back(153);	// thr 6
    run_3_20.push_back(154);	// thr 7
    run_3_20.push_back(155);	// thr 8
    run_3_20.push_back(156);	// thr 9
    run_3_20.push_back(157);	// thr 10
    run_3_20.push_back(158);	// thr 11
    run_3_20.push_back(161);	// thr 12

    // energy 2 GeV:
    run_2_20.push_back(162);	// thr 3
    run_2_20.push_back(163);	// thr 4
    run_2_20.push_back(164);	// thr 5
    run_2_20.push_back(165);	// thr 6
    run_2_20.push_back(166);	// thr 7
    run_2_20.push_back(167);	// thr 8
    run_2_20.push_back(168);	// thr 9
    run_2_20.push_back(169);	// thr 10
    run_2_20.push_back(170);	// thr 11
    run_2_20.push_back(171);	// thr 12


    // energy 120 GeV CERN, 150mm Data:
    //run120.push_back(752);	// thr 3
    //run120.push_back(753);	// thr 4
    run120.push_back(754);	// thr 5 // !!?? Push back 754 as place holders
    run120.push_back(754);	// thr 5
    run120.push_back(754);	// thr 5
    run120.push_back(755);	// thr 6
    run120.push_back(756);	// thr 7
    run120.push_back(757);	// thr 8
    run120.push_back(758);	// thr 9
    run120.push_back(759);	// thr 10        alternative: 760
    run120.push_back(760);	// thr 10        alternative: 759
    run120.push_back(761);	// thr 11
    run120.push_back(762);	// thr 12

  }

  // Runmode: 0 for all, 1 for clustersize, 2 for threshold, 3 for threshold and clustersize, 4 for noise, 5 for E plot, 9 for testing...
  // 6 for clustersie only, 7 for pointing, 10 for testing of GBLwidths, 11 optimising thickness and Highlanf factor

  Int_t runmode;

  // Open config file
  std::string s_config = "../conf/config.txt";
  ifstream conf(s_config.c_str());
  if(!conf)
  {
    std::cout << "../conf/config.txt file missing!" << std::endl;
    return 1;
  } else
  {
    conf >> runmode;
  }

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

    /*tex = new TLatex(0.0007590689,2.837255,"#circ");
      tex->SetTextColor(4);
      tex->SetTextSize(0.1);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(0.001444937,2.53555,"#circ");
      tex->SetTextColor(2);
      tex->SetTextSize(0.1);
      tex->SetLineWidth(2);
      tex->Draw(); 
      tex = new TLatex(0.004529811,2.13057,"#circ");
      tex->SetTextColor(1);
      tex->SetTextSize(0.1);
      tex->SetLineWidth(2);
      tex->Draw(); */

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
	  if(i==3) std::cout << "j = " << j < " :   ";


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
    for (Int_t j=0; j<threshcount; j++)
    {
      std::cout << " j = " << j << std::endl;


      thres.push_back(j+(13-threshcount));
      ethres.push_back(.0);

      effi_reader( run_2_150[j], 2.0, j ); 
      effi_reader( run_3_150[j], 3.0, j );
      effi_reader( run_4_150[j], 4.0, j );
      effi_reader( run_5_150[j], 5.0, j );
      effi_reader( run_6_150[j], 6.0, j );

      effi_reader( run_2_20[j], 2.0, j );
      effi_reader( run_3_20[j], 3.0, j );
      effi_reader( run_4_20[j], 4.0, j );
      effi_reader( run_5_20[j], 5.0, j );
      effi_reader( run_6_20[j], 6.0, j );

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
      h_effi_E_150[i]= new TGraphErrors(threshcount);
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
      h_effi_E_20[i]= new TGraphErrors(threshcount);
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

    /*TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,-13,113);
      Graph_Graph1->SetMinimum(3.643291);
      Graph_Graph1->SetMaximum(5.702298);
      Graph_Graph1->SetDirectory(0);
      Graph_Graph1->SetStats(0);
      res_meas_XY->SetHistogram(Graph_Graph1);*/

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


  if (runmode == 15)
  {
    std::cout << " This calculates the systematic uncertainties of the pointing resolution for certain assumed input uncerts." << std::endl;

    double ebeam = 6.;
    double kappa = .75;
    // Create telescope
    //
    telescopebuild = "AnaTel_wide.geom";
    AnaTel *tscope150 = new AnaTel(telescopebuild.c_str(), ebeam);
    double sigma_int = 3.24e-3;
    tscope150->SetResolution(sigma_int);
    tscope150->SetKappa(kappa);
    std::cout << " - Tscope created, dz = 150 mm, NO DUT" << std::endl;

    int plane = 3;
    int Biased = 1;


    double sig_0 = tscope150->GetPointingResGBL(plane,Biased);

    tscope150->SetBeam(ebeam*0.95);
    double sig_E_minus = tscope150->GetPointingResGBL(plane,Biased);

    tscope150->SetBeam(ebeam*1.05);
    double sig_E_plus = tscope150->GetPointingResGBL(plane,Biased);

    tscope150->SetBeam(ebeam);

    tscope150->SetKappa(0.9*kappa);
    double sig_kappa_minus = tscope150->GetPointingResGBL(plane,Biased);

    tscope150->SetKappa(1.1*kappa);
    double sig_kappa_plus = tscope150->GetPointingResGBL(plane,Biased);

    tscope150->SetKappa(kappa);

    std::cout << " sys. uncerts in E: sig_0 = " << sig_0*1e3 << " + (" << (sig_E_plus - sig_0)*1e3 << ") - (" << (sig_E_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " sig_E_plus = " << (sig_E_plus)*1e3 << "   sig_E_minus = " << sig_E_minus*1e3 << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_E_plus - sig_0) + fabs(sig_E_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;

    std::cout << " sys. uncerts in KAPPA: sig_0 = " << sig_0*1e3 << " + (" << (sig_kappa_plus - sig_0)*1e3 << ") - (" << fabs(sig_kappa_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " sig_kappa_plus = " << (sig_kappa_plus)*1e3 << "   sig_kappa_minus = " << sig_kappa_minus*1e3 << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_kappa_plus - sig_0) + fabs(sig_kappa_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;

    telescopebuild = "AnaTel_narrow.geom";
    AnaTel *tscope20 = new AnaTel(telescopebuild.c_str(), ebeam);
    tscope20->SetResolution(sigma_int);
    tscope20->SetKappa(kappa);
    std::cout << "\n - Tscope created, dz = 20 mm, NO DUT" << std::endl;

    sig_0 = tscope20->GetPointingResGBL(plane,Biased);

    tscope20->SetBeam(ebeam*0.95);
    sig_E_minus = tscope20->GetPointingResGBL(plane,Biased);

    tscope20->SetBeam(ebeam*1.05);
    sig_E_plus = tscope20->GetPointingResGBL(plane,Biased);

    tscope20->SetBeam(ebeam);

    tscope20->SetKappa(0.9*kappa);
    sig_kappa_minus = tscope20->GetPointingResGBL(plane,Biased);

    tscope20->SetKappa(1.1);
    sig_kappa_plus = tscope20->GetPointingResGBL(plane,Biased);

    tscope20->SetKappa(1.0);

    std::cout << " sys. uncerts in E: sig_0 = " << sig_0*1e3 << " + (" << (sig_E_plus - sig_0)*1e3 << ") - (" << fabs(sig_E_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_E_plus - sig_0) + fabs(sig_E_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;

    std::cout << " sys. uncerts in KAPPA: sig_0 = " << sig_0*1e3 << " + (" << (sig_kappa_plus - sig_0)*1e3 << ") - (" << (sig_kappa_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_kappa_plus - sig_0) + fabs(sig_kappa_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;
  }


  if (runmode == 14)
  {
    std::cout << " This calculates the systematic uncertainties for the estiamted residual width for certain assumed input uncerts." << std::endl;

    double ebeam = 6.;
    // Create telescope
    //
    telescopebuild = "AnaTel_wide.geom";
    AnaTel *tscope150 = new AnaTel(telescopebuild.c_str(), ebeam);
    double sigma_int = 3.24e-3;
    tscope150->SetResolution(sigma_int);
    std::cout << " - Tscope created, dz = 150 mm, NO DUT" << std::endl;

    int plane = 0;
    int Biased = 1;


    double sig_0 = tscope150->GetWidthGBL(plane,Biased);

    tscope150->SetBeam(ebeam*0.95);
    double sig_E_minus = tscope150->GetWidthGBL(plane,Biased);

    tscope150->SetBeam(ebeam*1.05);
    double sig_E_plus = tscope150->GetWidthGBL(plane,Biased);

    tscope150->SetBeam(ebeam);

    tscope150->SetKappa(0.9);
    double sig_kappa_minus = tscope150->GetWidthGBL(plane,Biased);

    tscope150->SetKappa(1.1);
    double sig_kappa_plus = tscope150->GetWidthGBL(plane,Biased);

    tscope150->SetKappa(1.0);

    std::cout << " sys. uncerts in E: sig_0 = " << sig_0*1e3 << " + (" << (sig_E_plus - sig_0)*1e3 << ") - (" << (sig_E_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_E_plus - sig_0) + fabs(sig_E_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;

    std::cout << " sys. uncerts in KAPPA: sig_0 = " << sig_0*1e3 << " + (" << (sig_kappa_plus - sig_0)*1e3 << ") - (" << fabs(sig_kappa_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_kappa_plus - sig_0) + fabs(sig_kappa_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;

    telescopebuild = "AnaTel_narrow.geom";
    AnaTel *tscope20 = new AnaTel(telescopebuild.c_str(), ebeam);
    tscope20->SetResolution(sigma_int);
    std::cout << "\n - Tscope created, dz = 20 mm, NO DUT" << std::endl;

    sig_0 = tscope20->GetWidthGBL(plane,Biased);

    tscope20->SetBeam(ebeam*0.95);
    sig_E_minus = tscope20->GetWidthGBL(plane,Biased);

    tscope20->SetBeam(ebeam*1.05);
    sig_E_plus = tscope20->GetWidthGBL(plane,Biased);

    tscope20->SetBeam(ebeam);

    tscope20->SetKappa(0.9);
    sig_kappa_minus = tscope20->GetWidthGBL(plane,Biased);

    tscope20->SetKappa(1.1);
    sig_kappa_plus = tscope20->GetWidthGBL(plane,Biased);

    tscope20->SetKappa(1.0);

    std::cout << " sys. uncerts in E: sig_0 = " << sig_0*1e3 << " + (" << (sig_E_plus - sig_0)*1e3 << ") - (" << fabs(sig_E_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_E_plus - sig_0) + fabs(sig_E_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;

    std::cout << " sys. uncerts in KAPPA: sig_0 = " << sig_0*1e3 << " + (" << (sig_kappa_plus - sig_0)*1e3 << ") - (" << (sig_kappa_minus - sig_0)*1e3 << ")" << std::endl;
    std::cout << " rel sys. uncerts: " << (fabs(sig_kappa_plus - sig_0) + fabs(sig_kappa_minus - sig_0))/2. / sig_0*100 << " %" << std::endl;
  }

  if (runmode == 13)
  {
    std::cout << " This plots the estimated residual width as a function of the assumed intrinsic resolution." << std::endl;

    double ebeam = 3.;
    // Create telescope
    //
    telescopebuild = "AnaTel_wide.geom";
    AnaTel *tscope150 = new AnaTel(telescopebuild.c_str(), ebeam);
    double sigma_int = 3.24e-3;
    tscope150->SetResolution(sigma_int);
    std::cout << " - Tscope created, dz = 150 mm, NO DUT" << std::endl;

    std::cout << " --- BIASED RESIDUALS ---" << std::endl;
    std::cout << "     Pointing reso estimate at plane 0: " << tscope150->GetPointingResGBL(0,1) 
      << "   Residual estimate = "              << tscope150->GetWidthGBL(0,1) << std::endl;

    // keep kappa fixed, vary sigma int
    const int steps_sigma = 41;
    double step_width_sigma = 0.1e-3;
    double start_sigma = 1.0e-3;
    std::vector<double> sigmas;

    for(unsigned int ii = 0; ii<steps_sigma; ii++) sigmas.push_back(start_sigma + ii*step_width_sigma);

    //srd::vector<double> rhat_sigma(31);


    TGraph* gr_rhat150 = new TGraph(31);
    gr_rhat150->SetNameTitle("","");
    double rhat = -1.;
    for (int ii = 0; ii < steps_sigma; ii++){
      std::cout << " int res is now " << sigmas.at(ii) << std::endl;
      tscope150->SetResolution(sigmas.at(ii));
      rhat = tscope150->GetRhat(2,1);
      std::cout << " rhat is " << rhat << std::endl;
      gr_rhat150->SetPoint(ii, sigmas.at(ii)*1e3, rhat*1e3);
      rhat = -1.;
    }

    gr_rhat150->GetXaxis()->SetLimits(0., sigmas.back()*1.1*1e3);
    gr_rhat150->SetMinimum(0.);
    gr_rhat150->GetXaxis()->SetTitle("intrinsic reso. #sigma_{int} [#mum]");
    gr_rhat150->GetYaxis()->SetTitle("intrinsic reso. error plane 3 [#mum]"); // i.e. estimated width of biased residual distribution

    gr_rhat150->SetMarkerStyle(20);

    TCanvas *cmg_rhat_150 = new TCanvas("cmg_rhat_150","cmg1_rhat_50",10,10,800,600);
    cmg_rhat_150->cd();

    gr_rhat150->Draw("apl");

    cmg_rhat_150->Write();


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
    for (Int_t j=0; j<threshcount; j++)
    {
      std::cout << " j = " << j << std::endl;


      thres.push_back(j+(13-threshcount));
      ethres.push_back(.0);

      CSchecker( run_2_150[j], 2.0, j ); 
      CSchecker( run_3_150[j], 3.0, j );
      CSchecker( run_4_150[j], 4.0, j );
      CSchecker( run_5_150[j], 5.0, j );
      CSchecker( run_6_150[j], 6.0, j );
      //CSchecker( run120[j],  120.0, j );

      CSchecker( run_2_20[j], 2.0, j );
      CSchecker( run_3_20[j], 3.0, j );
      CSchecker( run_4_20[j], 4.0, j );
      CSchecker( run_5_20[j], 5.0, j );
      CSchecker( run_6_20[j], 6.0, j );

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

    pulls20.open("pulls20-xx.txt");
    pulls150.open("pulls150-xx.txt");
    //pulls20.open("pulls20-12.txt");
    //pulls150.open("pulls150-12.txt");

    telescopebuild = "AnaTel_narrow.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    //fitter(118,6.);
    for( int i : run_6_20) fitter(i,6.0);
    for( int i : run_5_20) fitter(i,5.0);
    for( int i : run_4_20) fitter(i,4.0);
    for( int i : run_3_20) fitter(i,3.0);
    for( int i : run_2_20) fitter(i,2.0);


    telescopebuild = "AnaTel_wide.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    //fitter(64,6.);
    for( int i : run_6_150) fitter(i,6.0);
    for( int i : run_5_150) fitter(i,5.0);
    for( int i : run_4_150) fitter(i,4.0);
    for( int i : run_3_150) fitter(i,3.0);
    for( int i : run_2_150) fitter(i,2.0);
    //for( int i : run120)    fitter(i,120.0);


    //std::cout << "fill tgraphs" << std::endl;
    FillGraph();

    pulls20.close();
    pulls150.close();

  }

  if (runmode == 10)
  {


    double ebeam = 3.;
    // Create telescope
    //
    telescopebuild = "AnaTel_wide.geom";
    AnaTel *tscope150 = new AnaTel(telescopebuild.c_str(), ebeam);
    double sigma_int = 3.3e-3;
    tscope150->SetResolution(sigma_int);
    std::cout << " - Tscope created, dz = 150 mm, NO DUT" << std::endl;

    // Give original beam energy with spread // assume spread to be negligible
    std::cout << " E = " << ebeam << " , sigma_int " << sigma_int<< std::endl;

    // unbiased residuals
    std::cout << " --- UNBIASED RESIDUALS ---" << std::endl;
    std::cout << "     Pointing reso estimate at plane 0: " << tscope150->GetPointingRes(0,0) 
      << "   Residual estimate = "              << tscope150->GetWidth(0,0) << std::endl;

    std::cout << "     Pointing reso estimate at plane 3: " << tscope150->GetPointingRes(3,0) 
      << "   Residual estimate = "              << tscope150->GetWidth(3,0) << std::endl;

    for(int ii = 0; ii < 6; ii++){
      std::cout << "\n " << ii << std::endl;
      std::cout << "     GBL Pointing reso estimate at plane 0: " << tscope150->GetPointingResGBL(ii,0) 
	<< "   Residual estimate = " 		       << tscope150->GetWidthGBL(ii,0) << std::endl;
    }





    // now biased redidual estiamtes
    std::cout << "\n --- BIASED RESIDUALS ---" << std::endl;
    std::cout << "     Pointing reso estimate at plane 0: " << tscope150->GetPointingRes(0,1) 
      << "   Residual estimate = "              << tscope150->GetWidth(0,1) << std::endl;

    std::cout << "     Pointing reso estimate at plane 3: " << tscope150->GetPointingRes(3,1) 
      << "   Residual estimate = "              << tscope150->GetWidth(3,1) << std::endl;

    for(int ii = 0; ii < 6; ii++){
      std::cout << "\n " << ii << std::endl;
      std::cout << "     GBL Pointing reso estimate at plane 0: " << tscope150->GetPointingResGBL(ii,1) 
	<< "   Residual estimate = " 		       << tscope150->GetWidthGBL(ii,1)
	<< "   sigma_int estimate = " 		       << sqrt(tscope150->GetWidthGBL(ii,1)*tscope150->GetWidthGBL(ii,0)) << std::endl;
      std::cout << "     GBL kink Res Error estimate at plane 0: " << tscope150->GetkResError(ii,1)*1e3 << std::endl;
    }






    // Create telescope
    //
    telescopebuild = "AnaTel_narrow.geom";
    AnaTel *tscope20 = new AnaTel(telescopebuild.c_str(), ebeam);
    sigma_int = 3.3e-3;
    tscope20->SetResolution(sigma_int);
    std::cout << "\nTscope created, dz = 20 mm, NO DUT" << std::endl;

    // Give original beam energy with spread // assume spread to be negligible
    std::cout << " E = " << ebeam << " , sigma_int " << sigma_int<< std::endl;

    std::cout << "     Pointing reso estimate at plane 0: " << tscope20->GetPointingRes(0,1) 
      << "   Residual estimate = "              << tscope20->GetWidth(0,1) << std::endl;

    std::cout << "     Pointing reso estimate at plane 3: " << tscope20->GetPointingRes(3,1) 
      << "   Residual estimate = "              << tscope20->GetWidth(3,1) << std::endl;


    for(int ii = 0; ii < 6; ii++){
      std::cout << "\n " << ii << std::endl;
      std::cout << "     GBL Pointing reso estimate at plane 0: " << tscope20->GetPointingResGBL(ii,1) 
	<< "   Residual estimate = " 		       << tscope20->GetWidthGBL(ii,1)
	<< "   sigma_int estimate = " 		       << sqrt(tscope20->GetWidthGBL(ii,1)*tscope20->GetWidthGBL(ii,0)) << std::endl;
      std::cout << "     GBL kink Res Error estimate at plane 0: " << tscope20->GetkResError(ii,1)*1e3 << std::endl;
    }



    // Create telescope
    //
    telescopebuild = "AnaTel_middle.geom";
    AnaTel *tscope50 = new AnaTel(telescopebuild.c_str(), ebeam);
    sigma_int = 3.3e-3;
    tscope50->SetResolution(sigma_int);
    std::cout << "\nTscope created, dz = 50 mm, NO DUT" << std::endl;

    // Give original beam energy with spread // assume spread to be negligible
    std::cout << " E = " << ebeam << " , sigma_int " << sigma_int<< std::endl;

    std::cout << "     Pointing reso estimate at plane 0: " << tscope50->GetPointingRes(0,1) 
      << "   Residual estimate = "              << tscope50->GetWidth(0,1) << std::endl;

    std::cout << "     Pointing reso estimate at plane 3: " << tscope50->GetPointingRes(3,1) 
      << "   Residual estimate = "              << tscope50->GetWidth(3,1) << std::endl;

    for(int ii = 0; ii < 6; ii++){
      std::cout << "\n " << ii << std::endl;
      std::cout << "     GBL Pointing reso estimate at plane 0: " << tscope50->GetPointingResGBL(ii,1) 
	<< "   Residual estimate = " 		       << tscope50->GetWidthGBL(ii,1)
	<< "   sigma_int estimate = " 		       << sqrt(tscope50->GetWidthGBL(ii,1)*tscope50->GetWidthGBL(ii,0)) << std::endl;
      std::cout << "     GBL kink Res Error estimate at plane 0: " << tscope50->GetkResError(ii,1)*1e3 << std::endl;
    }

  }

  if (runmode == 9)
  {
    global_thresh=6;



    telescopebuild = "AnaTel_wide.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    fitter(63,6.);
    //fitter(73,5.);
    //fitter(84,4.);
    //fitter(96,3.);
    //fitter(106,2.);

    // CERN
    //fitter(755,120.);



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

    telescopebuild = "AnaTel_wide.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = 150.0*j;

    for(int i=0;i<run_6_150.size();i++)
    {
      global_thresh = i+3;
      fitter( run_6_150[i], 6.0 );
    }
    for(int i=0;i<run_5_150.size();i++)
    {
      global_thresh = i+3;
      fitter( run_5_150[i], 5.0 );
    }
    for(int i=0;i<run_4_150.size();i++)
    {
      global_thresh = i+3;
      fitter( run_4_150[i], 4.0 );
    }
    for(int i=0;i<run_3_150.size();i++)
    {
      global_thresh = i+3;
      fitter( run_3_150[i], 3.0 );
    }
    for(int i=0;i<run_2_150.size();i++)
    {
      global_thresh = i+3;
      fitter( run_2_150[i], 2.0 );
    }

    std::cout << " " << std::endl;
    std::cout << "Thin Geometry" << std::endl;
    std::cout << " " << std::endl;

    telescopebuild = "AnaTel_narrow.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = 20.0*j;

    for(int i=0;i<run_6_20.size();i++)
    {
      global_thresh = i+3;
      fitter( run_6_20[i], 6.0 );
    }
    for(int i=0;i<run_5_20.size();i++)
    {
      global_thresh = i+3;
      fitter( run_5_20[i], 5.0 );
    }
    for(int i=0;i<run_4_20.size();i++)
    {
      global_thresh = i+3;
      fitter( run_4_20[i], 4.0 );
    }
    for(int i=0;i<run_3_20.size();i++)
    {
      global_thresh = i+3;
      fitter( run_3_20[i], 3.0 );
    }
    for(int i=0;i<run_2_20.size();i++)
    {
      global_thresh = i+3;
      fitter( run_2_20[i], 2.0 );
    }

    telescopebuild = "AnaTel_wide.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = 150.0*j;

    for(int i=0;i<run120.size();i++)
    {
      global_thresh = i+3;
      fitter( run120[i], 120.0 );
    }

    std::cout << "fill tgraphs" << std::endl;
    FillGraph();

  }

  // Run according to clustersize, this should be at one threshold -> to be done -> Geometry needs implementing
  if(runmode == 1)
  {
    std::cout << "Running over clustersize" << std::endl;
    // Smallest clustersize is 1, largest is clustercount -1
    const Int_t clustercount = 5;
    Double_t clusterresult[clustercount];
    Double_t clustererror[clustercount];
    Double_t x[clustercount];
    for (Int_t j=1;j<clustercount;j++)
    {
      x[j] = j+0.0;
      std::ostringstream convert;
      convert << j;
      submask = convert.str();

      // this needs to be done:
      for(int i=0;i<run_3_150.size();i++)
	fitter( run_3_150[i], 3 );
      clusterresult[j] = m26_resolution*1000.0;
      clustererror[j] = m26_res_error*1000.0;
    }

    // As a comparison: put "clustersize = 0" as the normal uncut run.
    submask = "";
    for(int i=0;i<run_3_150.size();i++)
      fitter( run_3_150[i], 3 );
    clusterresult[0] = m26_resolution*1000.0;
    clustererror[0] = m26_res_error*1000.0;

    TCanvas *c1 = new TCanvas("c1","Resolution vs. Clustersize",10,10,800,600);
    c1->SetFillColor(0);
    c1->SetGrid();
    TGraph *gr = new TGraph(clustercount,x,clusterresult);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(21);
    gr->SetTitle("Resolution vs. Clustersize");
    gr->GetXaxis()->SetTitle("Clustersize");
    gr->GetYaxis()->SetTitle("#sigma_{M26} in #mum");
    gr->Draw("ACP");
    c1->Update();
    //   c1->GetFrame()->SetFillColor(0);
    //    c1->GetFrame()->SetBorderSize(0);
    c1->Modified();
    c1->Print("../bin/pics/clustersize_somename.eps");
    _outputFile->cd();
    c1->Write();
    c1->Close();
  }

  // Run over thresholds
  if(runmode == 2)
  {
    std::cout << "Threshold mode" << std::endl;
    Double_t threshresult_2_150[threshcount];
    Double_t threshresult_3_150[threshcount];
    Double_t threshresult_4_150[threshcount];
    Double_t threshresult_5_150[threshcount];
    Double_t threshresult_6_150[threshcount];
    Double_t threshresult_2_20[threshcount];
    Double_t threshresult_3_20[threshcount];
    Double_t threshresult_4_20[threshcount];
    Double_t threshresult_5_20[threshcount];
    Double_t threshresult_6_20[threshcount];
    Double_t threshresult120[threshcount] = {0.};


    Double_t thresherror_2_150[threshcount];
    Double_t thresherror_3_150[threshcount];
    Double_t thresherror_4_150[threshcount];
    Double_t thresherror_5_150[threshcount];
    Double_t thresherror_6_150[threshcount];
    Double_t thresherror_2_20[threshcount];
    Double_t thresherror_3_20[threshcount];
    Double_t thresherror_4_20[threshcount];
    Double_t thresherror_5_20[threshcount];
    Double_t thresherror_6_20[threshcount];
    Double_t thresherror120[threshcount];


    Double_t x[threshcount];
    Double_t xerrorthresh[threshcount];
    for (Int_t j=0; j<threshcount; j++)
    {
      std::cout << " j = " << j << std::endl;

      telescopebuild = "AnaTel_wide.geom";
      planedistance = 150;
      for(int jj=0;jj<nplanes;jj++)
	posx[jj] = planedistance*jj;


      x[j] = j+(3+(10-threshcount));
      global_thresh = x[j];
      std::cout << " threhold = " << global_thresh << std::endl;
      xerrorthresh[j] = 0.0;
      submask = "";


      fitter( run_2_150[j], 2.0 ); 
      threshresult_2_150[j] = m26_resolution*1000.0;
      thresherror_2_150[j] = global_plot_error; // FIXME need correct error estimate here

      fitter( run_3_150[j], 3.0 );
      threshresult_3_150[j] = m26_resolution*1000.0;
      thresherror_3_150[j] = global_plot_error;

      fitter( run_4_150[j], 4.0 );
      threshresult_4_150[j] = m26_resolution*1000.0;
      thresherror_4_150[j] = global_plot_error;

      fitter( run_5_150[j], 5.0 );
      threshresult_5_150[j] = m26_resolution*1000.0;
      thresherror_5_150[j] = global_plot_error;

      fitter( run_6_150[j], 6.0 );
      threshresult_6_150[j] = m26_resolution*1000.0;
      thresherror_6_150[j] = global_plot_error;


      //if( run120[j] == 752 || run120[j] == 753)
      //{
      //  threshresult120[j] = 1.0;
      //  thresherror120[j] = 0.0;
      //} else {
      fitter( run120[j], 120.0 );
      threshresult120[j] = m26_resolution*1000.0;
      thresherror120[j] = global_plot_error;
      //}

      telescopebuild = "AnaTel_narrow.geom";
      planedistance = 20;
      for(int jj=0;jj<nplanes;jj++)
	posx[jj] = planedistance*jj;



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

      fitter( run_6_20[j], 6.0 );
      threshresult_6_20[j] = m26_resolution*1000.0;
      thresherror_6_20[j] = global_plot_error;


    }

    std::cout << "fill tgraphs" << std::endl;
    FillGraph();

    TGraphErrors *gr2_150 = new TGraphErrors((threshcount),x,threshresult_2_150,xerrorthresh,thresherror_2_150);
    TGraphErrors *gr3_150 = new TGraphErrors((threshcount),x,threshresult_3_150,xerrorthresh,thresherror_3_150);
    TGraphErrors *gr4_150 = new TGraphErrors((threshcount),x,threshresult_4_150,xerrorthresh,thresherror_4_150);
    TGraphErrors *gr5_150 = new TGraphErrors((threshcount),x,threshresult_5_150,xerrorthresh,thresherror_5_150);
    TGraphErrors *gr6_150 = new TGraphErrors((threshcount),x,threshresult_6_150,xerrorthresh,thresherror_6_150);
    TGraphErrors *gr120   = new TGraphErrors((threshcount),x,threshresult120,xerrorthresh,thresherror120);

    TGraphErrors *gr2_20  = new TGraphErrors((threshcount),x,threshresult_2_20,xerrorthresh,thresherror_2_20);
    TGraphErrors *gr3_20  = new TGraphErrors((threshcount),x,threshresult_3_20,xerrorthresh,thresherror_3_20);
    TGraphErrors *gr4_20  = new TGraphErrors((threshcount),x,threshresult_4_20,xerrorthresh,thresherror_4_20);
    TGraphErrors *gr5_20  = new TGraphErrors((threshcount),x,threshresult_5_20,xerrorthresh,thresherror_5_20);
    TGraphErrors *gr6_20  = new TGraphErrors((threshcount),x,threshresult_6_20,xerrorthresh,thresherror_6_20);

    TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
    threshold->cd();


    //TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
    //gStyle->SetPadBorderMode(0);
    //gStyle->SetOptStat(0);
    threshold->SetFillColor(0);
    //gPad->SetGridx();
    //gPad->SetGridy();

    gr2_150->SetMarkerStyle(20);
    gr2_150->SetMarkerColor(kBlack);
    gr2_150->SetMarkerSize(3);
    gr3_150->SetMarkerStyle(21);
    gr3_150->SetMarkerColor(kGreen);
    gr3_150->SetMarkerSize(3);
    gr4_150->SetMarkerStyle(22);
    gr4_150->SetMarkerColor(kRed);
    gr4_150->SetMarkerSize(3);
    gr5_150->SetMarkerStyle(23);
    gr5_150->SetMarkerColor(kBlue);
    gr5_150->SetMarkerSize(3);
    gr6_150->SetMarkerStyle(29);
    gr6_150->SetMarkerColor(kOrange);
    gr6_150->SetMarkerSize(3);
    gr120->SetMarkerStyle(3);
    gr120->SetMarkerColor(kMagenta);
    gr120->SetMarkerSize(3);

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
    gr6_20->SetMarkerStyle(28);
    gr6_20->SetMarkerColor(kOrange);
    gr6_20->SetMarkerSize(3);

    TMultiGraph* mg = new TMultiGraph();

    mg->Add(gr2_150,"p");
    mg->Add(gr3_150,"p");
    mg->Add(gr4_150,"p");
    mg->Add(gr5_150,"p");
    mg->Add(gr6_150,"p");
    mg->Add(gr120,"p");
    mg->Add(gr2_20,"p");
    mg->Add(gr3_20,"p");
    mg->Add(gr4_20,"p");
    mg->Add(gr5_20,"p");
    mg->Add(gr6_20,"p");

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
    leg->AddEntry(gr2_150,"p = 2 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr3_150,"p = 3 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr4_150,"p = 4 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr5_150,"p = 5 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr6_150,"p = 6 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr120,"p = 120 GeV, #Delta_{z} = 150 mm","p");

    leg->AddEntry(gr2_20,"p = 2 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr3_20,"p = 3 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr4_20,"p = 4 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr5_20,"p = 5 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr6_20,"p = 6 GeV, #Delta_{z} = 20 mm","p");

    leg->Draw();

    // Output
    threshold->Print("pics/threshold_all.eps");
    _outputFile->cd();
    threshold->Write();
    threshold->Close();
  }

  // Run over thresholds (x) and clustersize (y), res is z, geometry needs reset and thin implementing
  if(runmode == 3)
  {
    std::cout << "Threshold and clustersize mode" << std::endl;
    //const Int_t threshcount = 12;
    const Int_t clustercount = 5;
    Double_t threshresult2[threshcount][clustercount];
    Double_t threshresult3[threshcount][clustercount];
    Double_t threshresult4[threshcount][clustercount];
    Double_t threshresult5[threshcount][clustercount];
    Double_t threshresulttotal[threshcount][clustercount];
    Double_t thresherror2[threshcount][clustercount];
    Double_t thresherror3[threshcount][clustercount];
    Double_t thresherror4[threshcount][clustercount];
    Double_t thresherror5[threshcount][clustercount];
    Double_t thresherrortotal[threshcount][clustercount];
    Double_t x[threshcount][clustercount];
    Double_t xerrorthresh[threshcount][clustercount];
    Double_t y[threshcount][clustercount];
    Double_t yerrorcluster[threshcount][clustercount];

    for (Int_t i=0;i<clustercount;i++)
    {
      for (Int_t j=0;j<(threshcount-2);j++)
      {
	std::cout << "I am now at clustersize " << i << " and threshold " << j+3 << " !" << std::endl;
	y[j][i] = i;
	x[j][i] = j+3;
	global_thresh = x[j][i];
	xerrorthresh[j][i] = 0.0;
	std::ostringstream convert;
	convert << i;
	submask = convert.str();
	if (i == 0)
	  submask = "";
	/*	fitter( run_2_150[j], 2.0 );
		threshresult2[j][i] = m26_resolution*1000.0;
		thresherror2[j][i] = m26_res_error*1000.0;
		fitter( run_3_150[j], 3.0 );
		threshresult3[j][i] = m26_resolution*1000.0;
		thresherror3[j][i] = m26_res_error*1000.0; */
	fitter( run_4_150[j], 4.4 );
	threshresult4[j][i] = m26_resolution*1000.0;
	thresherror4[j][i] = m26_res_error*1000.0;
	/*	fitter( run_5_150[j], 5.0 );
		threshresult5[j][i] = m26_resolution*1000.0;
		thresherror5[j][i] = m26_res_error*1000.0; */
      }
    }

    /*   TCanvas *c_thr_vs_clu_2 = new TCanvas("thr_vs_clu_2", "Threshold vs. Clustersize", 600, 400);
	 c_thr_vs_clu_2->cd();
	 TH2D *hist2D_thr_vs_clu_2 = new TH2D("hist2D_thr_vs_clu_2", "Histo_thr_vs_clu_2", 10, 3., 13., 5, 0., 5.);
	 for(Int_t i=0;i<clustercount;i++)
	 {
	 for(Int_t j=0;j<(threshcount-2);j++)
	 {
	 hist2D_thr_vs_clu_2->Fill((j+3),i,threshresult2[j][i]);
	 }
	 }
	 hist2D_thr_vs_clu_2->GetXaxis()->SetTitle("Threshold");
	 hist2D_thr_vs_clu_2->GetYaxis()->SetTitle("Clustersize");
	 hist2D_thr_vs_clu_2->GetZaxis()->SetTitle("#sigma_{M26}");
	 hist2D_thr_vs_clu_2->Draw("LEGO2");
	 c_thr_vs_clu_2->Print("pics/thr_vs_clu_2.eps");
	 c_thr_vs_clu_2->Write();

	 TCanvas *c_thr_vs_clu_3 = new TCanvas("thr_vs_clu_3", "Threshold vs. Clustersize", 600, 400);
	 c_thr_vs_clu_3->cd();
	 TH2D *hist2D_thr_vs_clu_3 = new TH2D("hist2D_thr_vs_clu_3", "Histo_thr_vs_clu_3", 10, 3., 13., 5, 0., 5.);
	 for(Int_t i=0;i<clustercount;i++)
	 {
	 for(Int_t j=0;j<(threshcount-2);j++)
	 {
	 hist2D_thr_vs_clu_3->Fill((j+3),i,threshresult3[j][i]);
	 }
	 }
	 hist2D_thr_vs_clu_3->GetXaxis()->SetTitle("Threshold");
	 hist2D_thr_vs_clu_3->GetYaxis()->SetTitle("Clustersize");
	 hist2D_thr_vs_clu_3->GetZaxis()->SetTitle("#sigma_{M26}");
	 hist2D_thr_vs_clu_3->Draw("LEGO2");
	 c_thr_vs_clu_3->Print("pics/thr_vs_clu_3.eps");
	 c_thr_vs_clu_3->Write();
	 */
    TCanvas *c_thr_vs_clu_4 = new TCanvas("thr_vs_clu_4", "Threshold vs. Clustersize", 600, 400);
    c_thr_vs_clu_4->cd();
    TH2D *hist2D_thr_vs_clu_4 = new TH2D("hist2D_thr_vs_clu_4", "Histo_thr_vs_clu_4", 10, 3., 13., 5, 0., 5.);
    for(Int_t i=0;i<clustercount;i++)
    {
      for(Int_t j=0;j<(threshcount-2);j++)
      {
	hist2D_thr_vs_clu_4->Fill((j+3),i,threshresult4[j][i]);
      }
    }
    hist2D_thr_vs_clu_4->GetXaxis()->SetTitle("Threshold");
    hist2D_thr_vs_clu_4->GetYaxis()->SetTitle("Clustersize");
    hist2D_thr_vs_clu_4->GetZaxis()->SetTitle("#sigma_{M26}");
    hist2D_thr_vs_clu_4->Draw("LEGO2");
    c_thr_vs_clu_4->Print("pics/thr_vs_clu_4.eps");
    c_thr_vs_clu_4->Write();

    /*  TCanvas *c_thr_vs_clu_5 = new TCanvas("thr_vs_clu_5", "Threshold vs. Clustersize", 600, 400);
	c_thr_vs_clu_5->cd();
	TH2D *hist2D_thr_vs_clu_5 = new TH2D("hist2D_thr_vs_clu_5", "Histo_thr_vs_clu_5", 10, 3., 13., 5, 0., 5.);
	for(Int_t i=0;i<clustercount;i++)
	{
	for(Int_t j=0;j<(threshcount-2);j++)
	{
	hist2D_thr_vs_clu_5->Fill((j+3),i,threshresult5[j][i]);
	}
	}
	hist2D_thr_vs_clu_5->GetXaxis()->SetTitle("Threshold");
	hist2D_thr_vs_clu_5->GetYaxis()->SetTitle("Clustersize");
	hist2D_thr_vs_clu_5->GetZaxis()->SetTitle("#sigma_{M26}");
	hist2D_thr_vs_clu_5->Draw("LEGO2");
	c_thr_vs_clu_5->Print("pics/thr_vs_clu_5.eps");
	c_thr_vs_clu_5->Write();



	TCanvas *c_thr_vs_clu_total = new TCanvas("thr_vs_clu_total", "Threshold vs. Clustersize", 600, 400);
	c_thr_vs_clu_total->cd();
	TH2D *hist2D_thr_vs_clu_total = new TH2D("hist2D_thr_vs_clu_total", "Histo_thr_vs_clu_total", 10, 3., 13., 5, 0., 5.);
	for(Int_t i=0;i<clustercount;i++)
	{
	for(Int_t j=0;j<(threshcount-2);j++)
	{
	threshresulttotal[j][i]=((threshresult2[j][i]+threshresult3[j][i]+threshresult4[j][i]+threshresult5[j][i])/4.);
	hist2D_thr_vs_clu_total->Fill((j+3),i,threshresulttotal[j][i]);
	}
	}
	hist2D_thr_vs_clu_total->GetXaxis()->SetTitle("Threshold");
	hist2D_thr_vs_clu_total->GetYaxis()->SetTitle("Clustersize");
	hist2D_thr_vs_clu_total->GetZaxis()->SetTitle("#sigma_{M26}");
	hist2D_thr_vs_clu_total->Draw("LEGO2");
	c_thr_vs_clu_total->Print("pics/thr_vs_clu_total.eps");
	c_thr_vs_clu_total->Write();

*/
    /*

       TGraphErrors *gr2 = new TGraphErrors((threshcount-2),x,threshresult2,xerrorthresh,thresherror2);
       TGraphErrors *gr3 = new TGraphErrors((threshcount-2),x,threshresult3,xerrorthresh,thresherror3);
       TGraphErrors *gr4 = new TGraphErrors((threshcount-2),x,threshresult4,xerrorthresh,thresherror4);
       TGraphErrors *gr5 = new TGraphErrors((threshcount-2),x,threshresult5,xerrorthresh,thresherror5);
       TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
       TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
       gStyle->SetPadBorderMode(0);
       gStyle->SetOptStat(0);
       threshold->SetFillColor(0);
       threshold->Divide(1,1);
       threshold->cd(1);
       gStyle->SetErrorX(0);
       gPad->SetLogx(0);
       gPad->SetLogy(0);
       gr2->SetMarkerStyle(20);
       gr2->SetMarkerColor(kBlack);
       gr2->SetMarkerSize(2);
       gr3->SetMarkerStyle(21);
       gr3->SetMarkerColor(kGreen);
       gr3->SetMarkerSize(2);
       gr4->SetMarkerStyle(22);
       gr4->SetMarkerColor(kRed);
       gr4->SetMarkerSize(2);
       gr5->SetMarkerStyle(23);
       gr5->SetMarkerColor(kBlue);
       gr5->SetMarkerSize(2);
       histo_cfg(h_axis, "Threshold (s/n)","#sigma_{M26} (#mum)","");
       th_axis->SetMinimum(0.0);
       th_axis->SetMaximum(5.0);
       th_axis->Draw("hist");
       gr2->Draw("P");
       gr3->Draw("P");
       gr4->Draw("P");
       gr5->Draw("P");

       TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
       leg->SetBorderSize(0);
       leg->SetFillColor(0);
       leg->SetFillStyle(0);
       leg->SetHeader("Energy:");
       leg->AddEntry(gr2,"p = 2 GeV","p");
       leg->AddEntry(gr3,"p = 3 GeV","p");
       leg->AddEntry(gr4,"p = 4.4 GeV","p");
       leg->AddEntry(gr5,"p = 5 GeV","p");
       leg->Draw();

    // Output
    threshold->Print("pics/threshold_all.eps"); */
  }


  // Effi is now done in runmode 17
  /*
  // Run over threshold to get efficiency and noise -> 120gev missing
  if(runmode == 4)
  {

  // Results go in here
  Double_t threshnoise2[threshcount];
  Double_t thresheffi2[threshcount];
  Double_t threshnoise3[threshcount];
  Double_t thresheffi3[threshcount];
  Double_t threshnoise4[threshcount];
  Double_t thresheffi4[threshcount];
  Double_t threshnoise5[threshcount];
  Double_t thresheffi5[threshcount];
  Double_t threshnoise6[threshcount];
  Double_t thresheffi6[threshcount];
  Double_t threshnoise2_error[threshcount];
  Double_t thresheffi2_error[threshcount];
  Double_t threshnoise3_error[threshcount];
  Double_t thresheffi3_error[threshcount];
  Double_t threshnoise4_error[threshcount];
  Double_t thresheffi4_error[threshcount];
  Double_t threshnoise5_error[threshcount];
  Double_t thresheffi5_error[threshcount];
  Double_t threshnoise6_error[threshcount];
  Double_t thresheffi6_error[threshcount];
  Double_t thin_threshnoise2[threshcount];
  Double_t thin_thresheffi2[threshcount];
  Double_t thin_threshnoise3[threshcount];
  Double_t thin_thresheffi3[threshcount];
  Double_t thin_threshnoise4[threshcount];
  Double_t thin_thresheffi4[threshcount];
  Double_t thin_threshnoise5[threshcount];
  Double_t thin_thresheffi5[threshcount];
  Double_t thin_threshnoise6[threshcount];
  Double_t thin_thresheffi6[threshcount];
  Double_t thin_threshnoise2_error[threshcount];
  Double_t thin_thresheffi2_error[threshcount];
  Double_t thin_threshnoise3_error[threshcount];
  Double_t thin_thresheffi3_error[threshcount];
  Double_t thin_threshnoise4_error[threshcount];
  Double_t thin_thresheffi4_error[threshcount];
  Double_t thin_threshnoise5_error[threshcount];
  Double_t thin_thresheffi5_error[threshcount];
  Double_t thin_threshnoise6_error[threshcount];
  Double_t thin_thresheffi6_error[threshcount];
  Double_t x[threshcount];
  Double_t xerror[threshcount] = {0.0};

  std::cout << " " << std::endl;
  std::cout << "Mode 4" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Running over all runs - efficiency and noise" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Wide Geometry" << std::endl;
  std::cout << " " << std::endl;

  telescopebuild = "AnaTel_wide.geom";
  planedistance = 150;
  for(int j=0;j<nplanes;j++)
  posx[j] = 150.0*j;

  for(int i=0;i<run_2_150.size();i++)
  {
  noise( run_2_150[i] );
  x[i] = i+3-0.1;
  threshnoise2[i] = avgnoise;
  thresheffi2[i] = avgeffi;
  threshnoise2_error[i] = avgnoise_error;
  thresheffi2_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_3_150.size();i++)
  {
    noise( run_3_150[i] );
    x[i] = i+3-0.05;
    threshnoise3[i] = avgnoise;
    thresheffi3[i] = avgeffi;
    threshnoise3_error[i] = avgnoise_error;
    thresheffi3_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_4_150.size();i++)
  {
    noise( run_4_150[i] );
    x[i] = i+3;
    threshnoise4[i] = avgnoise;
    thresheffi4[i] = avgeffi;
    threshnoise4_error[i] = avgnoise_error;
    thresheffi4_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_5_150.size();i++)
  {
    noise( run_5_150[i] );
    x[i] = i+3+0.05;
    threshnoise5[i] = avgnoise;
    thresheffi5[i] = avgeffi;
    threshnoise5_error[i] = avgnoise_error;
    thresheffi5_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_6_150.size();i++)
  {
    noise( run_6_150[i] );
    x[i] = i+3+0.01;
    threshnoise6[i] = avgnoise;
    thresheffi6[i] = avgeffi;
    threshnoise6_error[i] = avgnoise_error;
    thresheffi6_error[i] = avgeffi_error;
  }


  std::cout << " " << std::endl;
  std::cout << "Thin Geometry" << std::endl;
  std::cout << " " << std::endl;

  telescopebuild = "AnaTel_narrow.geom";
  planedistance = 20;

  for(int j=0;j<nplanes;j++)
    posx[j] = 20.0*j;

  for(int i=0;i<run_2_20.size();i++)
  {
    noise( run_2_20[i] );
    x[i] = i+3-0.1;
    thin_threshnoise2[i] = avgnoise;
    thin_thresheffi2[i] = avgeffi;
    thin_threshnoise2_error[i] = avgnoise_error;
    thin_thresheffi2_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_3_20.size();i++)
  {
    noise( run_3_20[i] );
    x[i] = i+3-0.05;
    thin_threshnoise3[i] = avgnoise;
    thin_thresheffi3[i] = avgeffi;
    thin_threshnoise3_error[i] = avgnoise_error;
    thin_thresheffi3_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_4_20.size();i++)
  {
    noise( run_4_20[i] );
    x[i] = i+3+0.;
    thin_threshnoise4[i] = avgnoise;
    thin_thresheffi4[i] = avgeffi;
    thin_threshnoise4_error[i] = avgnoise_error;
    thin_thresheffi4_error[i] = avgeffi_error;
  }


  for(int i=0;i<run_5_20.size();i++)
  {
    noise( run_5_20[i] );
    x[i] = i+3+0.05;
    thin_threshnoise5[i] = avgnoise;
    thin_thresheffi5[i] = avgeffi;
    thin_threshnoise5_error[i] = avgnoise_error;
    thin_thresheffi5_error[i] = avgeffi_error;
  }

  for(int i=0;i<run_6_20.size();i++)
  {
    noise( run_6_20[i] );
    x[i] = i+3+0.1;
    thin_threshnoise6[i] = avgnoise;
    thin_thresheffi6[i] = avgeffi;
    thin_threshnoise6_error[i] = avgnoise_error;
    thin_thresheffi6_error[i] = avgeffi_error;
  }

  std::cout << " - Loop over all runs done." << std::endl;


  // Create graphs with the information
  TGraphErrors *gr2n = new TGraphErrors(10,x,threshnoise2,xerror,threshnoise2_error);
  TGraphErrors *gr2e = new TGraphErrors(10,x,thresheffi2,xerror,thresheffi2_error);
  TGraphErrors *gr3n = new TGraphErrors(10,x,threshnoise3,xerror,threshnoise3_error);
  TGraphErrors *gr3e = new TGraphErrors(10,x,thresheffi3,xerror,thresheffi3_error);
  TGraphErrors *gr4n = new TGraphErrors(10,x,threshnoise4,xerror,threshnoise4_error);
  TGraphErrors *gr4e = new TGraphErrors(10,x,thresheffi4,xerror,thresheffi4_error);
  TGraphErrors *gr5n = new TGraphErrors(10,x,threshnoise5,xerror,threshnoise5_error);
  TGraphErrors *gr5e = new TGraphErrors(10,x,thresheffi5,xerror,thresheffi5_error);
  TGraphErrors *gr6n = new TGraphErrors(10,x,threshnoise6,xerror,threshnoise6_error);
  TGraphErrors *gr6e = new TGraphErrors(10,x,thresheffi6,xerror,thresheffi6_error);

  TGraphErrors *gr2nth = new TGraphErrors(10,x,thin_threshnoise2,xerror,thin_threshnoise2_error);
  TGraphErrors *gr2eth = new TGraphErrors(10,x,thin_thresheffi2,xerror,thin_thresheffi2_error);
  TGraphErrors *gr3nth = new TGraphErrors(10,x,thin_threshnoise3,xerror,thin_threshnoise3_error);
  TGraphErrors *gr3eth = new TGraphErrors(10,x,thin_thresheffi3,xerror,thin_thresheffi3_error);
  TGraphErrors *gr4nth = new TGraphErrors(10,x,thin_threshnoise4,xerror,thin_threshnoise4_error);
  TGraphErrors *gr4eth = new TGraphErrors(10,x,thin_thresheffi4,xerror,thin_thresheffi4_error);
  TGraphErrors *gr5nth = new TGraphErrors(10,x,thin_threshnoise5,xerror,thin_threshnoise5_error);
  TGraphErrors *gr5eth = new TGraphErrors(10,x,thin_thresheffi5,xerror,thin_thresheffi5_error);
  TGraphErrors *gr6nth = new TGraphErrors(10,x,thin_threshnoise6,xerror,thin_threshnoise6_error);
  TGraphErrors *gr6eth = new TGraphErrors(10,x,thin_thresheffi6,xerror,thin_thresheffi6_error);

  std::cout << " Start Plotting." << std::endl;

  TMultiGraph* mg = new TMultiGraph();

  // Let's plot this
  //TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
  TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);
  threshold->SetFillColor(0);
  //threshold->Divide(1,1);
  threshold->cd();
  gStyle->SetErrorX(0);
  gPad->SetLogx(0);
  gPad->SetLogy(0);

  // Set apperance
  gr2n->SetMarkerStyle(21);
  gr2n->SetMarkerColor(kBlack);
  gr2n->SetMarkerSize(2);
  gr2e->SetMarkerStyle(20);
  gr2e->SetMarkerColor(kBlack);
  gr2e->SetMarkerSize(3);
  gr3n->SetMarkerStyle(21);
  gr3n->SetMarkerColor(kRed);
  gr3n->SetMarkerSize(2);
  gr3e->SetMarkerStyle(20);
  gr3e->SetMarkerColor(kRed);
  gr3e->SetMarkerSize(3);
  gr4n->SetMarkerStyle(21);
  gr4n->SetMarkerColor(kGreen);
  gr4n->SetMarkerSize(2);
  gr4e->SetMarkerStyle(20);
  gr4e->SetMarkerColor(kGreen);
  gr4e->SetMarkerSize(3);
  gr5n->SetMarkerStyle(21);
  gr5n->SetMarkerColor(kBlue);
  gr5n->SetMarkerSize(2);
  gr5e->SetMarkerStyle(20);
  gr5e->SetMarkerColor(kBlue);
  gr5e->SetMarkerSize(3);
  gr6n->SetMarkerStyle(21);
  gr6n->SetMarkerColor(kOrange);
  gr6n->SetMarkerSize(2);
  gr6e->SetMarkerStyle(20);
  gr6e->SetMarkerColor(kOrange);
  gr6e->SetMarkerSize(3);


  gr2nth->SetMarkerStyle(25);
  gr2nth->SetMarkerColor(kBlack);
  gr2nth->SetMarkerSize(2);
  gr2eth->SetMarkerStyle(24);
  gr2eth->SetMarkerColor(kBlack);
  gr2eth->SetMarkerSize(3);
  gr3nth->SetMarkerStyle(25);
  gr3nth->SetMarkerColor(kRed);
  gr3nth->SetMarkerSize(2);
  gr3eth->SetMarkerStyle(24);
  gr3eth->SetMarkerColor(kRed);
  gr3eth->SetMarkerSize(3);
  gr4nth->SetMarkerStyle(25);
  gr4nth->SetMarkerColor(kGreen);
  gr4nth->SetMarkerSize(2);
  gr4eth->SetMarkerStyle(24);
  gr4eth->SetMarkerColor(kGreen);
  gr4eth->SetMarkerSize(3);
  gr5nth->SetMarkerStyle(25);
  gr5nth->SetMarkerColor(kBlue);
  gr5nth->SetMarkerSize(2);
  gr5eth->SetMarkerStyle(24);
  gr5eth->SetMarkerColor(kBlue);
  gr5eth->SetMarkerSize(3);
  gr6nth->SetMarkerStyle(25);
  gr6nth->SetMarkerColor(kOrange);
  gr6nth->SetMarkerSize(2);
  gr6eth->SetMarkerStyle(24);
  gr6eth->SetMarkerColor(kOrange);
  gr6eth->SetMarkerSize(3);

  //histo_cfg(h_axis, "Threshold (s/n)","N","");


  mg->Add(gr2e,"P");
  mg->Add(gr3e,"P");
  mg->Add(gr4e,"P");
  mg->Add(gr5e,"P");
  mg->Add(gr6e,"P");
  mg->Add(gr2eth,"P");
  mg->Add(gr3eth,"P");
  mg->Add(gr4eth,"P");
  mg->Add(gr5eth,"P");
  mg->Add(gr6eth,"P");

  mg->Add(gr2n,"P");
  mg->Add(gr3n,"P");
  mg->Add(gr4n,"P");
  mg->Add(gr5n,"P");
  mg->Add(gr6n,"P");
  mg->Add(gr2nth,"P");
  mg->Add(gr3nth,"P");
  mg->Add(gr4nth,"P");
  mg->Add(gr5nth,"P");
  mg->Add(gr6nth,"P");

  mg->Draw("ap");
  mg->SetMinimum(0.0);
  mg->SetMaximum(1.0);
  mg->GetXaxis()->SetTitle("M26 threshold #xi_{n}");
  mg->GetYaxis()->SetTitle("M26 efficiency");

  // The legend
  TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("Performance:");
  leg->SetNColumns(2);

  leg->AddEntry(gr2e,"2 GeV efficiency","p");
  leg->AddEntry(gr2n,"2 GeV noise","p");
  leg->AddEntry(gr3e,"3 GeV efficiency","p");
  leg->AddEntry(gr3n,"3 GeV noise","p");
  leg->AddEntry(gr4e,"4 GeV efficiency","p");
  leg->AddEntry(gr4n,"4 GeV noise","p");
  leg->AddEntry(gr5e,"5 GeV efficiency","p");
  leg->AddEntry(gr5n,"5 GeV noise","p");
  leg->AddEntry(gr6e,"6 GeV efficiency","p");
  leg->AddEntry(gr6n,"6 GeV noise","p");

  leg->AddEntry(gr2eth,"2 GeV efficiency 20mm","p");
  leg->AddEntry(gr2nth,"2 GeV noise 20mm","p");
  leg->AddEntry(gr3eth,"3 GeV efficiency 20mm","p");
  leg->AddEntry(gr3nth,"3 GeV noise 20mm","p");
  leg->AddEntry(gr4eth,"4 GeV efficiency 20mm","p");
  leg->AddEntry(gr4nth,"4 GeV noise 20mm","p");
  leg->AddEntry(gr5eth,"5 GeV efficiency 20mm","p");
  leg->AddEntry(gr5nth,"5 GeV noise 20mm","p");
  leg->AddEntry(gr6eth,"6 GeV efficiency 20mm","p");
  leg->AddEntry(gr6nth,"6 GeV noise 20mm","p");

  leg->Draw();

  // Output
  threshold->Print("pics/noiseeffi.eps");
  _outputFile->cd();
  threshold->Write();
  threshold->Close();

}*/

/*
// Do a measured resolutions vs E plot
if(runmode == 5)
{

std::cout << " " << std::endl;
std::cout << "Mode 5" << std::endl;
std::cout << " " << std::endl;

Double_t resolution[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
Double_t resolutionerror[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
Double_t energy[7] = { 2, 3, 4.4, 5, 2, 3, 5};
Double_t energyerror[7] = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

std::vector<int> energyruns;

energyruns.push_back(5043);
energyruns.push_back(5063);
energyruns.push_back(98);
energyruns.push_back(37);
energyruns.push_back(236);
energyruns.push_back(249);
energyruns.push_back(262);

for (int j=0;j<7;j++)
{

int runnumber = energyruns.at(j);

std::string help;
if(runnumber < 113 || runnumber > 700) help = inputDir150 + inputFile;
else help = inputDir20 + inputFile;
TString filename(help.c_str());
if (runnumber <= 999)
filename+="0";
if (runnumber <= 99)
filename+="0";
filename+=runnumber;
filename+="-fitter.root";
TFile *tfile = new TFile(filename);
//TH1D *h_m26[planescount];
tfile->cd();


TH1D *h_m26[2];

// Load file

h_m26[0] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftX"))->Clone();
h_m26[1] = (TH1D*) (tfile->Get("DUTHisto03/DUTshiftY"))->Clone();

float temp1 = 0.0;
float temp2 = 0.0;

for (int i=0;i<2;i++)
{
TF1 *f1 = new TF1("f1","gaus");
h_m26[i]->Fit(f1,"EMQI0","");
Double_t sigma_1 = f1->GetParameter(2);
Double_t error_1 = f1->GetParError(2);

temp1 += sigma_1;
temp2 += error_1;

std::cout << "mean 1 " << sigma_1 << std::endl;
}

temp1 = temp1 / 2.0;
temp2 = temp2 / 2.0;

resolution[j] = temp1*1000;
resolutionerror[j] = temp2*1000;

std::cout << "resolution is " << resolution[j] << " at j " << j << std::endl;

}
// Create graphs with the information
TGraphErrors *gr = new TGraphErrors(7,energy,resolution,energyerror,resolutionerror);

// Let's plot this
TH1D *h_axis = new TH1D("th_axis","th_axis",1, 0.0, 10.0);
TCanvas *energyplot = new TCanvas("energyplot","energyplot",10,10,800,600);
gStyle->SetPadBorderMode(0);
gStyle->SetOptStat(0);
energyplot->SetFillColor(0);
energyplot->Divide(1,1);
energyplot->cd(1);
gStyle->SetErrorX(0);
gPad->SetLogx(0);
gPad->SetLogy(0);

//histo_cfg(h_axis, "Threshold (s/n)","N","");
h_axis->SetMinimum(0.0);
h_axis->SetMaximum(10.0);
h_axis->Draw("hist");

gr->Draw("P");

energyplot->Print("pics/energyplot.eps");
_outputFile->cd();
energyplot->Write();
energyplot->Close();

}*/

// FIXME update
// plot clustersize
if(runmode == 6)
{

  // Results go in here
  Double_t threshcluster2[threshcount];
  Double_t threshcluster3[threshcount];
  Double_t threshcluster4[threshcount];
  Double_t threshcluster5[threshcount];
  Double_t threshclusterthreshcouunt0[threshcount];
  Double_t threshcluster2_error[threshcount];
  Double_t threshcluster3_error[threshcount];
  Double_t threshcluster4_error[threshcount];
  Double_t threshcluster5_error[threshcount];
  Double_t threshcluster120[threshcount];
  Double_t threshcluster120_error[threshcount];
  Double_t thin_threshcluster2[threshcount];
  Double_t thin_threshcluster3[threshcount];
  Double_t thin_threshcluster5[threshcount];
  Double_t thin_threshcluster2_error[threshcount];
  Double_t thin_threshcluster3_error[threshcount];
  Double_t thin_threshcluster5_error[threshcount];
  Double_t x[threshcount];
  Double_t xerror[threshcount] = {0.0};

  std::cout << " " << std::endl;
  std::cout << "Mode 6" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Running over all runs - cluster size" << std::endl;
  std::cout << " " << std::endl;
  std::cout << "Wide Geometry" << std::endl;
  std::cout << " " << std::endl;

  telescopebuild = "AnaTel_wide.geom";
  planedistance = 150;
  for(int j=0;j<nplanes;j++)
    posx[j] = 150.0*j;

  for(int i=0;i<run_2_150.size();i++)
  {
    getclusize( run_2_150[i] );
    x[i] = i+3;
    threshcluster2[i] = avgclustersize;
    threshcluster2_error[i] = avgclustersize_error;
  }

  for(int i=0;i<run_3_150.size();i++)
  {
    getclusize( run_3_150[i] );
    x[i] = i+3;
    threshcluster3[i] = avgclustersize;
    threshcluster3_error[i] = avgclustersize_error;
  }

  for(int i=0;i<run_4_150.size();i++)
  {
    getclusize( run_4_150[i] );
    x[i] = i+3;
    threshcluster4[i] = avgclustersize;
    threshcluster4_error[i] = avgclustersize_error;
  }

  for(int i=0;i<run_5_150.size();i++)
  {
    getclusize( run_5_150[i] );
    x[i] = i+3;
    threshcluster5[i] = avgclustersize;
    threshcluster5_error[i] = avgclustersize_error;
  }

  for(int i=0;i<run120.size();i++)
  {
    getclusize( run120[i] );
    x[i] = i+3;
    threshcluster120[i] = avgclustersize;
    threshcluster120_error[i] = avgclustersize_error;
  }

  std::cout << " " << std::endl;
  std::cout << "Thin Geometry" << std::endl;
  std::cout << " " << std::endl;

  telescopebuild = "AnaTel_narrow.geom";
  planedistance = 20;

  for(int j=0;j<nplanes;j++)
    posx[j] = 20.0*j;

  for(int i=0;i<run_2_20.size();i++)
  {
    getclusize( run_2_20[i] );
    x[i] = i+3;
    thin_threshcluster2[i] = avgclustersize;
    thin_threshcluster2_error[i] = avgclustersize_error;
  }

  for(int i=0;i<run_3_20.size();i++)
  {
    getclusize( run_3_20[i] );
    x[i] = i+3;
    thin_threshcluster3[i] = avgclustersize;
    thin_threshcluster3_error[i] = avgclustersize_error;
  }

  for(int i=0;i<run_5_20.size();i++)
  {
    getclusize( run_5_20[i] );
    x[i] = i+3;
    thin_threshcluster5[i] = avgclustersize;
    thin_threshcluster5_error[i] = avgclustersize_error;
  }

  // Create graphs with the information
  TGraphErrors *gr2n = new TGraphErrors(10,x,threshcluster2,xerror,threshcluster2_error);
  TGraphErrors *gr3n = new TGraphErrors(10,x,threshcluster3,xerror,threshcluster3_error);
  TGraphErrors *gr4n = new TGraphErrors(10,x,threshcluster4,xerror,threshcluster4_error);
  TGraphErrors *gr5n = new TGraphErrors(10,x,threshcluster5,xerror,threshcluster5_error);
  TGraphErrors *gr120n = new TGraphErrors(10,x,threshcluster120,xerror,threshcluster120_error);

  TGraphErrors *gr2nth = new TGraphErrors(10,x,thin_threshcluster2,xerror,thin_threshcluster2_error);
  TGraphErrors *gr3nth = new TGraphErrors(10,x,thin_threshcluster3,xerror,thin_threshcluster3_error);
  TGraphErrors *gr5nth = new TGraphErrors(10,x,thin_threshcluster5,xerror,thin_threshcluster5_error);

  // Let's plot this
  TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
  TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);
  threshold->SetFillColor(0);
  threshold->Divide(1,1);
  threshold->cd(1);
  gStyle->SetErrorX(0);
  gPad->SetLogx(0);
  gPad->SetLogy(0);

  // Set apperance
  gr2n->SetMarkerStyle(22);
  gr2n->SetMarkerColor(kRed);
  gr2n->SetMarkerSize(2);
  gr3n->SetMarkerStyle(22);
  gr3n->SetMarkerColor(kBlue);
  gr3n->SetMarkerSize(2);
  gr4n->SetMarkerStyle(22);
  gr4n->SetMarkerColor(kGreen);
  gr4n->SetMarkerSize(2);
  gr5n->SetMarkerStyle(22);
  gr5n->SetMarkerColor(kBlack);
  gr5n->SetMarkerSize(2);
  gr120n->SetMarkerStyle(22);
  gr120n->SetMarkerColor(kOrange);
  gr120n->SetMarkerSize(2);

  gr2nth->SetMarkerStyle(34);
  gr2nth->SetMarkerColor(kRed);
  gr2nth->SetMarkerSize(2);
  gr3nth->SetMarkerStyle(34);
  gr3nth->SetMarkerColor(kBlue);
  gr3nth->SetMarkerSize(2);
  gr5nth->SetMarkerStyle(34);
  gr5nth->SetMarkerColor(kBlack);
  gr5nth->SetMarkerSize(2);

  //histo_cfg(h_axis, "Threshold (s/n)","N","");
  h_axis->SetMinimum(0.0);
  h_axis->SetMaximum(3.0);
  h_axis->Draw("hist");

  gr2n->Draw("P");
  gr3n->Draw("P");
  gr4n->Draw("P");
  gr5n->Draw("P");
  gr120n->Draw("P");

  gr2nth->Draw("P");
  gr3nth->Draw("P");
  gr5nth->Draw("P");

  // The legend
  TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetHeader("Performance:");

  leg->AddEntry(gr2n,"2 GeV avg. cluster size","p");
  leg->AddEntry(gr3n,"3 GeV avg. cluster size","p");
  leg->AddEntry(gr4n,"4.4 GeV avg. cluster size","p");
  leg->AddEntry(gr5n,"5 GeV avg. cluster size","p");
  leg->AddEntry(gr120n,"120 GeV avg. cluster size","p");

  leg->AddEntry(gr2nth,"2 GeV avg. cluster size 20mm","p");
  leg->AddEntry(gr3nth,"3 GeV avg. cluster size 20mm","p");
  leg->AddEntry(gr5nth,"5 GeV avg. cluster size 20mm","p");

  leg->Draw();

  // Output
  threshold->Print("pics/clusize.eps");
  _outputFile->cd();
  threshold->Write();
  threshold->Close();

}

/*
// plot pointing
if(runmode == 7)
{

// Results go in here
Double_t threshcluster2[threshcount];
Double_t threshcluster2_error[threshcount];
Double_t threshcluster2w[threshcount];
Double_t threshcluster2w_error[threshcount];
Double_t x[threshcount];
Double_t xerror[threshcount] = {0.0};
Double_t xw[threshcount];
Double_t xwerror[threshcount] = {0.0};

std::cout << " " << std::endl;
std::cout << "Mode 7" << std::endl;
std::cout << " " << std::endl;

double thresh = 6-3;
// x -3 with x the threshold

getpointing( run_2_150[thresh] , 3.415567);
x[0] = 2;
xerror[0] = 0.2;
threshcluster2[0] = avgmeas;
threshcluster2_error[0] = avgmeas_error;

getpointing( run_3_150[thresh] , 3.468273);
x[1] = 3;
xerror[1] = 0.3;
threshcluster2[1] = avgmeas;
threshcluster2_error[1] = avgmeas_error;

getpointing( run_4_150[thresh] , 3.530031);
x[2] = 4.4;
xerror[2] = 0.44;
threshcluster2[2] = avgmeas;
threshcluster2_error[2] = avgmeas_error;

getpointing( run_5_150[thresh] , 3.444235);
x[3] = 5;
xerror[3] = 0.5;
threshcluster2[3] = avgmeas;
threshcluster2_error[3] = avgmeas_error;


// 3.194782 thr 6
// 3.17949 thr 7
// 3.31353 thr 8
// 3.4976 thr 9 
getpointing( run120[thresh+2] , 3.31353);
x[4] = 120;
xerror[4] = 12;
threshcluster2[4] = avgmeas;
threshcluster2_error[4] = avgmeas_error;

getpointing( run_2_20[thresh] , 3.490106);
xw[0] = 2;
xwerror[0] = 0.2;
threshcluster2w[0] = avgmeas;
threshcluster2w_error[0] = avgmeas_error;

getpointing( run_3_20[thresh] , 3.445178);
xw[1] = 3;
xwerror[1] = 0.3;
threshcluster2w[1] = avgmeas;
threshcluster2w_error[1] = avgmeas_error;

getpointing( run_5_20[thresh] , 3.420506);
xw[2] = 5;
xwerror[2] = 0.5;
threshcluster2w[2] = avgmeas;
threshcluster2w_error[2] = avgmeas_error;


//getpointing( testing[0] , 3.194782);

//getpointing( 293 , 3.01298);

//getpointing( 294 , 3.17294);


//getpointing( 295 , 3.3148);
//getpointing( 296 , 3.60078);

// all fixed but m26

getpointing( 293 , 2.51878);
getpointing( 294 , 2.59982);
getpointing( 295 , 2.74087);
getpointing( 296 , 2.96215);
getpointing( 297 , 3.17257);
getpointing( 298 , 3.33224);
getpointing( 299 , 3.51072);
getpointing( 300 , 3.60078);



x[5] = 12.5;
xerror[5] = 1.25;
threshcluster2[5] = avgmeas;
threshcluster2_error[5] = avgmeas_error;



// Create graphs with the information
TGraphErrors *gr2n = new TGraphErrors(6,x,threshcluster2,xerror,threshcluster2_error);
TGraphErrors *gr2w = new TGraphErrors(3,xw,threshcluster2w,xwerror,threshcluster2w_error);


// Let's plot this
TH1D *h_axis = new TH1D("th_axis","th_axis",100, 0.0, 130.0);
TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
gStyle->SetPadBorderMode(0);
gStyle->SetOptStat(0);
threshold->SetFillColor(0);
threshold->Divide(1,1);
threshold->cd(1);
gStyle->SetErrorX(0);
gPad->SetLogx(0);
gPad->SetLogy(0);

// Set apperance
gr2n->SetMarkerStyle(22);
gr2n->SetMarkerColor(kRed);
gr2n->SetMarkerSize(3);
gr2w->SetMarkerStyle(22);
gr2w->SetMarkerColor(kBlue);
gr2w->SetMarkerSize(3);


//histo_cfg(h_axis, "Beam Energy in GeV","#sigma_{Point}","");
h_axis->SetMinimum(0.0);
h_axis->SetMaximum(10.0);
h_axis->Draw("hist");

gr2n->Draw("P");
gr2w->Draw("P");


// The legend
TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
// leg->SetHeader("Performance:");

leg->AddEntry(gr2n,"#Delta_{z} = 150mm","p");
leg->AddEntry(gr2w,"#Delta_{z} = 20mm","p");
leg->Draw();

// Output
threshold->Print("pics/pointing.eps");
_outputFile->cd();
threshold->Write();
threshold->Close();

}*/

/*if(runmode == 8)
  {
  std::cout << "test mode" << std::endl;

// histoplot(37);

for (int i = 0;i<10;i++)
{
histoplot (run_2_150[i]);
histoplot (run_3_150[i]);
histoplot (run_4_150[i]);
histoplot (run_5_150[i]);
histoplot (run120[i]);
histoplot (run_5_20[i]);
histoplot (run_3_20[i]);
histoplot (run_2_20[i]);

}


_outputFile->cd();
delta0->Write();
delta1->Write();
delta2->Write();
delta3->Write();
delta4->Write();
delta5->Write();
delta6->Write();
delta7->Write();
delta8->Write();
delta9->Write();
delta10->Write();
delta11->Write();


}*/ 


_outputFile->Close();

// And we're done
std::cout << "\nDone.\n" << std::endl;

}
