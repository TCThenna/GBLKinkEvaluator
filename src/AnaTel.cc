/*   Exchange pointing resolution estimate by GBL calculation

Author: H. Jansen
email hendrik.jansen@desy.de
Date: 3.8.15
*/

#include "GblTrajectory.h"
#include "GblPoint.h"
#include "GblData.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "AnaTel.h"
#include "TMatrixD.h"
#include "TVectorD.h"
//#include "TF1.h"

TMatrixD Jac5( double ds )
{
  /*
     straight line, no B-field
     track = 
     q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  TMatrixD jac(5, 5);
  jac.UnitMatrix();
  jac[3][1] = ds; // x = xp * ds
  jac[4][2] = ds; // y = yp * ds
  return jac;
}

gbl::GblPoint AnaTel::getPoint(double step, double res, TVectorD wscat, bool IsPlane, bool has_meas) {

  // Propagate:
  TMatrixD jacPointToPoint = Jac5(step);
  gbl::GblPoint point(jacPointToPoint);

  // Add scatterer:
  TVectorD scat(2);
  scat.Zero(); // mean is zero
  point.addScatterer(scat, wscat);

  // Add measurement if requested:
  if(has_meas) {
    // measurement = residual
    TVectorD meas(2);
    meas.Zero(); // ideal

    // Precision = 1/resolution^2
    TVectorD measPrec(2);
    measPrec[0] = 1.0 / res / res;
    measPrec[1] = 1.0 / res / res;

    TMatrixD proL2m(2,2);
    proL2m.UnitMatrix();

    point.addMeasurement(proL2m, meas, measPrec);

  }


  _s += step;
  _sPoint.push_back(_s);
  if(IsPlane) _ID.push_back(_sPoint.size());

  return point;
}

gbl::GblPoint AnaTel::getPoint(double step, TVectorD wscat, bool IsPlane) {
  // This does not add a measurement - no resultion is given!
  return AnaTel::getPoint(step, 0.0, wscat, IsPlane, false);
}

AnaTel::AnaTel(std::string GeomFile, double _eBeam)
{

  AnaTel::SetBeam(_eBeam, 0.0);
  ifstream geometryFile;

  std::string help = "../geometries/";
  help += GeomFile;
  //std::cout << " Reading geometry file from: " << help << std::endl;
  geometryFile.open(help.c_str());

  //cout << "Reading telescope geometry description from " << GeomFile << endl;

  geometryFile >> _nTelPlanes  >> _iDUT ;

  // !!! Important: DUT position in geometry file numbered 0...N-1 !!!
  //     _iDUT = -1 denotes no DUT in the Setup

  _planePosition   = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

  _sumeps = 0.0;

  // Set kappa = 1 (i.e. use Highland as is)
  _kappa = 1.0;

  // Planes have to be ordered in position along the beam line !
  // This is not checked !!!

  //std::cout << " nTelPlanes = " << _nTelPlanes << std::endl;

  for(int ipl=0; ipl < _nTelPlanes; ipl++)
  {
    int iActive;
    double resolution;

    // All dimensions should be given in mm !!!

    geometryFile >> _planePosition[ipl]
      >> _planeThickness[ipl]
      >> _planeX0[ipl]
      >> iActive 
      >> resolution ;

    if(iActive)
    {
      _isActive[ipl] = true ;
      _planeResolution[ipl]=resolution;
    }
    else
    {
      _isActive[ipl] = false ;
      _planeResolution[ipl]=0.;
    }
  }

  // Print out geometry information
  //cout << "Telescope configuration with " << _nTelPlanes << " planes" << endl;


  for(int ipl=0; ipl < _nTelPlanes; ipl++)
  {
    std::stringstream ss ; 

    if ( ipl == _iDUT)
    {
      if(_isActive[ipl])
	ss << ipl << " : active  DUT   at" ;
      else
	ss << ipl << " : passive  DUT  at" ; 
    }
    else
    {
      if(_isActive[ipl])
	ss << ipl << " : active  plane at" ;
      else
	ss << ipl << " : passive plane at" ; 
    }


    ss << "  Z [mm] = " << _planePosition[ipl] 
      << " dZ [um] = " << _planeThickness[ipl]*1000. ;

    if(_isActive[ipl])
      ss << "  Res [um] = " << _planeResolution[ipl]*1000. ;

    //cout << ss.str()  << endl;
  }
  // Allocate arrays for track fitting

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];
  _tempScat1 = new double[_nTelPlanes];
  _tempScat2 = new double[_nTelPlanes];
  _tempScat3 = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _pointingResolution  = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];


  if(_iDUT>=0)
    _useDUT=true;
  else
    _useDUT=false;

  //_eBeam=0.; // not needed anymore as constructor requires to pass the energy already
  _useBeamConstraint=false ;

  // GBL

  _s = 0;
  _sPoint.clear();
  _listOfPoints.clear();

}


AnaTel::AnaTel(Int_t Npl)
{
  _nTelPlanes=Npl;

  _planePosition   = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];
  _tempScat1 = new double[_nTelPlanes];
  _tempScat2 = new double[_nTelPlanes];
  _tempScat3 = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _pointingResolution  = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];

  _iDUT=-1;
  _useDUT=false;

  _useBeamConstraint=false ;
}




void AnaTel::SetPlane(Int_t Ipl, Double_t Position,  Double_t Thickness,   Double_t X0, Double_t Resolution)
{

  _planePosition[Ipl]=Position;
  _planeThickness[Ipl]=Thickness;
  _planeX0[Ipl]=X0;

  if(Resolution>0)
  {
    _isActive[Ipl] = true ;
    _planeResolution[Ipl]=Resolution;
  }
  else
  {
    _isActive[Ipl] = false ;
    _planeResolution[Ipl]=0.;
  }


  std::stringstream ss ; 
  if(_isActive[Ipl])
    ss << "Active  plane at" ;
  else
    ss << "Passive plane at" ; 

  ss << "  Z [mm] = " << _planePosition[Ipl] 
    << " plane thickness [um] = " << _planeThickness[Ipl]*1000. ;

  if(_isActive[Ipl])
    ss << "  Res [um] = " << _planeResolution[Ipl]*1000. ;

  //std::cout << ss.str()  << std::endl;
} 


void  AnaTel::SetDUT(Int_t idut, Bool_t UseInFit)
{
  _iDUT=idut;
  _useDUT=UseInFit;

  if (_iDUT >=0 && _iDUT <  _nTelPlanes) 
  {
    /*if(_isActive[_iDUT] && UseInFit)
      std::cout << "Active plane " << idut << "  set as DUT, " << "used in fit" << std::endl;
    else
      if(_isActive[_iDUT])
	std::cout << "Active plane " << idut << "  set as DUT, " << "not used in fit" << std::endl;
      else
	std::cout  << "Passive plane " << idut << "  set as DUT "  << std::endl;
    */


    if(!_isActive[_iDUT] && UseInFit)
    {
      std::cout << "Warning: passive plane can not be used in fit "  << std::endl;
      _useDUT=false;
    }

  }
  else
    std::cout << "Setup without DUT "  << std::endl;

}

void  AnaTel::SetBeam(Double_t Energy, Double_t Spread)
{
  _eBeam = Energy;

  if(Spread>0.)
  {
    _useBeamConstraint= true ;
    _beamSpread = Spread;
  }
  else
  {
    _useBeamConstraint= false ;
    _beamSpread = 0.;
  }

}

void AnaTel::SetResolution(Int_t Ipl, Double_t res)
{
  _planeResolution[Ipl] = res;
}


void AnaTel::SetResolution(Double_t Resolution)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
    _planeResolution[ipl]= Resolution;
}

void AnaTel::SetResolution(Double_t * Resolution)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
    _planeResolution[ipl]= Resolution[ipl];
}

void AnaTel::SetThickness(Double_t thickness)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
    _planeThickness[ipl]= thickness;
}

void AnaTel::SetKappa(Double_t kappa)
{
  _kappa = kappa;
}

Int_t AnaTel::GetNplanes()
{
  return _nTelPlanes;
}

Double_t * AnaTel::GetPosition()
{
  return _planePosition;
}

Double_t * AnaTel::GetThickness()
{
  return _planeThickness;
}

Double_t * AnaTel::GetResolution()
{
  return _planeResolution;
}


Double_t AnaTel::GetPointingResGBL(Int_t Ipl, Bool_t UseInFit)
{
  
  std::vector<double> pointingResolutionLoc;

  gbl::GblTrajectory traj = BuildTscope(Ipl, UseInFit);
  TVectorD aCorrection(5);
  TMatrixDSym aCovariance(5);

  for(int iipl = 0; iipl < _nTelPlanes; iipl++){
    traj.getResults( _ID[iipl], aCorrection, aCovariance );

    pointingResolutionLoc.push_back(sqrt(aCovariance(3,3)));
    //std::cout << iipl << "  " << _ID[iipl] <<" pointingRes = " << pointingResolutionLoc.back() << ",  ";
  }
  //std::cout << std::endl;
  //


  return pointingResolutionLoc.at(Ipl);
}

Double_t AnaTel::GetRhat(Int_t Ipl, Bool_t UseInFit)
{
  
  std::vector<double> rhat;

  gbl::GblTrajectory traj = BuildTscope(Ipl, UseInFit);

  unsigned int ndim = 2;
  TVectorD aResiduals(ndim);
  TVectorD aMeasErrors(ndim);
  TVectorD aResErrors(ndim);
  TVectorD aDownWeights(ndim);

  for(int iipl = 0; iipl < _nTelPlanes; iipl++){
    traj.getMeasResults( _ID[iipl], ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );

    rhat.push_back(aResErrors[0]);
    //std::cout << iipl << "  " << _ID[iipl] <<" rhat = " << rhat.back() << ",  ";
  }
  //std::cout << std::endl;
  //


  return rhat.at(Ipl);
}

Double_t AnaTel::GetkResError(Int_t Ipl, Bool_t UseInFit)
{
  
  std::vector<double> v_kinkError;

  gbl::GblTrajectory traj = BuildTscope(Ipl, UseInFit);
  TVectorD aCorrection(5);
  TMatrixDSym aCovariance(5);

  unsigned int ndim = 2;
  TVectorD aResiduals(ndim);
  TVectorD aMeasErrors(ndim);
  TVectorD aResErrors(ndim);
  TVectorD aDownWeights(ndim);


  TVectorD aKinks(ndim);
  TVectorD aKinkErrors(ndim);
  TVectorD kResErrors(ndim);
  TVectorD kDownWeights(ndim);

  for(int iipl = 0; iipl < _nTelPlanes; iipl++){
    traj.getResults( _ID[iipl], aCorrection, aCovariance );
    traj.getMeasResults( _ID[iipl], ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
    traj.getScatResults( _ID[iipl], ndim, aKinks, aKinkErrors, kResErrors, kDownWeights );

    v_kinkError.push_back(kResErrors[0]);
    //std::cout << iipl << "  " << _ID[iipl] <<" kinkError = " << v_kinkError.back() << ",  ";
  }
  //std::cout << std::endl;
  //


  return v_kinkError.at(Ipl);
}

gbl::GblTrajectory AnaTel::BuildTscope(Int_t Ipl, Bool_t UseInFit)
{

  _s = 0;
  _sPoint.clear();
  _listOfPoints.clear();

  std::cout << "Building the GBL t'scope " << std::endl;
  _listOfPoints.reserve(3*6);

  unsigned int ipl = 0;
  double step = 0;
  bool AddDummies = false;
  _s = 0.;
  if(AddDummies) _s = -10.;
  double tetSi, tetDUT, tetAir, X0Si, X0_Air_frac;

  TVectorD wscatSi(2); 
  TVectorD wscatDUT(2);
  TVectorD wscatAir(2);

  // loop over all scatterers first to calculate sumeps
  _sumeps = 0.;
  for( int ipl = 0; ipl < _nTelPlanes; ++ipl ){
    if(_planeThickness[ipl] > 0.0) {
      X0Si = _planeThickness[ipl] / _planeX0[ipl] + 0.05/ 285.6; // Si + Kapton (Kapton ist 2 * 25  = 50 um thick)
      _sumeps += X0Si;
    } else _sumeps += _planeX0[ipl];

    if( ipl < _nTelPlanes-1) {
      double distplane = _planePosition[ipl+1] - _planePosition[ipl];
      X0_Air_frac =   distplane  / 304200.; 
      _sumeps += X0_Air_frac;
    }

  }

  std::cout << "sumeps = " << _sumeps << std::endl;

  //std::cout << "E = " << _eBeam << std::endl;

  while(1)
  {
    if(_planeThickness[ipl] > 0.0) {
      _X0_Si_frac = _planeThickness[ipl] / _X0_Si + 50.e-3 / _X0_Kapton; // Si + Kapton per plane
      tetSi  = _kappa* 0.0136 * sqrt(_X0_Si_frac) / _eBeam * ( 1 + 0.038*log(_sumeps) );
    } else
      tetSi  = _kappa* 0.0136 * sqrt(_planeX0[ipl]) / _eBeam * ( 1 + 0.038*log(_sumeps) ); // for DUT (dirty hack :/ )
    //std::cout << " tet = " << tetSi << std::endl;

    wscatSi[0] = 1.0 / ( tetSi * tetSi ); // weight
    wscatSi[1] = 1.0 / ( tetSi * tetSi );

    //wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); // weight // FIXME for now w/o DUT
    //wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );

    if(ipl == 0 && AddDummies){
      // one dummy point:

      gbl::GblPoint * point = new gbl::GblPoint( Jac5( step ) );
      _listOfPoints.push_back(*point);
      _s += step;
      _sPoint.push_back(_s);
      step = 10; // [mm]
      delete point;
    }

    if( (UseInFit || ipl != Ipl) && _planeResolution[ipl] > 0.0) {
      _listOfPoints.push_back(AnaTel::getPoint(step, _planeResolution[ipl], wscatSi, true, true)); // add plane
    }
    else  _listOfPoints.push_back(AnaTel::getPoint(step, wscatSi, true)); // add plane w/o mfalse

    if((ipl+1) == _nTelPlanes) break;
    _dz = _planePosition[ipl+1] - _planePosition[ipl];
    // Air between planes; Factor 0.5 as the air is devided into two scatterers
    _X0_Air_frac =   0.5*_dz  / _X0_Air; 
    tetAir = _kappa* 0.0136 * sqrt(_X0_Air_frac) / _eBeam * ( 1 + 0.038*log(_sumeps) );

    wscatAir[0] = 1.0 / ( tetAir * tetAir ); // weight
    wscatAir[1] = 1.0 / ( tetAir * tetAir );


    step = 0.21*_dz;
    _listOfPoints.push_back(AnaTel::getPoint(step, wscatAir, false)); // add air, not added after last plane
    step = 0.58*_dz;
    _listOfPoints.push_back(AnaTel::getPoint(step, wscatAir, false)); // add air, not added after last plane
    // Set step for distance to next plane
    step = 0.21*_dz; 

    ipl++;

  }

  // add dummy point
  if(AddDummies){
    step = 10; // [mm]
    gbl::GblPoint * point = new gbl::GblPoint( Jac5( step ) );
    _listOfPoints.push_back(*point);
    _s += step;
    _sPoint.push_back(_s);
    delete point;
  }

  // fit trajectory:

  gbl::GblTrajectory traj( _listOfPoints, 0 );

  double Chi2;
  int Ndf;
  double lostWeight;




  // debug:
  
  /*
  std::cout << " Is traj valid? " << traj.isValid() << std::endl;
  traj.printPoints();
  traj.printTrajectory();
  //traj.printData();
  
  std::cout << "traj with " << traj.getNumPoints() << " points:" << std::endl;
  for( int ipl = 0; ipl < _nTelPlanes; ++ipl ){
    std::cout << "  plane " << ipl << ", lab " << _ID[ipl];
    std::cout << "  z " << _sPoint[_ID[ipl]-1];
    std::cout << std::endl;
  }
  */
  
  //std::cout << " -- Perform GBL fit -- " << std::endl;
  traj.fit( Chi2, Ndf, lostWeight );
  

  return traj;
}

Double_t AnaTel::GetWidthGBL(Int_t Ipl, Bool_t UseInFit)
{

  if(_planeResolution[Ipl] <=0.)
  {
    std::cerr << " Can not calculate residual width for passive plane !" << std::endl ;
    return 0.;
  }

  Double_t PointingRes = GetPointingResGBL( Ipl, UseInFit);
  //std::cout << "      GBL PR for plane " << Ipl << " = " << PointingRes << std::endl;

  Double_t Width;

  //std::cout << " Current intrinsic resolution used in GetWidthGBL(): " << _planeResolution[Ipl] << std::endl;


  if(UseInFit)
    Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]-PointingRes*PointingRes);
  else
    Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]+PointingRes*PointingRes);

  return Width;
}

void AnaTel::PrintPlanes()
{
  std::cout << "ntelPlanes = " << _nTelPlanes << std::endl;
  std::cout << "E = " << _eBeam << std::endl;
  std::cout << "s = " << _s << std::endl;
  std::cout << "dz = " << _dz << std::endl;
  std::cout << " size _sPoint = " << _sPoint.size() << std::endl;
  std::cout << " size _LoP = " << _listOfPoints.size() << std::endl;

  for(unsigned int ipl = 0; ipl < _nTelPlanes; ipl++)
  {
    std::cout << ipl << " thickness = " << _planeThickness[ipl] 
                   << " position  = " << _planePosition[ipl]
		   << " planeX0   = " << _planeX0[ipl]
		   << " planeReso = " << _planeResolution[ipl] 
		   << std::endl;
  }
  
}

