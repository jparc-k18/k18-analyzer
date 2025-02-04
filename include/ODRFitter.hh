#ifndef ODR_FITTER_HH
#define ODR_FITTER_HH

#include <TVector3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TMath.h>
#include <iostream>

class ODRFitter
{
public:
  ODRFitter();
  ~ODRFitter();

  void Reset();

  void SetCentroid(Double_t x, Double_t y, Double_t z);
  void SetCentroid(TVector3 pos);

  void AddPoint(Double_t x, Double_t y, Double_t z, Double_t w = 1);

  void SetMatrixA(Double_t a00, Double_t a01, Double_t a02,
		  Double_t a11, Double_t a12, Double_t a22);

  void SetWeightSum(Double_t val);
  void SetNumPoints(Int_t n);

  Bool_t FitPlane();
  Bool_t FitLine();

  Bool_t Solve();
  void ChooseEigenValue(Int_t n);

  TVector3 GetCentroid() {return TVector3(fXCentroid, fYCentroid, fZCentroid);}
  TVector3 GetNormal() {return TVector3((*fNormal)[0], (*fNormal)[1], (*fNormal)[2]);}
  TVector3 GetDirection() {return TVector3((*fNormal)[0], (*fNormal)[1], (*fNormal)[2]);}

  Int_t GetNumPoints() {return fNumPoints;}
  Double_t GetWeightSum() {return fWeightSum;}
  Double_t GetRMSLine() {return fRMSLine;}
  Double_t GetRMSPlane() {return fRMSPlane;}

private:

  Int_t fNumPoints;
  Double_t fWeightSum;
  Double_t fSumOfPC2;

  Double_t fXCentroid;
  Double_t fYCentroid;
  Double_t fZCentroid;

  TVectorD *fNormal;

  TMatrixD *fMatrixA;
  TVectorD *fEigenValues;
  TMatrixD *fEigenVectors;

  Double_t fRMSLine;
  Double_t fRMSPlane;

};

#endif
  
