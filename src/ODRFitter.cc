#include "ODRFitter.hh"
#include <iostream>

ODRFitter::ODRFitter() {
  fNormal = new TVectorD(3);
  fMatrixA = new TMatrixD(3,3);
  fEigenValues = new TVectorD(3);
  fEigenVectors = new TMatrixD(3,3);

  Reset();
}

ODRFitter::~ODRFitter() {}

void ODRFitter::Reset() {
  fXCentroid = 0;  fYCentroid = 0;  fZCentroid = 0;

  fNumPoints = 0;
  fWeightSum = 0;
  fSumOfPC2 = 0;

  for (int i = 0 ; i < 3 ; i++) 
    for (int j = 0 ; j < 3 ; j++) 
      (*fMatrixA)[i][j] = 0;
  
  fRMSLine = -1;
  fRMSPlane = -1;

}

void ODRFitter::SetCentroid(Double_t x, Double_t y, Double_t z){

  fXCentroid = x;
  fYCentroid = y;
  fZCentroid = z;
}

void ODRFitter::AddPoint(Double_t x, Double_t y, Double_t z, Double_t w) {

  Double_t dx = x - fXCentroid;
  Double_t dy = y - fYCentroid;
  Double_t dz = z - fZCentroid;

  Double_t wx2 = w * dx * dx;
  Double_t wy2 = w * dy * dy;
  Double_t wz2 = w * dz * dz;

  (*fMatrixA)[0][0] += wx2;
  (*fMatrixA)[0][1] += w * dx * dy;
  (*fMatrixA)[0][2] += w * dz * dx;

  (*fMatrixA)[1][1] += wy2;
  (*fMatrixA)[1][2] += w * dy * dz;
  
  (*fMatrixA)[2][2] += wz2;

  fSumOfPC2 += wx2 + wy2 + wz2;
  fWeightSum += w;
  fNumPoints++;
}

void ODRFitter::SetMatrixA(Double_t a00, Double_t a01, Double_t a02,
			   Double_t a11, Double_t a12, Double_t a22){

  (*fMatrixA)[0][0] = a00;
  (*fMatrixA)[0][1] = a01;
  (*fMatrixA)[0][2] = a02;

  (*fMatrixA)[1][1] = a11;
  (*fMatrixA)[1][2] = a12;
  
  (*fMatrixA)[2][2] = a22;

  fSumOfPC2 += a00 + a11 + a22;
}

void ODRFitter::SetWeightSum(Double_t weightSum) {fWeightSum = weightSum;}
void ODRFitter::SetNumPoints(Int_t numPoints) {fNumPoints = numPoints;}

Bool_t ODRFitter::Solve(){

  (*fMatrixA)[1][0] = (*fMatrixA)[0][1];
  (*fMatrixA)[2][0] = (*fMatrixA)[0][2];
  (*fMatrixA)[2][1] = (*fMatrixA)[1][2];
  
  if ((*fMatrixA)[0][0] == 0 && (*fMatrixA)[1][1] == 0 &&(*fMatrixA)[2][2] == 0)
    return false;

  (*fEigenVectors) = fMatrixA->EigenVectors(*fEigenValues);

  return true;
}

void ODRFitter::ChooseEigenValue(Int_t n){
  (*fNormal) = TMatrixDColumn((*fEigenVectors),n);

  fRMSLine = (fSumOfPC2 - (*fEigenValues)[n]) / (fWeightSum - 2*fWeightSum/fNumPoints);
  fRMSLine = TMath::Sqrt(fRMSLine);

  fRMSPlane = (*fEigenValues)[n] / (fWeightSum - 2*fWeightSum/fNumPoints);
  if (fRMSPlane < 0) fRMSPlane = 0;
  fRMSPlane = TMath::Sqrt(fRMSPlane);
}

Bool_t ODRFitter::FitLine() {

  if (Solve() == false)
    return false;

  ChooseEigenValue(0);
  return true;
}

Bool_t ODRFitter::FitPlane() {

  if (Solve() == false)
    return false;

  ChooseEigenValue(2);
  return true;
}





  

  
