/*
  DetectorID.hh

  2012/1/24
*/

#ifndef DetectorID_h
#define DetectorID_h 1

//Detector ID
//Base chambers
const int DetIdBC1  =101;
const int DetIdBC2  =102;
const int DetIdBC3  =103;
const int DetIdBC4  =104;
const int DetIdSDC1 =105;
const int DetIdSDC2 =106;
const int DetIdHDC  =113; 
const int DetIdSDC3 =107;
const int DetIdSDC4 =108;
const int DetIdK6BDC=109;
const int DetIdBFT  =110;
const int DetIdSFT  =111;
const int DetIdSsd  =112;
//Base counters
const int DetIdBH1 =  3;
const int DetIdBH2 =  4;
const int DetIdTOF =  5;
const int DetIdAC  =  6;
const int DetIdLC  =  7;
const int DetIdGC  =  8;

const int DetIdBAC     =  9;
const int DetIdBMW     =  10;
const int DetIdMisc    =  11;
const int DetIdMatrix  =  12;
const int DetIdSFV     =  24;
const int DetIdTOFMT   =  25;


//Test counters
//const int DetIdBAC =  9;
//const int DetIdTGT =  10;
//const int DetIdBMW =  10;
//Misc
//const int DetIdMisc  =  11;
//const int DetIdMatrix=  12;

//For SksMinus
//const int DetIdSFV =  24;

//For E27
const int DetIdRC  = 21;
const int DetIdAC1  = 22;
const int DetIdBVH  = 23;

//Scaler
const int DetIdScaler = 210;


//DC Number of Plane
const int PlMinBcIn  =  1;
const int PlMaxBcIn  = 12;
const int PlMinBcOut = 13;
const int PlMaxBcOut = 24;

const int PlMinSdcIn  =  1;
const int PlMaxSdcIn  =  10;

const int PlMinSdcOut = 31;
const int PlMaxSdcOut = 42;

const int PlOffsBc = 100;
const int PlOffsSdcOut = 30;
const int PlOffsSsd = 130;

const int NumOfLayersBc     = 6;
const int NumOfLayersSdc    = 6;
const int NumOfLayersSFT    = 3;
const int NumOfLayersSdc2   = 6;
const int NumOfLayersSsd    = 2;

const int NumOfLayersBcIn  = PlMaxBcIn  - PlMinBcIn  + 1;
const int NumOfLayersBcOut = PlMaxBcOut - PlMinBcOut + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;

//BFT config
const int NumOfPlaneBFT  = 2;
const int NumOfSegBFT    = 160;

// SFT config
// SFT X layer has U D plane.
// SFT UV layers have only U plnane.
enum SFT_PLANE{SFT_X1, SFT_X2, SFT_V, SFT_U};
const int NumOfPlaneSFT  = 4;
const int NumOfSegSFT_X  = 256;
const int NumOfSegSFT_UV = 320;

const int IdK18Target = 130;

//Hodo Segments
const int NumOfSegGC  =   1;
const int NumOfSegBH1 =  11;
const int NumOfSegBH2 =   5;
const int NumOfSegTOF =  32;
const int NumOfSegAC  =  30;
const int NumOfSegAC1  =  18;
const int NumOfSegBVH  =  4;
const int NumOfSegLC   = 28;
const int NumOfSegSFV  = 18;
const int NumOfSegMisc = 16;
const int NumOfSegBMW  = 1;

const int NumOfSegBAC =  6;
const int NumOfSegTGT =   3;

const int NumOfMisc    =  20;
const int NumOfMatrix  =  10;


//for E27
const int NumOfSegRC = 42;
const int NumOfLayersRC = 7; 
const int PlaneId_RC = 6;  

const int NumOfScaler = 32;




const int pwo_type[4]={'B', 'E', 'C', 'L'};


//Number of Wires
const int MaxWireBC1  =  256;
const int MaxWireBC2  =  256;
const int MaxWireBC3  =  64;
const int MaxWireBC4  =  64;
const int MaxWireSsd  = 768;
//const int MaxWireBC4  =  64;

const int MaxWireSDC1  =  96;
const int MaxWireSDC2  =  96;
const int MaxWireHDC   =  112;

const int MaxWireSDC3   =  120;
const int MaxWireSDC3X  =  108;
const int MaxWireSDC3U  =  120;
const int MaxWireSDC3V  =  120;

const int MaxWireSDC4   =  120;
const int MaxWireSDC4X  =  108;
const int MaxWireSDC4U  =  120;
const int MaxWireSDC4V  =  120;

// E05 trigger flag bit assignment (tK18/tHodo trigflag[20]).
// Decoded from K18LogBook41.pdf p.199 (Run 10054 trigflag TDC panel
// grid, 2015.10.28): panel order matches per-bit firing rates in
// the actual data (clock=100%; KK/KPi rates match the OR terms of
// the production trigger formula on the same page: BH2*(1/100^2|1)
// + KIn*TOF/100 + (K,pi)/3 + (K,K)). Bits 18-19 are always 1 in
// data but are not labeled in the logbook -- likely spill/final-OR
// bookkeeping, not a physics condition.
enum TrigFlag
{
  TrigKIn,        // 0: K in
  TrigPiIn,       // 1: pi in
  TrigTofLcAcBar, // 2: TOFxLCxACbar
  TrigTofLc,      // 3: TOFxLC
  TrigUnbiasPS,   // 4: unbias prescaled
  TrigPiBeam,     // 5: pi beam
  TrigKPBeam,     // 6: K+p beam
  TrigPBeam,      // 7: p beam
  TrigKK,         // 8: (K,K)
  TrigKPi,        // 9: (K,pi)
  TrigKPiPS,      // 10: (K,pi) prescaled
  TrigMtxWide,    // 11: wide matrix
  TrigUnTofPS,    // 12: unbias TOF prescaled
  TrigClock,      // 13: clock (always 1)
  TrigKKPS,       // 14: (K,K) prescaled
  TrigScaler,     // 15: scaler
  TrigTofAcLc,    // 16: TOFxACxLC
  TrigMtxTight,   // 17: tight matrix
  TrigSpare18,    // 18: unlabeled, always 1 in data
  TrigSpare19,    // 19: unlabeled, always 1 in data
  NumOfSegTrig
};

#endif
