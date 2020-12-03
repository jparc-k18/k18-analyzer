// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>

// Counters ___________________________________________________________
const int DetIdBH1      =  1;
const int DetIdBH2      =  2;
const int DetIdBAC      =  3;
const int DetIdPVAC     =  4;
const int DetIdFAC      =  5;
const int DetIdSCH      =  6;
const int DetIdTOF      =  7;


const int DetIdLAC      = 12;
const int DetIdWC       = 13;
const int NumOfSegBH1   = 11;
const int NumOfSegBH2   =  5;
const int NumOfSegBAC   =  2;
const int NumOfSegSCH   = 64;
const int NumOfSegTOF   = 24;
const int NumOfSegPVAC  =  1;
const int NumOfSegFAC   =  1;
const int NumOfSegLAC   = 15;
const int NumOfSegWC    =  2; //for E40 data analysis

// Misc _______________________________________________________________
const int DetIdTrig       = 21;
const int DetIdScaler     = 22;
const int DetIdMsT        = 25;
const int DetIdMtx        = 26;
const int DetIdFpgaBH2Mt  = 29;
const int DetIdVmeRm      = 81;
const int DetIdMsTRM      = 82;
const int DetIdHulRM      = 83;
const int NumOfSegTrig    = 32;
const int NumOfSegScaler  = 96;
const int SpillEndFlag    = 27; // 0-based
const int NumOfPlaneVmeRm = 2;

enum eTriggerFlag
  {
    kBh21K      =  0,
    kBh22K      =  1,
    kBh23K      =  2,
    kBh24K      =  3,
    kBh25K      =  4,
    kBh26K      =  5,
    kBh27K      =  6,
    kBh28K      =  7,
    kBh2K       =  8,
    kElseOr     =  9,
    kBeam       = 10,
    kBeamTof    = 11,
    kBeamPi     = 12,
    kBeamP      = 13,
    kCoin1      = 14,
    kCoin2      = 15,
    kE03        = 16,
    kBh2KPs     = 17,
    kBeamPs     = 18,
    kBeamTofPs  = 19,
    kBeamPiPs   = 10,
    kBeamPPs    = 21,
    kCoin1Ps    = 22,
    kCoin2Ps    = 23,
    kE03Ps      = 24,
    kClock      = 25,
    kReserve2   = 26,
    kSpillEnd   = 27,
    kMatrix     = 28,
    kMstAccept  = 29,
    kMstClear   = 30,
    kTofTiming  = 31
  };

const int DetIdVmeCalib      = 999;
const int NumOfPlaneVmeCalib =   5;
const int NumOfSegVmeCalib   =  32;

// Trackers ___________________________________________________________
const int DetIdBC3  = 103;
const int DetIdBC4  = 104;
const int DetIdSDC1 = 105;
const int DetIdSDC2 = 106;
const int DetIdSDC3 = 107;
const int DetIdSDC4 = 108;
const int DetIdBFT  = 110;

const int PlMinBcIn        =   1;
const int PlMaxBcIn        =  12;
const int PlMinBcOut       =  13;
const int PlMaxBcOut       =  24;
const int PlMinSdcIn       =   1;
const int PlMaxSdcIn       =  10;
const int PlMinSdcOut      =  31;
const int PlMaxSdcOut      =  38;
const int PlOffsBc         = 100;
const int PlOffsSdcIn      =   0;
const int PlOffsSdcOut     =  30;
const int PlOffsVP         =  20;
const int PlOffsFht        =  80;

const int NumOfLayersBc     = 6;
const int NumOfLayersSDC1   = 6;
const int NumOfLayersSDC2   = 4;
const int NumOfLayersSDC3   = 4;
const int NumOfLayersSDC4   = 4;
const int NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const int NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;
const int NumOfLayersVP     = 5;

const int MaxWireBC3      =  64;
const int MaxWireBC4      =  64;

const int MaxWireSDC1     =  64;
const int MaxWireSDC2X    =  70;
const int MaxWireSDC2Y    =  40;
const int MaxWireSDC3     = 128;
const int MaxWireSDC4X    =  96;
const int MaxWireSDC4Y    =  64;

const int NumOfPlaneBFT   =   2;
const int NumOfSegBFT     = 160;

// HulRm -----------------------------------------------
const int NumOfHulRm   = 4;

// Matrix ----------------------------------------------
const int NumOfSegSFT_Mtx = 48;

// MsT -------------------------------------------------
enum TypesMst{typeHrTdc, typeLrTdc, typeFlag, NumOfTypesMst};
const int NumOfMstHrTdc = 32;
const int NumOfMstLrTdc = 64;
const int NumOfMstFlag  = 7;
enum dTypesMst
  {
    mstClear,
    mstAccept,
    finalClear,
    cosolationAccept,
    fastClear,
    level2,
    noDecision,
    size_dTypsMsT
  };

// Scaler ----------------------------------------------
const int NumOfScaler  = 2;

#endif
