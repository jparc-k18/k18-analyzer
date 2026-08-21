/*
  EvDisp.hh

  2012/1/24
*/

#ifndef EvDisp_h
#define EvDisp_h 1

#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGeometry.h"
#include "TMixture.h"
#include "TBRIK.h"
#include "TTRD1.h"
#include "TTRD2.h"
#include "TTUBS.h"
#include "TTUBS.h"
#include "TRotMatrix.h"
#include "TNode.h"
#include "TPolyMarker3D.h"
#include "TPolyLine3D.h"
#include "TView.h"
#include "TPad.h"
#include "TButton.h"
#include "TMarker3DBox.h"
#include "TPave.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TGraph.h"

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "ThreeVector.hh"
#include "DetectorID.hh"

#define MaxTrack      10

extern  TApplication *theApp;

class DCLocalTrack;

class EvDisp : public TObject {
private:
  TGeometry *gevdisp_;
  TNode     *node_;
  TPad      *tp_;
  TCanvas   *tc_;

  TGeometry *gevdisp_vtx_;
  TNode     *node_vtx_;
  TPad      *tp_vtx_;
  TCanvas   *tc_vtx_;

  TBRIK     *world_;
  TBRIK     *world_vtx_;

  /*
  TTUBE     *Target_Tube_;
  TNode     *Target_Node_;
  */
  //TNode     *Target_Node_vtx_;

  TBRIK     *Bh2Wall_Brik_;
  TNode     *Bh2Wall_Node_;

  TBRIK     *Bh2Seg_Brik_[NumOfSegBH2];
  TNode     *Bh2Seg_Node_[NumOfSegBH2];

  TBRIK     *Bh2Wall_Brik_vtx_;
  TNode     *Bh2Wall_Node_vtx_;

  TBRIK     *Bh2Seg_Brik_vtx_[NumOfSegBH2];
  TNode     *Bh2Seg_Node_vtx_[NumOfSegBH2];


  TTUBS     *Yoke1_Tubs_;
  TNode     *Yoke1_Node_;

  TTRD1     *Yoke2_Trd_;
  TNode     *Yoke2_Node_;

  TBRIK     *Yoke3_Brik_;
  TNode     *Yoke3_Node_;

  TTRD1     *Yoke4_Trd_;
  TNode     *Yoke4_Node_;

  TTRD1     *Yoke5_Trd_;
  TNode     *Yoke5_Node_;
  /*
  TNode     *Bdc3x_Node_[MaxWireBDC];
  TNode     *Bdc3xp_Node_[MaxWireBDC];
  TNode     *Bdc3u_Node_[MaxWireBDC];
  TNode     *Bdc3up_Node_[MaxWireBDC];
  TNode     *Bdc3v_Node_[MaxWireBDC];
  TNode     *Bdc3vp_Node_[MaxWireBDC];

  TNode     *Bdc3x_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3xp_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3u_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3up_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3v_Node_vtx_[MaxWireBDC];
  TNode     *Bdc3vp_Node_vtx_[MaxWireBDC];

  TNode     *Bdc4x_Node_[MaxWireBDC];
  TNode     *Bdc4xp_Node_[MaxWireBDC];
  TNode     *Bdc4u_Node_[MaxWireBDC];
  TNode     *Bdc4up_Node_[MaxWireBDC];
  TNode     *Bdc4v_Node_[MaxWireBDC];
  TNode     *Bdc4vp_Node_[MaxWireBDC];

  TNode     *Bdc4x_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4xp_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4u_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4up_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4v_Node_vtx_[MaxWireBDC];
  TNode     *Bdc4vp_Node_vtx_[MaxWireBDC];
  */

  TNode     *Sdc1u1_Node_[MaxWireSDC1];
  TNode     *Sdc1u2_Node_[MaxWireSDC1];
  TNode     *Sdc1v1_Node_[MaxWireSDC1];
  TNode     *Sdc1v2_Node_[MaxWireSDC1];

  TNode     *Sdc1u1_Node_vtx_[MaxWireSDC1];
  TNode     *Sdc1u2_Node_vtx_[MaxWireSDC1];
  TNode     *Sdc1v1_Node_vtx_[MaxWireSDC1];
  TNode     *Sdc1v2_Node_vtx_[MaxWireSDC1];

  TNode     *Sdc2v1_Node_[MaxWireSDC2];
  TNode     *Sdc2v2_Node_[MaxWireSDC2];
  TNode     *Sdc2u1_Node_[MaxWireSDC2];
  TNode     *Sdc2u2_Node_[MaxWireSDC2];
  TNode     *Sdc2x1_Node_[MaxWireSDC2];
  TNode     *Sdc2x2_Node_[MaxWireSDC2];

  TNode     *Sdc2v1_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2v2_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2u1_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2u2_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2x1_Node_vtx_[MaxWireSDC2];
  TNode     *Sdc2x2_Node_vtx_[MaxWireSDC2];
  
  TNode     *Sdc3v1_Node_[MaxWireSDC3V];
  TNode     *Sdc3x1_Node_[MaxWireSDC3X];
  TNode     *Sdc3u1_Node_[MaxWireSDC3U];
  TNode     *Sdc3v2_Node_[MaxWireSDC3V];
  TNode     *Sdc3x2_Node_[MaxWireSDC3X];
  TNode     *Sdc3u2_Node_[MaxWireSDC3U];

  TNode     *Sdc4v1_Node_[MaxWireSDC4V];
  TNode     *Sdc4x1_Node_[MaxWireSDC4X];
  TNode     *Sdc4u1_Node_[MaxWireSDC4U];
  TNode     *Sdc4v2_Node_[MaxWireSDC4V];
  TNode     *Sdc4x2_Node_[MaxWireSDC4X];
  TNode     *Sdc4u2_Node_[MaxWireSDC4U];


  TBRIK     *TofWall_Brik_;
  TNode     *TofWall_Node_;

  TBRIK     *TofSeg_Brik_;
  TNode     *TofSeg_Node_[NumOfSegTOF];
  /*
  TBRIK     *Ac1_Brik_;
  TNode     *Ac1_Node_;

  TBRIK     *Ac2_Brik_;
  TNode     *Ac2_Node_;
  */
  TBRIK     *LcWall_Brik_;
  TNode     *LcWall_Node_;

  TBRIK     *LcSeg_Brik_;
  TNode     *LcSeg_Node_[NumOfSegLC];

  mutable TPolyMarker3D   *InitStepMark_;

  //mutable TPolyLine3D     *LocalTrackBdcOut_[MaxTrack];
  mutable TPolyLine3D     *LocalTrackSdcIn_[MaxTrack];
  mutable TPolyLine3D     *LocalTrackSdcOut_[MaxTrack];
  //mutable TPolyMarker3D   *VtxPoint_[MaxTrack];

  mutable TPolyMarker3D   *SksStepMark_;

  //mutable int nTrackBdcOut_;
  mutable int nTrackSdcIn_;
  mutable int nTrackSdcOut_;
  //mutable int nVertex_;

public:
  EvDisp(void);
  ~EvDisp(void);
  static EvDisp & GetInstance( void );
  void Initialize(void);
  void ConstructTarget(void);
  void ConstructTargetVtx(void);
  void ConstructBH2(void);
  void ConstructBH2Vtx(void);
  void ConstructSKS(void);
  void ConstructBdcOut(void);
  void ConstructBdcOutVtx(void);
  void ConstructSdcIn(void);
  void ConstructSdcInVtx(void);
  void ConstructSdcOut(void);
  void ConstructTOF(void);
  void ConstructRangeCounter(void);
  void ConstructRangeCounterVtx(void);
  void DrawInitTrack(int nStep, ThreeVector *StepPoint) const;
  void DrawInitTrack(void) const;
  void DrawHitWire(int lnum, int hit_wire, bool range_check=true, bool tdc_check=true) const;
  void DrawTrackWire(int lnum, int hit_wire, int it) const;
  void DrawHitBH2(int seg, int Tu, int Td) const;
  void DrawHitHodoscope(int lnum, int seg, int Tu, int Td) const;
  void DrawHitRangeCounter(int lnum, int seg, int Tu, int Td) const;
  void DrawBdcOutLocalTrack(ThreeVector globalPos0, 
			    ThreeVector globalPos1) const;
  void DrawBdcOutLocalTrack(DCLocalTrack *tp) const;
  void DrawSdcInLocalTrack(ThreeVector globalPos0, 
			    ThreeVector globalPos1) const;
  void DrawSdcInLocalTrack(DCLocalTrack *tp) const;
  void DrawSdc1LocalTrack(DCLocalTrack *tp) const;
  void DrawSdcOutLocalTrack(ThreeVector globalPos0, 
			    ThreeVector globalPos1) const;
  void DrawSdcOutLocalTrack(DCLocalTrack *tp) const;
  void DrawVertex(ThreeVector vtxPoint, double cost) const;
  void DrawSksTrack(int nStep, ThreeVector *StepPoint) const;
  void EndOfEvent(void) const;
  void ResetVisibility() const;
  void calcRotMatrix(double TA, double RA1, double RA2, double *rotMat);
  void get_command(void) const;

private:
  static EvDisp *evDisp_;
  static TApplication *theApp;
};


#endif
