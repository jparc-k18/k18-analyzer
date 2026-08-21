/*
  EvDisp.cc

  2012/1/24
*/

//Unit is mm.

#include "EvDisp.hh"
#include "DCGeomMan.hh"

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <string>
#include <cstdlib>

#include "DCLocalTrack.hh"

const int MaxChar = 200;

const double Deg2Rad = acos(-1)/180.;
const double Rad2Deg = 180./acos(-1);

EvDisp *EvDisp::evDisp_ = 0;
TApplication *EvDisp::theApp =0;

EvDisp::EvDisp(void)
{
  InitStepMark_ = NULL;
  /*  
  for (int i=0; i<MaxTrack; i++)
    LocalTrackBdcOut_[i] = NULL;
  */
  for (int i=0; i<MaxTrack; i++)
    LocalTrackSdcIn_[i] = NULL;

  for (int i=0; i<MaxTrack; i++)
    LocalTrackSdcOut_[i] = NULL;
  /*
  for (int i=0; i<MaxTrack; i++)
    VtxPoint_[i] = NULL;
  */
  SksStepMark_ = NULL;
}

EvDisp::~EvDisp(void)
{
  if (InitStepMark_)
    delete InitStepMark_;
  /*
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBdcOut_[i])
      delete LocalTrackBdcOut_[i];
  */
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn_[i])
      delete LocalTrackSdcIn_[i];

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcOut_[i])
      delete LocalTrackSdcOut_[i];
  /*
  for (int i=0; i<MaxTrack; i++)
    if (VtxPoint_[i])
      delete VtxPoint_[i];
  */
  if (SksStepMark_)
    delete SksStepMark_;

  delete node_;
  /*
  delete Target_Node_;
  delete Target_Tube_;
  */
  delete Yoke1_Node_;
  delete Yoke1_Tubs_;

  delete Yoke2_Node_;
  delete Yoke2_Trd_;

  delete Yoke3_Node_;
  delete Yoke3_Brik_;

  delete Yoke4_Node_;
  delete Yoke4_Trd_;

  delete Yoke5_Node_;
  delete Yoke5_Trd_;
  /*
  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3x_Node_[i];
    delete Bdc3x_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3xp_Node_[i];
    delete Bdc3xp_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3u_Node_[i];
    delete Bdc3u_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3up_Node_[i];
    delete Bdc3up_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc3v_Node_[i];
    delete Bdc3v_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) { 
    delete Bdc3vp_Node_[i];
    delete Bdc3vp_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4x_Node_[i];
    delete Bdc4x_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4xp_Node_[i];
    delete Bdc4xp_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4u_Node_[i];
    delete Bdc4u_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4up_Node_[i];
    delete Bdc4up_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4v_Node_[i];
    delete Bdc4v_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireBDC; i++) {
    delete Bdc4vp_Node_[i];
    delete Bdc4vp_Node_vtx_[i];
  }
  */

  for (int i=0; i<MaxWireSDC1; i++) {
    delete Sdc1u1_Node_[i];
    delete Sdc1u2_Node_[i];
    delete Sdc1v1_Node_[i];
    delete Sdc1v2_Node_[i];

    delete Sdc1u1_Node_vtx_[i];
    delete Sdc1u2_Node_vtx_[i];
    delete Sdc1v1_Node_vtx_[i];
    delete Sdc1v2_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireSDC2; i++) {
    delete Sdc2u1_Node_[i];
    delete Sdc2u2_Node_[i];
    delete Sdc2v1_Node_[i];
    delete Sdc2v2_Node_[i];
    delete Sdc2x1_Node_[i];
    delete Sdc2x2_Node_[i];

    delete Sdc2u1_Node_vtx_[i];
    delete Sdc2u2_Node_vtx_[i];
    delete Sdc2v1_Node_vtx_[i];
    delete Sdc2v2_Node_vtx_[i];
    delete Sdc2x1_Node_vtx_[i];
    delete Sdc2x2_Node_vtx_[i];
  }

  for (int i=0; i<MaxWireSDC3V; i++) {
    delete Sdc3v1_Node_[i];
    delete Sdc3v2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC3X; i++) {
    delete Sdc3x1_Node_[i];
    delete Sdc3x2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC3U; i++) {
    delete Sdc3u1_Node_[i];
    delete Sdc3u2_Node_[i];
  }

  for (int i=0; i<MaxWireSDC4V; i++) {
    delete Sdc4v1_Node_[i];
    delete Sdc4v2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC4X; i++) {
    delete Sdc4x1_Node_[i];
    delete Sdc4x2_Node_[i];
  }
  for (int i=0; i<MaxWireSDC4U; i++) {
    delete Sdc4u1_Node_[i];
    delete Sdc4u2_Node_[i];
  }


  for (int i=0; i<NumOfSegBH2; i++) {
    delete Bh2Seg_Node_[i];
    delete Bh2Seg_Brik_[i];
    delete Bh2Seg_Node_vtx_[i];
    delete Bh2Seg_Brik_vtx_[i];
  }
  delete Bh2Wall_Brik_;
  delete Bh2Wall_Node_;
  delete Bh2Wall_Brik_vtx_;
  delete Bh2Wall_Node_vtx_;

  for (int i=0; i<NumOfSegTOF; i++) 
    delete TofSeg_Node_[i];
  delete TofSeg_Brik_;

  delete TofWall_Brik_;
  delete TofWall_Node_;
  /*
  delete Ac1_Brik_;
  delete Ac1_Node_;

  delete Ac2_Brik_;
  delete Ac2_Node_;
  */
  for (int i=0; i<NumOfSegLC; i++) 
    delete LcSeg_Node_[i];
  delete LcSeg_Brik_;

  delete LcWall_Brik_;
  delete LcWall_Node_;


  delete gevdisp_;
  delete tp_;
  delete tc_;


  delete gevdisp_vtx_;
  delete tp_vtx_;
  delete tc_vtx_;

}

EvDisp & EvDisp::GetInstance( void )
{
  if( !evDisp_ ){
    evDisp_ = new EvDisp();
  }
  if( !theApp ){
    theApp=new TApplication( "App", 0, 0 );
  }
  return *evDisp_;
}

void EvDisp::Initialize(void)
{
  //nTrackBdcOut_=0;
  //nTrackSdcIn_=0;
  nTrackSdcOut_=0;
  //nVertex_=0;
  
  gevdisp_ = new TGeometry("evdisp","K1.8 Event Display");

  ThreeVector worldSize(1000.0, 1000.0, 1000.0); /*mm*/
  world_ = new TBRIK("world_","world","void",
		      worldSize.x(), worldSize.y(), worldSize.z());

  node_ = new TNode("node_","node","world_",0.0,0.0,0.0);
  gevdisp_->GetNode("node_")->SetVisibility(0);

  ConstructBH2();
  std::cout << "Finish Construction of BH2" << std::endl;
  //ConstructTarget();
  //std::cout << "Finish Construction of Target" << std::endl;
  ConstructSKS();
  std::cout << "Finish Construction of SKS" << std::endl;
  /*
  ConstructBdcOut();
  std::cout << "Finish Construction of BdcOut" << std::endl;
  */
  ConstructSdcIn();
  std::cout << "Finish Construction of SdcIn" << std::endl;
  ConstructSdcOut();
  std::cout << "Finish Construction of SdcOut" << std::endl;

  ConstructTOF();
  std::cout << "Finish Construction of Tof" << std::endl;


  tc_ = new TCanvas("canvas","E559 Event Display",700,700);
  tp_ = new TPad("pad","E559 Event",0.0,0.0,1.0,1.0,10);
  tp_->Draw();

  tp_->cd();
  gevdisp_->Draw();
  //tp_->GetView()->SetParallel();
  tc_->cd();

  tc_->Update();

  gevdisp_vtx_ = new TGeometry("evdisp_vertex","E559 Event Display");

  world_vtx_ = new TBRIK("world_vtx_","world","void",
		      worldSize.x(), worldSize.y(), worldSize.z());

  node_vtx_ = new TNode("node_vtx_","node_vtx","world_vtx_",0.0,0.0,0.0);
  gevdisp_vtx_->GetNode("node_vtx_")->SetVisibility(0);


  ConstructBH2Vtx();
  std::cout << "Finish Construction of BH2Vtx" << std::endl;

  /*
  ConstructTargetVtx();
  std::cout << "Finish Construction of Target at vertex display" << std::endl;
  ConstructBdcOutVtx();
  std::cout << "Finish Construction of BdcOut at vertex display" << std::endl;
  */
  ConstructSdcInVtx();
  std::cout << "Finish Construction of SdcIn at vertex display" << std::endl;

  tc_vtx_ = new TCanvas("canvas_vtx","E559 Event Display",700,700);
  tp_vtx_ = new TPad("pad_vtx","E559 Event",0.0,0.0,1.0,1.0,10);
  tp_vtx_->Draw();

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tp_vtx_->GetView()->ZoomIn();
  //tp_vtx_->GetView()->SetParallel();
  tc_vtx_->cd();

  tc_vtx_->Update();

  ResetVisibility();
}
#if 0
void EvDisp::ConstructTarget(void)
{
  static const std::string funcname = "EvDisp::ConstructSKS";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //-----Target
  double TargetRmin = 0.0; 
  double TargetRmax = 33.9; 
  double TargetZ    = 110.0/2.0; 

  Target_Tube_ = new TTUBE("Target_Tube_","Target_Tube_", "void", 
			   TargetRmin, TargetRmax, TargetZ);

  double rotMatTarget[9];
  calcRotMatrix(90.0, 90.0, -50.0, rotMatTarget);
  TRotMatrix *rotTarget = new TRotMatrix("rotTarget","rotTarget",rotMatTarget);

  int lnum = geomMan.GetDetectorId("Target");
  //ThreeVector GlobalPos = geomMan.GetGlobalPosition(lnum);
  ThreeVector GlobalPos = geomMan.Local2GlobalPos(lnum, 
	  ThreeVector(0.0, 0.0, geomMan.GetLocalZ(lnum)+30.0));

  Target_Node_ = new TNode("Target_Node_", "Target_Node", "Target_Tube_",
				   GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
				   "rotTarget", "void");
}

void EvDisp::ConstructTargetVtx(void)
{
  static const std::string funcname = "EvDisp::ConstructSKS";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  //-----Target
  double TargetRmin = 0.0; 
  double TargetRmax = 33.9; 
  double TargetZ    = 110.0/2.0; 

  TTUBE *Target_Tube_vtx_ = new TTUBE("Target_Tube_vtx_","Target_Tube_vtx_", "void", 
			   TargetRmin, TargetRmax, TargetZ);

  //-----Target
  double rotMatTarget[9];
  calcRotMatrix(90.0, 90.0, -50.0, rotMatTarget);
  TRotMatrix *rotTarget = new TRotMatrix("rotTarget","rotTarget",rotMatTarget);

  int lnum = geomMan.GetDetectorId("Target");
  //ThreeVector GlobalPos = geomMan.GetGlobalPosition(lnum);
  ThreeVector GlobalPos = geomMan.Local2GlobalPos(lnum, 
	  ThreeVector(0.0, 0.0, geomMan.GetLocalZ(lnum)+30.0));

  Target_Node_vtx_ = new TNode("Target_Node_vtx_", "Target_Node_vtx", "Target_Tube_vtx_",
				   GlobalPos.x(), GlobalPos.y(), GlobalPos.z(),
				   "rotTarget", "void");
}
#endif

void EvDisp::ConstructBH2(void)
{
  static const std::string funcname = "EvDisp::ConstructBH2";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BH2
  double rotMatBh2[9];
  double Bh2WallX = 130.0/2.0;
  double Bh2WallY = 14.0/2.0;
  double Bh2WallZ = 60.0/2.0;

  double Bh2SegX[7] = {30./2., 20./2., 15./2., 15./2., 15./2., 20./2., 30./2.};
  double Bh2SegY[7] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};
  double Bh2SegZ[7] = {60./2., 60./2., 60./2., 60./2., 60./2., 60./2., 60./2.};

  double localPosX[7] = {-51.5, -28.5, -13., 0., 13., 28.5, 51.5};
  double localPosZ[7] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};

  lnum = geomMan.GetDetectorId("BH2(global)");
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBh2);

  TRotMatrix *rotBh2 = new TRotMatrix("rotBh2","rotBh2",rotMatBh2);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector bh2WallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector bh2WallGlobalPos = geomMan.Local2GlobalPos(lnum, bh2WallLocalPos);
  sprintf(object_name, "Bh2Wall_Brik_");
  Bh2Wall_Brik_ = new TBRIK(object_name, object_name, "void", 
			    Bh2WallX, Bh2WallY, Bh2WallZ);    
  sprintf(node_name, "Bh2Wall_Node_");
  Bh2Wall_Node_ = new TNode(node_name, node_name, object_name,
			    bh2WallGlobalPos.x(),
			    bh2WallGlobalPos.y(), 
			    bh2WallGlobalPos.z(),
			    "rotBh2", "void");
  Bh2Wall_Node_->SetVisibility(0);

  Bh2Wall_Node_->cd();

  for (int i=0; i<NumOfSegBH2; i++) {
    sprintf(object_name, "Bh2Seg_Brik_%d", i);
    Bh2Seg_Brik_[i] = new TBRIK(object_name, object_name, "void", 
				Bh2SegX[i], Bh2SegY[i], Bh2SegZ[i]);    

    ThreeVector bh2SegLocalPos = 
      ThreeVector(localPosX[i], 0.0, localPosZ[i]);
    sprintf(node_name, "Bh2Seg_Node_%d", i);
    Bh2Seg_Node_[i] = new TNode(node_name, node_name, object_name,
		bh2SegLocalPos.x(), bh2SegLocalPos.y(), bh2SegLocalPos.z());
  }
  node_->cd();

}

void EvDisp::ConstructBH2Vtx(void)
{
  static const std::string funcname = "EvDisp::ConstructBH2Vtx";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BH2
  double rotMatBh2[9];
  double Bh2WallX = 130.0/2.0;
  double Bh2WallY = 14.0/2.0;
  double Bh2WallZ = 60.0/2.0;

  double Bh2SegX[7] = {30./2., 20./2., 15./2., 15./2., 15./2., 20./2., 30./2.};
  double Bh2SegY[7] = {5./2., 5./2., 5./2., 5./2., 5./2., 5./2., 5./2.};
  double Bh2SegZ[7] = {60./2., 60./2., 60./2., 60./2., 60./2., 60./2., 60./2.};

  double localPosX[7] = {-51.5, -28.5, -13., 0., 13., 28.5, 51.5};
  double localPosZ[7] = {4.5, -4.5, 4.5, -4.5, 4.5, -4.5, 4.5};

  lnum = geomMan.GetDetectorId("BH2(global)");
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBh2);

  TRotMatrix *rotBh2 = new TRotMatrix("rotBh2","rotBh2",rotMatBh2);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector bh2WallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector bh2WallGlobalPos = geomMan.Local2GlobalPos(lnum, bh2WallLocalPos);
  sprintf(object_name, "Bh2Wall_Brik_");
  Bh2Wall_Brik_vtx_ = new TBRIK(object_name, object_name, "void", 
			    Bh2WallX, Bh2WallY, Bh2WallZ);    
  sprintf(node_name, "Bh2Wall_Node_");
  Bh2Wall_Node_vtx_ = new TNode(node_name, node_name, object_name,
			    bh2WallGlobalPos.x(),
			    bh2WallGlobalPos.y(), 
			    bh2WallGlobalPos.z(),
			    "rotBh2", "void");
  Bh2Wall_Node_vtx_->SetVisibility(0);

  Bh2Wall_Node_vtx_->cd();

  for (int i=0; i<NumOfSegBH2; i++) {
    sprintf(object_name, "Bh2Seg_Brik_vtx_%d", i);
    Bh2Seg_Brik_vtx_[i] = new TBRIK(object_name, object_name, "void", 
				Bh2SegX[i], Bh2SegY[i], Bh2SegZ[i]);    

    ThreeVector bh2SegLocalPos = 
      ThreeVector(localPosX[i], 0.0, localPosZ[i]);
    sprintf(node_name, "Bh2Seg_Node_vtx_%d", i);
    Bh2Seg_Node_vtx_[i] = new TNode(node_name, node_name, object_name,
		bh2SegLocalPos.x(), bh2SegLocalPos.y(), bh2SegLocalPos.z());
  }
  node_vtx_->cd();

}

void EvDisp::ConstructSKS(void)
{

  static const std::string funcname = "EvDisp::ConstructSKS";  

  //-----RotationMatrix
  TRotMatrix *rot = new TRotMatrix("rot","rot",90.0,0.0,0.0,0.0,90.0,-90.0);

  double Yoffset=1024.52;
  //-----Yoke1
  ThreeVector Yoke1Pos(0.0, -1000.0+Yoffset, 0.0); /*mm*/
  double Yoke1Rmin = 2005.0; 
  double Yoke1Rmax = 3405.0; 
  double Yoke1Z    = 500.0/2.0; 
  double Yoke1Phi1 = 44.488;  /*deg*/
  double Yoke1Phi2 = 91.024;  /*deg*/
  
  Yoke1_Tubs_ = new TTUBS("Yoke1_Tubs_", "SKS Yoke1", "void",
			 Yoke1Rmin, Yoke1Rmax, Yoke1Z,
			 Yoke1Phi1, Yoke1Phi1+Yoke1Phi2);

  node_->cd();
  Yoke1_Node_ = new TNode("Yoke1_Node_", "Yoke1_Node", "Yoke1_Tubs_",
			  Yoke1Pos.x(), Yoke1Pos.y(), Yoke1Pos.z());

  //-----Yoke2
  ThreeVector Yoke2Pos(0.0, -1058.0+Yoffset, 0.0); /*mm*/
  double Yoke2dX1 = 1180.0/2.0;
  double Yoke2dX2 = 1655.0/2.0;
  double Yoke2dY  = 500.0/2.0; 
  double Yoke2dZ  = 390.0/2.0;
  
  Yoke2_Trd_ = new TTRD1("Yoke2_Trd_", "SKS Yoke2", "void",
			 Yoke2dX1, Yoke2dX2, Yoke2dY, Yoke2dZ);

  node_->cd();
  Yoke2_Node_ = new TNode("Yoke2_Node_", "Yoke2_Node", "Yoke2_Trd_",
			  Yoke2Pos.x(), Yoke2Pos.y(), Yoke2Pos.z(),
			  "rot", "void");

  //-----Yoke3
  ThreeVector Yoke3Pos(0.0, -1855.5+Yoffset, 0.0); /*mm*/
  double Yoke3X = 1515.0/2.0;
  double Yoke3Y = 1205.0/2.0;
  double Yoke3Z  = 500.0/2.0; 

  Yoke3_Brik_ = new TBRIK("Yoke3_Brik_", "SKS Yoke3", "void",
			 Yoke3X, Yoke3Y, Yoke3Z);

  node_->cd();
  Yoke3_Node_ = new TNode("Yoke3_Node_", "Yoke3_Node", "Yoke3_Brik_",
			  Yoke3Pos.x(), Yoke3Pos.y(), Yoke3Pos.z());

  //-----Yoke4
  ThreeVector Yoke4Pos(0.0, -8.0+Yoffset, -504.0); /*mm*/
  double Yoke4dX1 = 3296.64/2.0;
  double Yoke4dX2 = 1656.0/2.0;
  double Yoke4dY  = 500.0/2.0; 
  double Yoke4dZ  = 2258.0/2.0;
  
  Yoke4_Trd_ = new TTRD1("Yoke4_Trd_", "SKS Yoke4", "void",
			 Yoke4dX1, Yoke4dX2, Yoke4dY, Yoke4dZ);

  node_->cd();
  Yoke4_Node_ = new TNode("Yoke4_Node_", "Yoke4_Node", "Yoke4_Trd_",
			  Yoke4Pos.x(), Yoke4Pos.y(), Yoke4Pos.z(),
			  "rot", "void");

  //-----Yoke5
  ThreeVector Yoke5Pos(0.0, -8.0+Yoffset, 504.0); /*mm*/
  double Yoke5dX1 = 3296.64/2.0;
  double Yoke5dX2 = 1656.0/2.0;
  double Yoke5dY  = 500.0/2.0; 
  double Yoke5dZ  = 2258.0/2.0;
  
  Yoke5_Trd_ = new TTRD1("Yoke5_Trd_", "SKS Yoke5", "void",
			 Yoke5dX1, Yoke5dX2, Yoke5dY, Yoke5dZ);

  node_->cd();
  Yoke5_Node_ = new TNode("Yoke5_Node_", "Yoke5_Node", "Yoke5_Trd_",
			  Yoke5Pos.x(), Yoke5Pos.y(), Yoke5Pos.z(),
			  "rot", "void");

  return;

}
#if 0
void EvDisp::ConstructBdcOut(void)
{
  static const std::string funcname = "EvDisp::ConstructBdcOut";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  const int OffsetGlobal2Local = 30;

  int lnum;
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BDC3X
  double Bdc3xRmin = 0.0; 
  double Bdc3xRmax = 0.01; 
  double Bdc3xZ    = 400.0/2.0; 
  //-- x
  lnum = geomMan.GetDetectorId("BDC3-x-1(global)");
  sprintf(object_name, "Bdc3x_Tube");
  TTUBE *Bdc3x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3x_Node_%d", wire);
    Bdc3x_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- xp
  lnum = geomMan.GetDetectorId("BDC3-x-2(global)");
  sprintf(object_name, "Bdc3xp_Tube");
  TTUBE *Bdc3xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3xp_Node_%d", wire);
    Bdc3xp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-----BDC3V
  double Bdc3vRmin = 0.0; 
  double Bdc3vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  double Bdc3vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3v);
  TRotMatrix *rotBdc3v = new TRotMatrix("rotBdc3v","rotBdc3v",rotMatBdc3v);

  //-- vp
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  sprintf(object_name, "Bdc3vp_Tube");
  TTUBE *Bdc3vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3vp_Node_%d", wire);
    Bdc3vp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3v", "void");
  }

  //-- v
  lnum = geomMan.GetDetectorId("BDC3-v-2(global)");
  sprintf(object_name, "Bdc3v_Tube");
  TTUBE *Bdc3v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3v_Node_%d", wire);
    Bdc3v_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3v", "void");
  }

  //-----BDC3U
  double Bdc3uRmin = 0.0; 
  double Bdc3uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  double Bdc3uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3u);
  TRotMatrix *rotBdc3u = new TRotMatrix("rotBdc3u","rotBdc3u",rotMatBdc3u);

  //-- up
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  sprintf(object_name, "Bdc3up_Tube");
  TTUBE *Bdc3up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3up_Node_%d", wire);
    Bdc3up_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3u", "void");
  }

  //-- u
  lnum = geomMan.GetDetectorId("BDC3-u-2(global)");
  sprintf(object_name, "Bdc3u_Tube");
  TTUBE *Bdc3u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3u_Node_%d", wire);
    Bdc3u_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3u", "void");
  }


  //-----BDC4U
  double Bdc4uRmin = 0.0; 
  double Bdc4uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  double Bdc4uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4u);
  TRotMatrix *rotBdc4u = new TRotMatrix("rotBdc4u","rotBdc4u",rotMatBdc4u);

  //-- u
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  sprintf(object_name, "Bdc4u_Tube");
  TTUBE *Bdc4u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4u_Node_%d", wire);
    Bdc4u_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4u", "void");
  }

  //-- up
  lnum = geomMan.GetDetectorId("BDC4-u-2(global)");
  sprintf(object_name, "Bdc4up_Tube");
  TTUBE *Bdc4up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4up_Node_%d", wire);
    Bdc4up_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4u", "void");
  }

  //-----BDC4V
  double Bdc4vRmin = 0.0; 
  double Bdc4vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  double Bdc4vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4v);
  TRotMatrix *rotBdc4v = new TRotMatrix("rotBdc4v","rotBdc4v",rotMatBdc4v);

  //-- v
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  sprintf(object_name, "Bdc4v_Tube");
  TTUBE *Bdc4v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4v_Node_%d", wire);
    Bdc4v_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4v", "void");
  }

  //-- vp
  lnum = geomMan.GetDetectorId("BDC4-v-2(global)");
  sprintf(object_name, "Bdc4vp_Tube");
  TTUBE *Bdc4vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4vp_Node_%d", wire);
    Bdc4vp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4v", "void");
  }

  //-----BDC4X
  double Bdc4xRmin = 0.0; 
  double Bdc4xRmax = 0.01; 
  double Bdc4xZ    = 400.0/2.0; 

  //-- xp
  lnum = geomMan.GetDetectorId("BDC4-x-1(global)");
  sprintf(object_name, "Bdc4xp_Tube");
  TTUBE *Bdc4xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4xp_Node_%d", wire);
    Bdc4xp_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- x
  lnum = geomMan.GetDetectorId("BDC4-x-2(global)");
  sprintf(object_name, "Bdc4x_Tube");
  TTUBE *Bdc4x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4x_Node_%d", wire);
    Bdc4x_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }
}
#endif

#if 0
void EvDisp::ConstructBdcOutVtx(void)
{
  static const std::string funcname = "EvDisp::ConstructBdcOutVtx";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  const int OffsetGlobal2Local = 30;
  
  int lnum;
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----BDC3X
  double Bdc3xRmin = 0.0; 
  double Bdc3xRmax = 0.01; 
  double Bdc3xZ    = 400.0/2.0; 
  //-- x
  lnum = geomMan.GetDetectorId("BDC3-x-1(global)");
  sprintf(object_name, "Bdc3x_Tube");
  TTUBE *Bdc3x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3x_Node_vtx_%d", wire);
    Bdc3x_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- xp
  lnum = geomMan.GetDetectorId("BDC3-x-2(global)");
  sprintf(object_name, "Bdc3xp_Tube");
  TTUBE *Bdc3xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3xRmin, Bdc3xRmax, Bdc3xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3xp_Node_vtx_%d", wire);
    Bdc3xp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-----BDC3V
  double Bdc3vRmin = 0.0; 
  double Bdc3vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  double Bdc3vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3v);
  TRotMatrix *rotBdc3v = new TRotMatrix("rotBdc3v","rotBdc3v",rotMatBdc3v);

  //-- vp
  lnum = geomMan.GetDetectorId("BDC3-v-1(global)");
  sprintf(object_name, "Bdc3vp_Tube");
  TTUBE *Bdc3vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3vp_Node_vtx_%d", wire);
    Bdc3vp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3v", "void");
  }

  //-- v
  lnum = geomMan.GetDetectorId("BDC3-v-2(global)");
  sprintf(object_name, "Bdc3v_Tube");
  TTUBE *Bdc3v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3vRmin, Bdc3vRmax, Bdc3vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3v_Node_vtx_%d", wire);
    Bdc3v_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3v", "void");
  }

  //-----BDC3U
  double Bdc3uRmin = 0.0; 
  double Bdc3uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  double Bdc3uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc3u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc3u);
  TRotMatrix *rotBdc3u = new TRotMatrix("rotBdc3u","rotBdc3u",rotMatBdc3u);

  //-- up
  lnum = geomMan.GetDetectorId("BDC3-u-1(global)");
  sprintf(object_name, "Bdc3up_Tube");
  TTUBE *Bdc3up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3up_Node_vtx_%d", wire);
    Bdc3up_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc3u", "void");
  }

  //-- u
  lnum = geomMan.GetDetectorId("BDC3-u-2(global)");
  sprintf(object_name, "Bdc3u_Tube");
  TTUBE *Bdc3u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc3uRmin, Bdc3uRmax, Bdc3uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc3u_Node_vtx_%d", wire);
    Bdc3u_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc3u", "void");
  }


  //-----BDC4U
  double Bdc4uRmin = 0.0; 
  double Bdc4uRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  double Bdc4uZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4u[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4u);
  TRotMatrix *rotBdc4u = new TRotMatrix("rotBdc4u","rotBdc4u",rotMatBdc4u);

  //-- u
  lnum = geomMan.GetDetectorId("BDC4-u-1(global)");
  sprintf(object_name, "Bdc4u_Tube");
  TTUBE *Bdc4u_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4u_Node_vtx_%d", wire);
    Bdc4u_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4u", "void");
  }

  //-- up
  lnum = geomMan.GetDetectorId("BDC4-u-2(global)");
  sprintf(object_name, "Bdc4up_Tube");
  TTUBE *Bdc4up_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4uRmin, Bdc4uRmax, Bdc4uZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4up_Node_vtx_%d", wire);
    Bdc4up_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4u", "void");
  }

  //-----BDC4V
  double Bdc4vRmin = 0.0; 
  double Bdc4vRmax = 0.01; 
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  double Bdc4vZ    = 400.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatBdc4v[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatBdc4v);
  TRotMatrix *rotBdc4v = new TRotMatrix("rotBdc4v","rotBdc4v",rotMatBdc4v);

  //-- v
  lnum = geomMan.GetDetectorId("BDC4-v-1(global)");
  sprintf(object_name, "Bdc4v_Tube");
  TTUBE *Bdc4v_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4v_Node_vtx_%d", wire);
    Bdc4v_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
				    "rotBdc4v", "void");
  }

  //-- vp
  lnum = geomMan.GetDetectorId("BDC4-v-2(global)");
  sprintf(object_name, "Bdc4vp_Tube");
  TTUBE *Bdc4vp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4vRmin, Bdc4vRmax, Bdc4vZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4vp_Node_vtx_%d", wire);
    Bdc4vp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
			     wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z(),
			     "rotBdc4v", "void");
  }

  //-----BDC4X
  double Bdc4xRmin = 0.0; 
  double Bdc4xRmax = 0.01; 
  double Bdc4xZ    = 400.0/2.0; 

  //-- xp
  lnum = geomMan.GetDetectorId("BDC4-x-1(global)");
  sprintf(object_name, "Bdc4xp_Tube");
  TTUBE *Bdc4xp_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4xp_Node_vtx_%d", wire);
    Bdc4xp_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }

  //-- x
  lnum = geomMan.GetDetectorId("BDC4-x-2(global)");
  sprintf(object_name, "Bdc4x_Tube");
  TTUBE *Bdc4x_Tube = new TTUBE(object_name, object_name, "void", 
			   Bdc4xRmin, Bdc4xRmax, Bdc4xZ);

  for (wire=1; wire<= MaxWireBDC; wire++) {
    localPos = geomMan.calcWirePosition(lnum-OffsetGlobal2Local, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Bdc4x_Node_vtx_%d", wire);
    Bdc4x_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				    wireGlobalPos.x(), wireGlobalPos.y(), wireGlobalPos.z());
  }
}
#endif

#if 1
void EvDisp::ConstructSdcIn(void)
{
  static const std::string funcname = "EvDisp::ConstructSdcIn";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  
  int lnum; /* layer number */
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----SDC1
  //-- u
  lnum = 1;
  double Sdc1uRmin = 0.0; 
  double Sdc1uRmax = 0.01; 
  double Sdc1uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u1);
  TRotMatrix *rotSdc1u1 = new TRotMatrix("rotSdc1u1","rotSdc1u1",rotMatSdc1u1);
  
  sprintf(object_name, "Sdc1u_Tube");
  TTUBE *Sdc1u_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1uRmin, Sdc1uRmax, Sdc1uZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u1_Node_%d", wire);
    Sdc1u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				     wireGlobalPos.x(), 
				     wireGlobalPos.y(), 
				     wireGlobalPos.z(),
				     "rotSdc1u1", "void");
  }
  //-- u2
  lnum = 2;
  double rotMatSdc1u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u2);
  TRotMatrix *rotSdc1u2 = new TRotMatrix("rotSdc1u2","rotSdc1u2",rotMatSdc1u2);
  
  for (wire=1; wire <= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u2_Node_%d", wire);
    Sdc1u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc1u2", "void");
  }

  //-- v1
  lnum = 3;
  double Sdc1vRmin = 0.0; 
  double Sdc1vRmax = 0.01; 
  double Sdc1vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v1);
  TRotMatrix *rotSdc1v1 = new TRotMatrix("rotSdc1v1","rotSdc1v1",rotMatSdc1v1);
  
  sprintf(object_name, "Sdc1v_Tube");
  TTUBE *Sdc1v_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1vRmin, Sdc1vRmax, Sdc1vZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v1_Node_%d", wire);
    Sdc1v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				     wireGlobalPos.x(), 
				     wireGlobalPos.y(), 
				     wireGlobalPos.z(),
				     "rotSdc1v1", "void");
  }
  //-- v2
  lnum = 4;
  double rotMatSdc1v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v2);
  TRotMatrix *rotSdc1v2 = new TRotMatrix("rotSdc1v2","rotSdc1v2",rotMatSdc1v2);
  
  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v2_Node_%d", wire);
    Sdc1v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc1v2", "void");
  }

  //-----SDC2
  //-- v
  lnum = 5;
  double Sdc2vRmin = 0.0; 
  double Sdc2vRmax = 0.01; 
  double Sdc2vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v1);
  TRotMatrix *rotSdc2v1 = new TRotMatrix("rotSdc2v1","rotSdc2v1",rotMatSdc2v1);

  sprintf(object_name, "Sdc2v1_Tube");
  TTUBE *Sdc2v_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2vRmin, Sdc2vRmax, Sdc2vZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v1_Node_%d", wire);
    Sdc2v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2v1", "void");
  }
  // v2
  lnum = 6;
  double rotMatSdc2v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v2);
  TRotMatrix *rotSdc2v2 = new TRotMatrix("rotSdc2v2","rotSdc2v2",rotMatSdc2v2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v2_Node_%d", wire);
    Sdc2v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2v2", "void");
  }

  //-- u1
  lnum = 7;
  double Sdc2uRmin = 0.0; 
  double Sdc2uRmax = 0.01; 
  double Sdc2uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u1);
  TRotMatrix *rotSdc2u1 = new TRotMatrix("rotSdc2u1","rotSdc2u1",rotMatSdc2u1);

  sprintf(object_name, "Sdc2u1_Tube");
  TTUBE *Sdc2u_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2uRmin, Sdc2uRmax, Sdc2uZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u1_Node_%d", wire);
    Sdc2u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2u1", "void");
  }
  // u2
  lnum = 8;
  double rotMatSdc2u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u2);
  TRotMatrix *rotSdc2u2 = new TRotMatrix("rotSdc2u2","rotSdc2u2",rotMatSdc2u2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u2_Node_%d", wire);
    Sdc2u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z(),
				   "rotSdc2u2", "void");
  }

  double Sdc2xRmin = 0.0; 
  double Sdc2xRmax = 0.01; 
  double Sdc2xZ    = 200.0/2.0; 
  //-- x1
  lnum = 9;

  sprintf(object_name, "Sdc2x_Tube");
  TTUBE *Sdc2x_Tube = new TTUBE(object_name, object_name, "void", 
			  Sdc2xRmin, Sdc2xRmax, Sdc2xZ);
  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x1_Node_%d", wire);
    Sdc2x1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(),
				   wireGlobalPos.z());
  }

  //-- x2
  lnum = 10;

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x2_Node_%d", wire);
    Sdc2x2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

}
#endif 



void EvDisp::ConstructSdcInVtx(void)
{
  static const std::string funcname = "EvDisp::ConstructSdcInVtx";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  
  int lnum; /* layer number */
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----SDC1
  //-- u
  lnum = 1;
  double Sdc1uRmin = 0.0; 
  double Sdc1uRmax = 0.01; 
  double Sdc1uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u1);
  TRotMatrix *rotSdc1u1 = new TRotMatrix("rotSdc1u1","rotSdc1u1",rotMatSdc1u1);
  
  sprintf(object_name, "Sdc1u_Tube");
  TTUBE *Sdc1u_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1uRmin, Sdc1uRmax, Sdc1uZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u1_Node_vtx_%d", wire);
    Sdc1u1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1u1", "void");
  }
  //-- u2
  lnum = 2;
  double rotMatSdc1u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1u2);
  TRotMatrix *rotSdc1u2 = new TRotMatrix("rotSdc1u2","rotSdc1u2",rotMatSdc1u2);
  
  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1u2_Node_vtx_%d", wire);
    Sdc1u2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1u2", "void");
  }

  //-- v1
  lnum = 3;
  double Sdc1vRmin = 0.0; 
  double Sdc1vRmax = 0.01; 
  double Sdc1vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc1v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v1);
  TRotMatrix *rotSdc1v1 = new TRotMatrix("rotSdc1v1","rotSdc1v1",rotMatSdc1v1);
  
  sprintf(object_name, "Sdc1v_Tube");
  TTUBE *Sdc1v_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc1vRmin, Sdc1vRmax, Sdc1vZ);

  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v1_Node_vtx_%d", wire);
    Sdc1v1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1v1", "void");
  }
  //-- v2
  lnum = 4;
  double rotMatSdc1v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
 		geomMan.GetRotAngle2(lnum), rotMatSdc1v2);
  TRotMatrix *rotSdc1v2 = new TRotMatrix("rotSdc1v2","rotSdc1v2",rotMatSdc1v2);
  
  for (wire=1; wire<= MaxWireSDC1; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc1v2_Node_vtx_%d", wire);
    Sdc1v2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z(),
				       "rotSdc1v2", "void");
  }

  //-----SDC2
  //-- v
  lnum = 5;
  double Sdc2vRmin = 0.0; 
  double Sdc2vRmax = 0.01; 
  double Sdc2vZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2v1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v1);
  TRotMatrix *rotSdc2v1 = new TRotMatrix("rotSdc2v1","rotSdc2v1",rotMatSdc2v1);

  sprintf(object_name, "Sdc2v1_Tube");
  TTUBE *Sdc2v_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2vRmin, Sdc2vRmax, Sdc2vZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v1_Node_vtx_%d", wire);
    Sdc2v1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2v1", "void");
  }
  // v2
  lnum = 6;
  double rotMatSdc2v2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2v2);
  TRotMatrix *rotSdc2v2 = new TRotMatrix("rotSdc2v2","rotSdc2v2",rotMatSdc2v2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2v2_Node_vtx_%d", wire);
    Sdc2v2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2v2", "void");
  }

  //-- u1
  lnum = 7;
  double Sdc2uRmin = 0.0; 
  double Sdc2uRmax = 0.01; 
  double Sdc2uZ    = 200.0/cos(geomMan.GetTiltAngle(lnum)*Deg2Rad)/2.0; 
  double rotMatSdc2u1[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u1);
  TRotMatrix *rotSdc2u1 = new TRotMatrix("rotSdc2u1","rotSdc2u1",rotMatSdc2u1);

  sprintf(object_name, "Sdc2u1_Tube");
  TTUBE *Sdc2u_Tube = new TTUBE(object_name, object_name, "void", 
				Sdc2uRmin, Sdc2uRmax, Sdc2uZ);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u1_Node_vtx_%d", wire);
    Sdc2u1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2u1", "void");
  }
  // u2
  lnum = 8;
  double rotMatSdc2u2[9];
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc2u2);
  TRotMatrix *rotSdc2u2 = new TRotMatrix("rotSdc2u2","rotSdc2u2",rotMatSdc2u2);

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2u2_Node_vtx_%d", wire);
    Sdc2u2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z(),
				       "rotSdc2u2", "void");
  }

  double Sdc2xRmin = 0.0; 
  double Sdc2xRmax = 0.01; 
  double Sdc2xZ    = 200.0/2.0; 
  //-- x1
  lnum = 9;

  sprintf(object_name, "Sdc2x_Tube");
  TTUBE *Sdc2x_Tube = new TTUBE(object_name, object_name, "void", 
			  Sdc2xRmin, Sdc2xRmax, Sdc2xZ);
  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x1_Node_vtx_%d", wire);
    Sdc2x1_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(),
				       wireGlobalPos.y(),
				       wireGlobalPos.z());
  }

  //-- x2
  lnum = 10;

  for (wire=1; wire<= MaxWireSDC2; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc2x2_Node_vtx_%d", wire);
    Sdc2x2_Node_vtx_[wire-1] = new TNode(node_name, node_name, object_name,
				       wireGlobalPos.x(), 
				       wireGlobalPos.y(), 
				       wireGlobalPos.z());
  }

}

#if 1
void EvDisp::ConstructSdcOut(void)
{
  static const std::string funcname = "EvDisp::ConstructSdcOut";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  
  int lnum; /* layer number */
  int wire;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----SDC3
  //-- v1
  double Sdc3vRmin = 0.0; 
  double Sdc3vRmax = 0.01; 
  double Sdc3vZ    = 1240.0/2.0; 
  double rotMatSdc3v1[9];
  lnum = 31;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3v1);
  TRotMatrix *rotSdc3v1 = new TRotMatrix("rotSdc3v1","rotSdc3v1",rotMatSdc3v1);

  sprintf(object_name, "Sdc3v1_Tube");
  TTUBE *Sdc3v1_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc3vRmin, Sdc3vRmax, Sdc3vZ);

  for (wire=1; wire<= MaxWireSDC3V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3v1_Node_%d", wire);
    Sdc3v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3v1", "void");
  }

  //-- x1
  double Sdc3xRmin = 0.0; 
  double Sdc3xRmax = 0.01; 
  double Sdc3xZ    = 1240.0/2.0; 
  lnum = 32;
  sprintf(object_name, "Sdc3x1_Tube");
  TTUBE *Sdc3x1_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc3xRmin, Sdc3xRmax, Sdc3xZ);

  for (wire=1; wire<= MaxWireSDC3X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3x1_Node_%d", wire);
    Sdc3x1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u1
  double Sdc3uRmin = 0.0; 
  double Sdc3uRmax = 0.01; 
  double Sdc3uZ    = 1240.0/2.0; 
  double rotMatSdc3u1[9];
  lnum = 33;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3u1);
  TRotMatrix *rotSdc3u1 = new TRotMatrix("rotSdc3u1","rotSdc3u1",rotMatSdc3u1);

  sprintf(object_name, "Sdc3u1_Tube");
  TTUBE *Sdc3u1_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc3uRmin, Sdc3uRmax, Sdc3uZ);

  for (wire=1; wire<= MaxWireSDC3U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3u1_Node_%d", wire);
    Sdc3u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3u1", "void");
  }

  //-- v2
  double rotMatSdc3v2[9];
  lnum = 34;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3v2);
  TRotMatrix *rotSdc3v2 = new TRotMatrix("rotSdc3v2","rotSdc3v2",rotMatSdc3v2);

  sprintf(object_name, "Sdc3v2_Tube");
  TTUBE *Sdc3v2_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc3vRmin, Sdc3vRmax, Sdc3vZ);

  for (wire=1; wire<= MaxWireSDC3V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3v2_Node_%d", wire);
    Sdc3v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3v2", "void");
  }

  //-- x2
  lnum = 35;
  sprintf(object_name, "Sdc3x2_Tube");
  TTUBE *Sdc3x2_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc3xRmin, Sdc3xRmax, Sdc3xZ);

  for (wire=1; wire<= MaxWireSDC3X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3x2_Node_%d", wire);
    Sdc3x2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u2
  double rotMatSdc3u2[9];
  lnum = 36;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc3u2);
  TRotMatrix *rotSdc3u2 = new TRotMatrix("rotSdc3u2","rotSdc3u2",rotMatSdc3u2);
  sprintf(object_name, "Sdc3u2_Tube");
  TTUBE *Sdc3u2_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc3uRmin, Sdc3uRmax, Sdc3uZ);

  for (wire=1; wire<= MaxWireSDC3U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc3u2_Node_%d", wire);
    Sdc3u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc3u2", "void");
  }

  //-----SDC4
  //-- v1
  double Sdc4vRmin = 0.0; 
  double Sdc4vRmax = 0.01; 
  double Sdc4vZ    = 1240.0/2.0; 
  double rotMatSdc4v1[9];
  lnum = 37;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4v1);
  TRotMatrix *rotSdc4v1 = new TRotMatrix("rotSdc4v1","rotSdc4v1",rotMatSdc4v1);

  sprintf(object_name, "Sdc4v1_Tube");
  TTUBE *Sdc4v1_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc4vRmin, Sdc4vRmax, Sdc4vZ);

  for (wire=1; wire<= MaxWireSDC4V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4v1_Node_%d", wire);
    Sdc4v1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4v1", "void");
  }

  //-- x1
  double Sdc4xRmin = 0.0; 
  double Sdc4xRmax = 0.01; 
  double Sdc4xZ    = 1240.0/2.0; 
  lnum = 38;
  sprintf(object_name, "Sdc4x1_Tube");
  TTUBE *Sdc4x1_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc4xRmin, Sdc4xRmax, Sdc4xZ);

  for (wire=1; wire<= MaxWireSDC4X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4x1_Node_%d", wire);
    Sdc4x1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u1
  double Sdc4uRmin = 0.0; 
  double Sdc4uRmax = 0.01; 
  double Sdc4uZ    = 1240.0/2.0; 
  double rotMatSdc4u1[9];
  lnum = 39;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4u1);
  TRotMatrix *rotSdc4u1 = new TRotMatrix("rotSdc4u1","rotSdc4u1",rotMatSdc4u1);

  sprintf(object_name, "Sdc4u1_Tube");
  TTUBE *Sdc4u1_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc4uRmin, Sdc4uRmax, Sdc4uZ);

  for (wire=1; wire<= MaxWireSDC4U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4u1_Node_%d", wire);
    Sdc4u1_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4u1", "void");
  }

  //-- v2
  double rotMatSdc4v2[9];
  lnum = 40;

  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4v2);
  TRotMatrix *rotSdc4v2 = new TRotMatrix("rotSdc4v2","rotSdc4v2",rotMatSdc4v2);

  sprintf(object_name, "Sdc4v2_Tube");
  TTUBE *Sdc4v2_Tube = new TTUBE(object_name, object_name, "void", 
				  Sdc4vRmin, Sdc4vRmax, Sdc4vZ);

  for (wire=1; wire<= MaxWireSDC4V; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4v2_Node_%d", wire);
    Sdc4v2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4v2", "void");
  }

  //-- x2
  lnum = 41;
  sprintf(object_name, "Sdc4x2_Tube");
  TTUBE *Sdc4x2_Tube = new TTUBE(object_name, object_name, "void", 
			   Sdc4xRmin, Sdc4xRmax, Sdc4xZ);

  for (wire=1; wire<= MaxWireSDC4X; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4x2_Node_%d", wire);
    Sdc4x2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(),
				   wireGlobalPos.y(), 
				   wireGlobalPos.z());
  }

  //-- u1
  double rotMatSdc4u2[9];
  lnum = 42;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatSdc4u2);
  TRotMatrix *rotSdc4u2 = new TRotMatrix("rotSdc4u2","rotSdc4u2",rotMatSdc4u2);
  sprintf(object_name, "Sdc4u2_Tube");
  TTUBE *Sdc4u2_Tube = new TTUBE(object_name, object_name, "void", 
				 Sdc4uRmin, Sdc4uRmax, Sdc4uZ);

  for (wire=1; wire<= MaxWireSDC4U; wire++) {
    localPos = geomMan.calcWirePosition(lnum, wire);
    ThreeVector wireLocalPos = ThreeVector(localPos, 0.0, 0.0);
    ThreeVector wireGlobalPos = geomMan.Local2GlobalPos(lnum, wireLocalPos);

    sprintf(node_name, "Sdc4u2_Node_%d", wire);
    Sdc4u2_Node_[wire-1] = new TNode(node_name, node_name, object_name,
				   wireGlobalPos.x(), 
				   wireGlobalPos.y(), 
				   wireGlobalPos.z(),
				   "rotSdc4u2", "void");
  }

}
#endif

#if 1
void EvDisp::ConstructTOF(void)
{
  static const std::string funcname = "EvDisp::ConstructTOF";  

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();

  int lnum;
  double localPos;
  char object_name[MaxChar];
  char node_name[MaxChar];

  //-----TOF
  double rotMatTof[9];
  double TofWallX = 70.0*32.0/2.0;
  double TofWallY = 30.0/2.0;
  double TofWallZ = 1000.0/2.0;

  double TofSegX = 70.0/2.0;
  double TofSegY = 30.0/2.0;
  double TofSegZ = 1000.0/2.0;

  lnum = 51;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatTof);
  TRotMatrix *rotTof = new TRotMatrix("rotTof","rotTof",rotMatTof);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector tofWallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector tofWallGlobalPos = geomMan.Local2GlobalPos(lnum, tofWallLocalPos);
  sprintf(object_name, "TofWall_Brik_");
  TofWall_Brik_ = new TBRIK(object_name, object_name, "void", 
			    TofWallX, TofWallY, TofWallZ);    
  sprintf(node_name, "TofWall_Node_");
  TofWall_Node_ = new TNode(node_name, node_name, object_name,
	    tofWallGlobalPos.x(), tofWallGlobalPos.y(), tofWallGlobalPos.z(),
	    "rotTof", "void");

  TofWall_Node_->cd();

  sprintf(object_name, "TofSeg_Brik_");
  TofSeg_Brik_ = new TBRIK(object_name, object_name, "void", 
			   TofSegX, TofSegY, TofSegZ);    
  for (int i=0; i<NumOfSegTOF; i++) {
    ThreeVector tofSegLocalPos = 
      ThreeVector((double)(i-14.5)*TofSegX*2.0, 0.0, 0.0);
    sprintf(node_name, "TofSeg_Node_%d", i);
    TofSeg_Node_[i] = new TNode(node_name, node_name, object_name,
		tofSegLocalPos.x(), tofSegLocalPos.y(), tofSegLocalPos.z());
  }
  node_->cd();
#if 0
  //-----AC1
  double rotMatAc1[9];
  double Ac1X = 1050.0/2.0;
  double Ac1Y = 400.0/2.0;
  double Ac1Z = 1200.0/2.0;

  lnum = 52;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatAc1);
  TRotMatrix *rotAc1 = new TRotMatrix("rotAc1","rotAc1",rotMatAc1);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector Ac1LocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector Ac1GlobalPos = geomMan.Local2GlobalPos(lnum, Ac1LocalPos);

  sprintf(object_name, "Ac1_Brik_");
  Ac1_Brik_ = new TBRIK(object_name, object_name, "void", Ac1X, Ac1Y, Ac1Z);
  sprintf(node_name, "Ac1_Node_");
  Ac1_Node_ = new TNode(node_name, node_name, object_name,
	    Ac1GlobalPos.x(), Ac1GlobalPos.y(), Ac1GlobalPos.z(),
	    "rotAc1", "void");

  //-----AC2
  double rotMatAc2[9];
  double Ac2X = 1400.0/2.0;
  double Ac2Y = 400.0/2.0;
  double Ac2Z = 1400.0/2.0;

  lnum = 53;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatAc2);
  TRotMatrix *rotAc2 = new TRotMatrix("rotAc2","rotAc2",rotMatAc2);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector Ac2LocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector Ac2GlobalPos = geomMan.Local2GlobalPos(lnum, Ac2LocalPos);

  sprintf(object_name, "Ac2_Brik_");
  Ac2_Brik_ = new TBRIK(object_name, object_name, "void", Ac2X, Ac2Y, Ac2Z);
  sprintf(node_name, "Ac2_Node_");
  Ac2_Node_ = new TNode(node_name, node_name, object_name,
	    Ac2GlobalPos.x(), Ac2GlobalPos.y(), Ac2GlobalPos.z(),
	    "rotAc2", "void");

#endif
  //-----LC
  double rotMatLc[9];
  double LcWallX = 100.0*28.0/2.0;
  double LcWallY = 40.0/2.0;
  double LcWallZ = 1400.0/2.0;

  double LcSegX = 100.0/2.0;
  double LcSegY = 40.0/2.0;
  double LcSegZ = 1400.0/2.0;

  lnum = 54;
  calcRotMatrix(geomMan.GetTiltAngle(lnum), geomMan.GetRotAngle1(lnum), 
		geomMan.GetRotAngle2(lnum), rotMatLc);
  TRotMatrix *rotLc = new TRotMatrix("rotLc","rotLc",rotMatLc);
  localPos = geomMan.calcWirePosition(lnum, 0);
  ThreeVector lcWallLocalPos = ThreeVector(localPos, 0.0, 0.0);
  ThreeVector lcWallGlobalPos = geomMan.Local2GlobalPos(lnum, lcWallLocalPos);

  sprintf(object_name, "LcWall_Brik_");
  LcWall_Brik_ = new TBRIK(object_name, object_name, "void", 
			    LcWallX, LcWallY, LcWallZ);    
  sprintf(node_name, "LcWall_Node_");
  LcWall_Node_ = new TNode(node_name, node_name, object_name,
	    lcWallGlobalPos.x(), lcWallGlobalPos.y(), lcWallGlobalPos.z(),
	    "rotLc", "void");

  LcWall_Node_->cd();

  sprintf(object_name, "LcSeg_Brik_");
  LcSeg_Brik_ = new TBRIK(object_name, object_name, "void", 
			  LcSegX, LcSegY, LcSegZ);    
  for (int i=0; i<NumOfSegLC; i++) {
    ThreeVector lcSegLocalPos = 
      ThreeVector((double)(i-13.5)*LcSegX*2.0, 0.0, 0.0);

    sprintf(node_name, "LcSeg_Node_%d", i);
    LcSeg_Node_[i] = new TNode(node_name, node_name, object_name,
		lcSegLocalPos.x(), lcSegLocalPos.y(), lcSegLocalPos.z());
  }
  node_->cd();

}
#endif

#if 0
void EvDisp::DrawInitTrack(int nStep, ThreeVector *StepPoint) const
{
  if ( InitStepMark_ )
    delete InitStepMark_;

  InitStepMark_ = new TPolyMarker3D(nStep);
  for (int i=0; i<nStep; i++) {
    InitStepMark_->SetPoint(i, StepPoint[i].x(), StepPoint[i].y(),  StepPoint[i].z());
  }
  InitStepMark_->SetMarkerSize(1);
  InitStepMark_->SetMarkerColor(kCyan);
  InitStepMark_->SetMarkerStyle(6);

  tp_->cd();
  InitStepMark_->Draw();

  tc_->Update();
}

void EvDisp::DrawInitTrack(void) const
{
  tp_->cd();
  if (InitStepMark_)
    InitStepMark_->Draw();

  tc_->Update();

}
#endif

void EvDisp::DrawHitWire(int lnum, int hit_wire, bool range_check, bool tdc_check) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

    switch (lnum) {
    // SDC1,2 1-10
    case 1:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1u1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1u1_Node_vtx_%d", hit_wire);
      break;
    case 2:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1u2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1u2_Node_vtx_%d", hit_wire);
      break;
    case 3:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1v1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1v1_Node_vtx_%d", hit_wire);
      break;
    case 4:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1)
	return;
      sprintf(node_name, "Sdc1v2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1v2_Node_vtx_%d", hit_wire);
      break;
    case 5:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2v1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2v1_Node_vtx_%d", hit_wire);
      break;
    case 6:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2v2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2v2_Node_vtx_%d", hit_wire);
      break;
    case 7:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2u1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2u1_Node_vtx_%d", hit_wire);
      break;
    case 8:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2u2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2u2_Node_vtx_%d", hit_wire);
      break;
    case 9:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2x1_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2x1_Node_vtx_%d", hit_wire);
      break;
    case 10:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2x2_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2x2_Node_vtx_%d", hit_wire);
      break;
    // SDC3,4 31-42
    case 31:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3V)
	return;
      sprintf(node_name, "Sdc3v1_Node_%d", hit_wire);
      break;
    case 32:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3X)
	return;
      sprintf(node_name, "Sdc3x1_Node_%d", hit_wire);
      break;
    case 33:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3U)
	return;
      sprintf(node_name, "Sdc3u1_Node_%d", hit_wire);
      break;
    case 34:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3V)
	return;
      sprintf(node_name, "Sdc3v2_Node_%d", hit_wire);
      break;
    case 35:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3X)
	return;
      sprintf(node_name, "Sdc3x2_Node_%d", hit_wire);
      break;
    case 36:
      if (hit_wire<=0 || hit_wire>MaxWireSDC3U)
	return;
      sprintf(node_name, "Sdc3u2_Node_%d", hit_wire);
      break;
    case 37:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4V)
	return;
      sprintf(node_name, "Sdc4v1_Node_%d", hit_wire);
      break;
    case 38:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4X)
	return;
      sprintf(node_name, "Sdc4x1_Node_%d", hit_wire);
      break;
    case 39:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4U)
	return;
      sprintf(node_name, "Sdc4u1_Node_%d", hit_wire);
      break;
    case 40:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4V)
	return;
      sprintf(node_name, "Sdc4v2_Node_%d", hit_wire);
      break;
    case 41:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4X)
	return;
      sprintf(node_name, "Sdc4x2_Node_%d", hit_wire);
      break;
    case 42:
      if (hit_wire<=0 || hit_wire>MaxWireSDC4U)
	return;
      sprintf(node_name, "Sdc4u2_Node_%d", hit_wire);
      break;
      /*
    // BDC3,4 113-124
    case 113:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3x_Node_vtx_%d", hit_wire);
      break;
    case 114:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3xp_Node_vtx_%d", hit_wire);
      break;
    case 115:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3vp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3vp_Node_vtx_%d", hit_wire);
      break;
    case 116:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3v_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3v_Node_vtx_%d", hit_wire);
      break;
    case 117:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3up_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3up_Node_vtx_%d", hit_wire);
      break;
    case 118:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc3u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc3u_Node_vtx_%d", hit_wire);
      break;
    case 119:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4u_Node_vtx_%d", hit_wire);
      break;
    case 120:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4up_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4up_Node_vtx_%d", hit_wire);
      break;
    case 121:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4v_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4v_Node_vtx_%d", hit_wire);
      break;
    case 122:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4vp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4vp_Node_vtx_%d", hit_wire);
      break;
    case 123:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4xp_Node_vtx_%d", hit_wire);
      break;
    case 124:
      if (hit_wire<=0 || hit_wire>MaxWireBDC)
	return;
      sprintf(node_name, "Bdc4x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Bdc4x_Node_vtx_%d", hit_wire);
      break;
      */
    default:
      std::cerr << "EvDisp::DrawHitWire No such plane ID " << lnum << std::endl;
      return;
    }

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  if (range_check && tdc_check) 
    gevdisp_->GetNode(node_name)->SetLineColor(kBlue);
  else if (range_check && !tdc_check) 
    gevdisp_->GetNode(node_name)->SetLineColor(28);
  else 
    gevdisp_->GetNode(node_name)->SetLineColor(kBlack);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();


  if ((lnum>=1 && lnum<=10) || (lnum>=113 && lnum<=124 )) {
    gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
    if (range_check && tdc_check) 
      gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlue);
    else if (range_check && !tdc_check) 
      gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(28);
    else 
      gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlack);

    tp_vtx_->cd();
    gevdisp_vtx_->Draw();
    tc_vtx_->Update();
  }

}
#if 0
void EvDisp::DrawTrackWire(int lnum, int hit_wire, int it) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

    switch (lnum) {
    // SDC1,2 1-11
    case 1:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1X)
	return;
      sprintf(node_name, "Sdc1x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1x_Node_vtx_%d", hit_wire);
      break;
    case 2:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1X)
	return;
      sprintf(node_name, "Sdc1xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1xp_Node_vtx_%d", hit_wire);
      break;
    case 3:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1Y)
	return;
      sprintf(node_name, "Sdc1y_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1y_Node_vtx_%d", hit_wire);
      break;
    case 4:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1Y)
	return;
      sprintf(node_name, "Sdc1yp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1yp_Node_vtx_%d", hit_wire);
      break;
    case 5:
      if (hit_wire<=0 || hit_wire>MaxWireSDC1U)
	return;
      sprintf(node_name, "Sdc1u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc1u_Node_vtx_%d", hit_wire);
      break;
    case 6:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2vp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2vp_Node_vtx_%d", hit_wire);
      break;
    case 7:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2v_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2v_Node_vtx_%d", hit_wire);
      break;
    case 8:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2up_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2up_Node_vtx_%d", hit_wire);
      break;
    case 9:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2u_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2u_Node_vtx_%d", hit_wire);
      break;
    case 10:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2xp_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2xp_Node_vtx_%d", hit_wire);
      break;
    case 11:
      if (hit_wire<=0 || hit_wire>MaxWireSDC2)
	return;
      sprintf(node_name, "Sdc2x_Node_%d", hit_wire);
      sprintf(node_name_vtx, "Sdc2x_Node_vtx_%d", hit_wire);
      break;
    default:
      std::cerr << "EvDisp::DrawHitWire No such plane ID " << lnum << std::endl;
      return;
    }

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  gevdisp_->GetNode(node_name)->SetLineColor(5+it);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

  gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
  gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(5+it);

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tc_vtx_->Update();

}
#endif

void EvDisp::DrawHitBH2(int seg, int Tu, int Td) const
{
  char node_name[MaxChar];
  char node_name_vtx[MaxChar];

  if (seg<0 || seg>=NumOfSegBH2)
    return;

  sprintf(node_name, "Bh2Seg_Node_%d", seg);
  sprintf(node_name_vtx, "Bh2Seg_Node_vtx_%d", seg);

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  if (Tu>0 && Td>0)
    gevdisp_->GetNode(node_name)->SetLineColor(kBlue);
  else
    gevdisp_->GetNode(node_name)->SetLineColor(kGreen);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

  gevdisp_vtx_->GetNode(node_name_vtx)->SetVisibility(1);
  if (Tu>0 && Td>0)
    gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kBlue);
  else
    gevdisp_vtx_->GetNode(node_name_vtx)->SetLineColor(kGreen);

  tp_vtx_->cd();
  gevdisp_vtx_->Draw();
  tc_vtx_->Update();

}
void EvDisp::DrawHitHodoscope(int lnum, int seg, int Tu, int Td) const
{
  char node_name[MaxChar];

  switch (lnum) {
    // Tof 51
  case 51:
    if (seg<0 || seg>=NumOfSegTOF)
      return;
    sprintf(node_name, "TofSeg_Node_%d", seg);
    break;
    // Tof 54
  case 54:
    if (seg<0 || seg>=NumOfSegLC)
      return;
    sprintf(node_name, "LcSeg_Node_%d", seg);
    break;
  default:
    std::cerr << "EvDisp::DrawHitHodoscope No such plane ID " << lnum << std::endl;
    return;
  }

  gevdisp_->GetNode(node_name)->SetVisibility(1);
  if (Tu>0 && Td>0)
    gevdisp_->GetNode(node_name)->SetLineColor(kBlue);
  else
    gevdisp_->GetNode(node_name)->SetLineColor(kGreen);

  tp_->cd();
  gevdisp_->Draw();
  tc_->Update();

}


#if 0
void EvDisp::DrawBdcOutLocalTrack(ThreeVector globalPos0, 
				  ThreeVector globalPos1) const
{
  if (nTrackBdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBdcOutLocalTrack nTrackBdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackBdcOut_[nTrackBdcOut_] = new TPolyLine3D(2);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineColor(kRed);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineWidth(1);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackBdcOut_[nTrackBdcOut_]->Draw();
  tc_->Update();
  nTrackBdcOut_++;

}

void EvDisp::DrawBdcOutLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackBdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawBdcOutLocalTrack nTrackBdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  int IdBdc3x = DCGeomMan::GetInstance().GetDetectorId("BDC3-x-1");
  double zBdc3x = DCGeomMan::GetInstance().GetLocalZ(IdBdc3x);
  double x0=tp->GetX(zBdc3x), y0=tp->GetY(zBdc3x);

  int IdK18Target = DCGeomMan::GetInstance().GetDetectorId("K18Target");
  double zK18Target = DCGeomMan::GetInstance().GetLocalZ(IdK18Target);
  double x1=tp->GetX(zK18Target), y1=tp->GetY(zK18Target);

  ThreeVector globalPosTgt = DCGeomMan::GetInstance().GetGlobalPosition(0);
  
  ThreeVector pos0 = ThreeVector(zBdc3x-zK18Target, x0, y0);
  pos0.RotateZ(130.0*Deg2Rad);
  ThreeVector globalPos0 = pos0+globalPosTgt;
  
  ThreeVector pos1 = ThreeVector(0.0, x1, y1);
  pos1.RotateZ(130.0*Deg2Rad);
  ThreeVector globalPos1 = pos1+globalPosTgt;
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackBdcOut_[nTrackBdcOut_] = new TPolyLine3D(2);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineColor(kRed);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetLineWidth(1);
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackBdcOut_[nTrackBdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackBdcOut_[nTrackBdcOut_]->Draw();
  tc_->Update();

  tp_vtx_->cd();
  LocalTrackBdcOut_[nTrackBdcOut_]->Draw();
  tc_vtx_->Update();

  nTrackBdcOut_++;

}
#endif


void EvDisp::DrawSdcInLocalTrack(ThreeVector globalPos0, 
				 ThreeVector globalPos1) const
{
  if (nTrackSdcIn_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcIn_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackSdcIn_[nTrackSdcIn_] = new TPolyLine3D(2);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineColor(kRed);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineWidth(1);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcIn_[nTrackSdcIn_]->Draw();
  tc_->Update();
  nTrackSdcIn_++;

}

void EvDisp::DrawSdcInLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackSdcIn_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcIn_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  double x0=tp->GetX0(), y0=tp->GetY0();
  
  int IdSdc2x2 = 10;
  double zSdc2x2 = DCGeomMan::GetInstance().GetLocalZ(IdSdc2x2);
  double x1=tp->GetX(zSdc2x2), y1=tp->GetY(zSdc2x2);
  
  ThreeVector globalPosTgt = DCGeomMan::GetInstance().GetGlobalPosition(0);
  
  ThreeVector pos0 = ThreeVector(0.0, x0, y0);
  pos0.RotateZ(130.0*Deg2Rad);
  ThreeVector globalPos0 = pos0+globalPosTgt;
  
  ThreeVector pos1 = ThreeVector(zSdc2x2, x1, y1);
  pos1.RotateZ(130.0*Deg2Rad);
  ThreeVector globalPos1 = pos1+globalPosTgt;
#if 0  
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif
  LocalTrackSdcIn_[nTrackSdcIn_] = new TPolyLine3D(2);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineColor(kRed);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetLineWidth(1);
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcIn_[nTrackSdcIn_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcIn_[nTrackSdcIn_]->Draw();
  tc_->Update();

  tp_vtx_->cd();
  LocalTrackSdcIn_[nTrackSdcIn_]->Draw();
  tc_vtx_->Update();

  nTrackSdcIn_++;

}


#if 0
void EvDisp::DrawSdcOutLocalTrack(ThreeVector globalPos0, 
				 ThreeVector globalPos1) const
{
  if (nTrackSdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
#if 0 
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif

  LocalTrackSdcOut_[nTrackSdcOut_] = new TPolyLine3D(2);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineColor(kRed);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineWidth(1);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcOut_[nTrackSdcOut_]->Draw();
  tc_->Update();
  nTrackSdcOut_++;
}
#endif

void EvDisp::DrawSdcOutLocalTrack(DCLocalTrack *tp) const
{
  if (nTrackSdcOut_>=MaxTrack) {
    std::cerr << "EvDisp::DrawSdcInLocalTrack nTrackSdcOut_ is greater than MaxTrack"
	      << std::endl;
    return;
  }
  /*
  int IdTof = DCGeomMan::GetInstance().GetTofId();
  double zTof = DCGeomMan::GetInstance().GetLocalZ( IdTof ); 
  double x0=tp->GetX(zTof), y0=tp->GetY(zTof);
  */
  int IdLc = DCGeomMan::GetInstance().GetLcId();
  double zLc = DCGeomMan::GetInstance().GetLocalZ( IdLc ); 
  double x0=tp->GetX(zLc), y0=tp->GetY(zLc);
  
  int IdSdc3x1 = DCGeomMan::GetInstance().GetDetectorId("SDC3-x-1");
  double zSdc3x1 = DCGeomMan::GetInstance().GetLocalZ(IdSdc3x1);
  double x1=tp->GetX(zSdc3x1), y1=tp->GetY(zSdc3x1);


  ThreeVector pos0 = ThreeVector(x0, y0, 0.0);
  //ThreeVector globalPos0 = DCGeomMan::GetInstance().Local2GlobalPos(IdTof, pos0);
  ThreeVector globalPos0 = DCGeomMan::GetInstance().Local2GlobalPos(IdLc, pos0);

  ThreeVector pos1 = ThreeVector(x1, y1, 0.0);
  ThreeVector globalPos1 = DCGeomMan::GetInstance().Local2GlobalPos(IdSdc3x1, pos1);
#if 0
  std::cout << "Pos0 (x,y,z) = (" 
	    << globalPos0.x() << ", "
	    << globalPos0.y() << ", "
	    << globalPos0.z() << ")" << std::endl;

  std::cout << "Pos1 (x,y,z) = (" 
	    << globalPos1.x() << ", "
	    << globalPos1.y() << ", "
	    << globalPos1.z() << ")" << std::endl;
#endif

  LocalTrackSdcOut_[nTrackSdcOut_] = new TPolyLine3D(2);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineColor(kRed);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetLineWidth(1);
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(0, globalPos0.x(), globalPos0.y(), globalPos0.z());
  LocalTrackSdcOut_[nTrackSdcOut_]->SetPoint(1, globalPos1.x(), globalPos1.y(), globalPos1.z());
  tp_->cd();
  LocalTrackSdcOut_[nTrackSdcOut_]->Draw();
  tc_->Update();
  nTrackSdcOut_++;
}



void EvDisp::DrawSksTrack(int nStep, ThreeVector *StepPoint) const
{
  if (SksStepMark_)
    SksStepMark_->Delete();

  SksStepMark_ = new TPolyMarker3D(nStep);
  for (int i=0; i<nStep; i++) {
    SksStepMark_->SetPoint(i, StepPoint[i].x(), StepPoint[i].y(),  StepPoint[i].z());
  }
  SksStepMark_->SetMarkerSize(1);
  SksStepMark_->SetMarkerColor(kBlue);
  SksStepMark_->SetMarkerStyle(6);

  tp_->cd();
  SksStepMark_->Draw();

  tc_->Update();
}

#if 0
void EvDisp::DrawVertex(ThreeVector vtxPoint, double cost) const
{
  if (nVertex_>=MaxTrack) {
    std::cerr << "EvDisp::DrawVertex nVertex_ is greater than MaxTrack"
	      << std::endl;
    return;
  }

  VtxPoint_[nVertex_] = new TPolyMarker3D(1);
  VtxPoint_[nVertex_]->SetPoint(0, vtxPoint.x(), vtxPoint.y(), vtxPoint.z());
  VtxPoint_[nVertex_]->SetMarkerSize(1);
  if (cost<0.995)
    VtxPoint_[nVertex_]->SetMarkerColor(6); //Pink
  else
    VtxPoint_[nVertex_]->SetMarkerColor(3); //Green
  VtxPoint_[nVertex_]->SetMarkerStyle(29); // Star   

  tp_->cd();
  VtxPoint_[nVertex_]->Draw();
  if (nVertex_==0) {
    tp_->GetView()->ZoomIn();
    tp_->GetView()->ZoomIn();
    tp_->GetView()->ZoomIn();
  }
  tc_->Update();

  tp_vtx_->cd();
  VtxPoint_[nVertex_]->Draw();
  if (nVertex_==0) {
    tp_vtx_->GetView()->ZoomIn();
    tp_vtx_->GetView()->ZoomIn();
    tp_vtx_->GetView()->ZoomIn();
    tp_vtx_->GetView()->ZoomIn();
  }
  tc_vtx_->Update();

  nVertex_++;
}
#endif


void EvDisp::EndOfEvent(void) const
{

  //nTrackBdcOut_=0;
  nTrackSdcIn_=0;
  nTrackSdcOut_=0;
  //nVertex_=0;

  if (InitStepMark_) {
    delete InitStepMark_;
    InitStepMark_ = NULL;
  }
  /*
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackBdcOut_[i]) {
      delete LocalTrackBdcOut_[i];
      LocalTrackBdcOut_[i] = NULL;
    }
  */
  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcIn_[i]) {
      delete LocalTrackSdcIn_[i];
      LocalTrackSdcIn_[i] = NULL;
    }

  for (int i=0; i<MaxTrack; i++)
    if (LocalTrackSdcOut_[i]) {
      delete LocalTrackSdcOut_[i];
      LocalTrackSdcOut_[i] = NULL;
    }
  /*
  for (int i=0; i<MaxTrack; i++)
    if (VtxPoint_[i]) {
      delete VtxPoint_[i];
      VtxPoint_[i] = NULL;
    }
  */
  if (SksStepMark_) {
    delete SksStepMark_;
    SksStepMark_ = NULL;
  }

  ResetVisibility();

}

void EvDisp::ResetVisibility(void) const
{
  /*
  for (int i=0; i<MaxWireBDC; i++) {
    Bdc3x_Node_[i]->SetVisibility(0);
    Bdc3xp_Node_[i]->SetVisibility(0);
    Bdc3u_Node_[i]->SetVisibility(0);
    Bdc3up_Node_[i]->SetVisibility(0);
    Bdc3v_Node_[i]->SetVisibility(0);
    Bdc3vp_Node_[i]->SetVisibility(0);
    Bdc4x_Node_[i]->SetVisibility(0);
    Bdc4xp_Node_[i]->SetVisibility(0);
    Bdc4u_Node_[i]->SetVisibility(0);
    Bdc4up_Node_[i]->SetVisibility(0);
    Bdc4v_Node_[i]->SetVisibility(0);
    Bdc4vp_Node_[i]->SetVisibility(0);
  }
  */

  for (int i=0; i<MaxWireSDC1; i++) {
    Sdc1u1_Node_[i]->SetVisibility(0);
    Sdc1u2_Node_[i]->SetVisibility(0);
    Sdc1v1_Node_[i]->SetVisibility(0);
    Sdc1v2_Node_[i]->SetVisibility(0);

    Sdc1u1_Node_vtx_[i]->SetVisibility(0);
    Sdc1u2_Node_vtx_[i]->SetVisibility(0);
    Sdc1v1_Node_vtx_[i]->SetVisibility(0);
    Sdc1v2_Node_vtx_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC2; i++) {
    Sdc2u1_Node_[i]->SetVisibility(0);
    Sdc2u2_Node_[i]->SetVisibility(0);
    Sdc2v1_Node_[i]->SetVisibility(0);
    Sdc2v2_Node_[i]->SetVisibility(0);
    Sdc2x1_Node_[i]->SetVisibility(0);
    Sdc2x2_Node_[i]->SetVisibility(0);

    Sdc2u1_Node_vtx_[i]->SetVisibility(0);
    Sdc2u2_Node_vtx_[i]->SetVisibility(0);
    Sdc2v1_Node_vtx_[i]->SetVisibility(0);
    Sdc2v2_Node_vtx_[i]->SetVisibility(0);
    Sdc2x1_Node_vtx_[i]->SetVisibility(0);
    Sdc2x2_Node_vtx_[i]->SetVisibility(0);
  }

  for (int i=0; i<MaxWireSDC3V; i++) {
    Sdc3v1_Node_[i]->SetVisibility(0);
    Sdc3v2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC3X; i++) {
    Sdc3x1_Node_[i]->SetVisibility(0);
    Sdc3x2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC3U; i++) {
    Sdc3u1_Node_[i]->SetVisibility(0);
    Sdc3u2_Node_[i]->SetVisibility(0);
  }

  for (int i=0; i<MaxWireSDC4V; i++) {
    Sdc4v1_Node_[i]->SetVisibility(0);
    Sdc4v2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC4X; i++) {
    Sdc4x1_Node_[i]->SetVisibility(0);
    Sdc4x2_Node_[i]->SetVisibility(0);
  }
  for (int i=0; i<MaxWireSDC4U; i++) {
    Sdc4u1_Node_[i]->SetVisibility(0);
    Sdc4u2_Node_[i]->SetVisibility(0);
  }

  /*
  for (int i=0; i<MaxWireBDC; i++) {
    Bdc3x_Node_vtx_[i]->SetVisibility(0);
    Bdc3xp_Node_vtx_[i]->SetVisibility(0);
    Bdc3u_Node_vtx_[i]->SetVisibility(0);
    Bdc3up_Node_vtx_[i]->SetVisibility(0);
    Bdc3v_Node_vtx_[i]->SetVisibility(0);
    Bdc3vp_Node_vtx_[i]->SetVisibility(0);
    Bdc4x_Node_vtx_[i]->SetVisibility(0);
    Bdc4xp_Node_vtx_[i]->SetVisibility(0);
    Bdc4u_Node_vtx_[i]->SetVisibility(0);
    Bdc4up_Node_vtx_[i]->SetVisibility(0);
    Bdc4v_Node_vtx_[i]->SetVisibility(0);
    Bdc4vp_Node_vtx_[i]->SetVisibility(0);
  }
  */
  for (int i=0; i<NumOfSegBH2; i++) {
    Bh2Seg_Node_[i]->SetLineColor(kBlack);
    Bh2Seg_Node_vtx_[i]->SetLineColor(kBlack);
  }

  for (int i=0; i<NumOfSegTOF; i++) 
    TofSeg_Node_[i]->SetLineColor(kBlack);

  for (int i=0; i<NumOfSegLC; i++) 
    LcSeg_Node_[i]->SetLineColor(kBlack);



}


void EvDisp::calcRotMatrix(double TA, double RA1, double RA2, double *rotMat)
{
  double ct1=cos(RA1*Deg2Rad), st1=sin(RA1*Deg2Rad);
  double ct2=cos(RA2*Deg2Rad), st2=sin(RA2*Deg2Rad);
  double ct0=cos(TA*Deg2Rad), st0=sin(TA*Deg2Rad);

  double rotMat1[3][3], rotMat2[3][3];

  /* rotation matrix which is same as DCGeomRecord.cc*/
  rotMat1[0][0] =  ct0*ct2-st0*ct1*st2;
  rotMat1[0][1] = -st0*ct2-ct0*ct1*st2;
  rotMat1[0][2] =  st1*st2;

  rotMat1[1][0] =  ct0*st2+st0*ct1*ct2;
  rotMat1[1][1] = -st0*st2+ct0*ct1*ct2;
  rotMat1[1][2] = -st1*ct2;

  rotMat1[2][0] =  st0*st1;
  rotMat1[2][1] =  ct0*st1;
  rotMat1[2][2] =  ct1;

  /* rotation matrix which rotate -90 deg at x axis*/
  rotMat2[0][0] =  1.0;
  rotMat2[0][1] =  0.0;
  rotMat2[0][2] =  0.0;

  rotMat2[1][0] =  0.0;
  rotMat2[1][1] =  0.0;
  rotMat2[1][2] =  1.0;

  rotMat2[2][0] =  0.0;
  rotMat2[2][1] = -1.0;
  rotMat2[2][2] =  0.0;

  for (int i=0; i<9; i++)
    rotMat[i]=0.0;

  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	//rotMat[3*i+j] += rotMat1[i][k]*rotMat2[k][j];
	rotMat[i+3*j] += rotMat1[i][k]*rotMat2[k][j];
      }
    }
  }

}


void EvDisp::get_command(void) const
{
  char ch;
  char data[100];
  static int stat=0;
  static int Nevent=0;
  static int ev=0;
  
  if (stat == 1 && Nevent > 0 && ev<Nevent) {
    ev++;
    return;
  } 
  if (ev==Nevent) {
    stat=0;
    ev=0;
  }

  if (stat == 0) {
    printf("q|n|p>");

    /* get command */
    scanf("%c",&ch);
    if (ch!='\n')
      while(getchar() != '\n');

    switch (ch) {
    case 'q': exit(0);
    case 'n':
      stat = 1;
      do {
	printf("event#>");
	scanf("%s",data);
      } while ((Nevent=atoi(data))<=0);
      std::cout << "Continue " << Nevent << "event" << std::endl;
      break;
    case 'p':
      theApp->Run(kTRUE);
      break;
    }
  }
}
