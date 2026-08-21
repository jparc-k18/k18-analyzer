
#include "DopplerCorr.hh"
#include "Parameter.hh"
#include <cmath>


const double HLMass = FLMass;


DopplerCorr::DopplerCorr()
{
  ThreeVector null(0,0,0);
  Vertex=null;
  RecoilMom=null;
}

void DopplerCorr::Initialize( ThreeVector recp, ThreeVector vert )
{
  Vertex=vert;
  RecoilMom=recp;
}

double DopplerCorr::CalcDopplerFactor( const int slot )
{
  ThreeVector Offset( 0.0, 0.0, 0.0 );
  ThreeVector geposi=GetGePosition(slot)+Offset;
  ThreeVector gamdir=geposi-Vertex;
  double gamcost=(gamdir*RecoilMom)/(gamdir.Mag()*RecoilMom.Mag());
  double recoilE=sqrt(RecoilMom.Mag()*RecoilMom.Mag()+HLMass*HLMass);
  double recoilbeta=RecoilMom.Mag()/recoilE;
  double gam=1/(sqrt(1-recoilbeta*recoilbeta));

  return gam*(1-recoilbeta*gamcost);
}

double DopplerCorr::CalcGtheta(const int slot)
{
  ThreeVector Offset( 0.0, 0.0, 0.0 );
  ThreeVector geposi=GetGePosition(slot)+Offset;
  ThreeVector gamdir=geposi-Vertex;
  double gamcost=(gamdir*RecoilMom)/(gamdir.Mag()*RecoilMom.Mag());
  return acos(gamcost)*180./acos(-1.);
}

double DopplerCorr::CalcBeta(void)
{
  double recoilE=sqrt(RecoilMom.Mag()*RecoilMom.Mag()+HLMass*HLMass);
  double recoilbeta=RecoilMom.Mag()/recoilE;
  return recoilbeta;
}


ThreeVector DopplerCorr::GetGePosition(const int slot)
{
  ThreeVector geposi;

  // Geposition 
  if(slot==0){
    ThreeVector gepos( -129.5, 105.0, -194.5 );
    geposi=gepos;
  }
  if(slot==1){
    ThreeVector gepos( 0.0, 155.0, -194.5 );
    geposi=gepos;
  }
  if(slot==2){
    ThreeVector gepos( 129.5, 105., -194.5 );
    geposi=gepos;
  }
  if(slot==3){
    ThreeVector gepos( -258., 35.0, -65. );
    geposi=gepos;
  }	  
  if(slot==4){
    ThreeVector gepos( -258., 35.0, 65. );
    geposi=gepos;
  }
  if(slot==5){
    ThreeVector gepos( -129., 105., -65. );
    geposi=gepos;
  }
  if(slot==6){
    ThreeVector gepos( -129., 105.0, 65. );
    geposi=gepos;
  }
  if(slot==7){
    ThreeVector gepos( 0.0, 175.0, -65. );
    geposi=gepos;
  }
  if(slot==8){
    ThreeVector gepos( 0.0, 175.0, 65.0 );
    geposi=gepos;
  }
  if(slot==9){
    ThreeVector gepos( 129.0, 105.0, -65. );
    geposi=gepos;
  }
  if(slot==10){
    ThreeVector gepos( 129.0, 105.0, 65.0 );
    geposi=gepos;
  }
  if(slot==11){
    ThreeVector gepos( 258.0, 35.0, -65.0 );
    geposi=gepos;
  }
  if(slot==12){
    ThreeVector gepos( 258.0, 35., 65.0);
    geposi=gepos;
  }
  if(slot==13){
    ThreeVector gepos( -129.5, 105., 194.5 );
    geposi=gepos;
  }
  if(slot==14){
    ThreeVector gepos( 0.0, 155.0, 194.5 );
    geposi=gepos;
  }
  if(slot==15){
    ThreeVector gepos( 129.5, 105.0, 194.5 );
    geposi=gepos;
  }
  if(slot==16){
    ThreeVector gepos( -129.5, -105.0, -194.5 );
    geposi=gepos;
  }
  if(slot==17){
    ThreeVector gepos( 0.0, -155.0, -194.5 );
    geposi=gepos;
  }
  if(slot==18){
    ThreeVector gepos( 129.5, -105., -194.5 );
    geposi=gepos;
  }
  if(slot==19){
    ThreeVector gepos( -258., -35.0, -65. );
    geposi=gepos;
  }	  
  if(slot==20){
    ThreeVector gepos( -258., -35.0, 65. );
    geposi=gepos;
  }
  if(slot==21){
    ThreeVector gepos( -129., -105., -65. );
    geposi=gepos;
  }
  if(slot==22){
    ThreeVector gepos( -129., -105.0, 65. );
    geposi=gepos;
  }
  if(slot==23){
    ThreeVector gepos( 0.0, -175.0, -65. );
    geposi=gepos;
  }
  if(slot==24){
    ThreeVector gepos( 0.0, -175.0, 65.0 );
    geposi=gepos;
  }
  if(slot==25){
    ThreeVector gepos( 129.0, -105.0, -65. );
    geposi=gepos;
  }
  if(slot==26){
    ThreeVector gepos( 129.0, -105.0, 65.0 );
    geposi=gepos;
  }
  if(slot==27){
    ThreeVector gepos( 258.0, -35.0, -65.0 );
    geposi=gepos;
  }
  if(slot==28){
    ThreeVector gepos( 258.0, -35.0, 65.0);
    geposi=gepos;
  }
  if(slot==29){
    ThreeVector gepos( -129.5, -105., 194.5 );
    geposi=gepos;
  }
  if(slot==30){
    ThreeVector gepos( 0.0, -155.0, 194.5 );
    geposi=gepos;
  }
  if(slot==31){
    ThreeVector gepos( 129.5, -105.0, 194.5 );
    geposi=gepos;
  }
  
  return geposi;
}

