//Parameters

// Momentum offset
const double pB_offset = 1.00;
const double pS_offset = 0.00;
const double pK18_offset = 0.00;
//const double pSKS_offset = 0.0;
//const double pSKS_offset = -0.0122897;
const double pSKS_offset = 0.0;

// Position offset
const double x_off = -2.8739;
const double y_off = -4.50123;
const double u_off =-0.0000924759;
const double v_off =-0.00180667;

// BH2 Cut Parameters
const double MinDeltaEBH2 = 0.5; //for Pion
const double MaxDeltaEBH2 = 3.0; //for Pion

// BH1 Cut Parameters
const double MinDeltaEBH1 = 0.5; //for Pion
const double MaxDeltaEBH1 = 3.0; //for Pion

// BH1-BH2
const double MinBeamToF  = -5.0;
const double MaxBeamToF  =  5.0;

// BH2 - target offset
const double offset =1.3108917166;

// TOF Cut Parameters
const double MinDeltaETof = -100.; //for Kaon
const double MaxDeltaETof = 100.0; //for Kaon
const double MinTimeTof =  -100.0;
const double MaxTimeTof =  100.0;

// Muon Cut Parameters
const double MinDeltaELc = -100.0; //for Kaon
const double MaxDeltaELc = 100.0; //for Kaon
const double MinTimeLc =  -100.0;
const double MaxTimeLc =  100.0;

//SdcIn vs TOF&LC cut
const double ModFactorTof  = 0;
const double ModFactorLc   = 0;
const double MinModTimeTof = 0.;
const double MaxModTimeTof = 0.;
const double MinModTimeLc  = 0.;
const double MaxModTimeLc  = 0.;

// TOF-SdcOut
const double MinXDifTof  =  -100.;
const double MaxXDifTof  =   100.;
const double MinYDifTof  =  -100.;
const double MaxYDifTof  =   100.;

// LC-SdcOut
const double MinXDifLc  =  -100.;
const double MaxXDifLc  =   100.;
const double MinYDifLc  =  -100.;
const double MaxYDifLc  =   100.;

//Multi cut of DCs
const double MaxMultiHitBcIn   = 5.0;
const double MaxMultiHitBcOut  = 15.;
const double MaxMultiHitSdcIn  = 10.;
const double MaxMultiHitSdcOut = 3.0;

const double MaxChisqrBcIn   = 20.;
const double MaxChisqrBcOut  = 20.;
const double MaxChisqrSdcIn  = 20.;
const double MaxChisqrSdcOut = 20.;
const double MaxChisqrK18    = 30.;
const double MaxChisqrSks    = 30.;

const int MaxNTrackBcIn   = 30;
const int MaxNTrackBcOut  = 10;
const int MaxNTrackSdcIn  = 10;
const int MaxNTrackSdcOut = 20;
const int MaxNTrackK18    = 20;
const int MaxNTrackSks    = 20;

const double MinTgtXSdcIn  =  -100.;
const double MaxTgtXSdcIn  =   100.;
const double MinTgtYSdcIn  =  -50.;
const double MaxTgtYSdcIn  =   50.;

const double MinTgtXBcOut  =  -100.;
const double MaxTgtXBcOut  =   100.;
const double MinTgtYBcOut  =  -50.;
const double MaxTgtYBcOut  =   50.;

const double MinMassSquare =  -2.;
const double MaxMassSquare =  2.;

const double MinVertexZ = -200.0;
const double MaxVertexZ =  200.0;

//Lorents Vector Constant
const double AtomicMassUnit  = 0.93149432;
const double PionMass        = 0.1395701;
const double KaonMass        = 0.493677;
const double ProtonMass      = 0.93827200;
const double NeutronMass     = 0.93956563;
const double CarbonMass      = 12.*AtomicMassUnit;
const double DeutronMass     = 2.*AtomicMassUnit+0.01313672;
const double DeltaC11        = 0.0106502;
const double CoreMass        = 11.*AtomicMassUnit+DeltaC11;
const double LambdaMass      = 1.115648;
const double SigmaPMass      = 1.18937;
const double ThetaMass       = 1.530;
const double MuonMass        = 0.1056583715;   
const double He4Mass         = 4.0026032*AtomicMassUnit;
const double F19Mass         = 18.9984032*AtomicMassUnit;
const double He3Mass         = 3.0160293*AtomicMassUnit;
const double F18Mass         = 18.0009377*AtomicMassUnit;

const double TargetMass1 = ProtonMass;
const double TargetMass2 = CarbonMass;
const double TargetMass3 = F19Mass;
const double TargetMass4 = He4Mass;

const double CoreMass1 = ProtonMass;
const double CoreMass2 = CoreMass;
const double CoreMass3 = F18Mass;
const double CoreMass4 = He3Mass;

const double BLambdaC12 = -0.011;
const double BLambdaF19 = -0.011;
const double BLambdaHe4 = -0.011;

const double BSigmaC12 = -0.011;
const double BSigmaF19 = -0.011;

const double CLMass = CoreMass + LambdaMass + BLambdaC12;
const double FLMass = F18Mass + LambdaMass + BLambdaF19;
const double HeLMass = He3Mass + LambdaMass + BLambdaHe4;
 
const double CSMass = CoreMass + SigmaPMass + BSigmaC12;
const double FSMass = F18Mass + SigmaPMass + BSigmaF19;


