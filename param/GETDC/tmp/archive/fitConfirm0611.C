
void fitConfirm0611()
{
  std::fstream parfile;
  parfile.open("TDCauto0611.param",std::fstream::in);

  int ch_flag[32]={0};;
  double par[32][9]={0};

  if(parfile.is_open())
    {
      std::cout<<"Reading parameter file ..."<<std::endl;
      for(int i=0;i<32;i++)
      {
	parfile>>ch_flag[i]>>par[i][0]>>par[i][1]>>par[i][2]>>par[i][3]>>par[i][4]>>par[i][5]>>par[i][6]>>par[i][7];
      }
      parfile.close();
    }
  std::string run="8225";
  std::string date="0611";
  std::stringstream rootfile;
  rootfile<<"~/directory/User_jaeyong/rootfile/HyperballJ0"<<run<<"_p18_"<<date<<".root";

  TFile* roofi = new TFile(rootfile.str().c_str(),"R");
  TTree* tr = dynamic_cast<TTree*>(gDirectory->Get("tree"));

  TCanvas *c1 = new TCanvas ("c1","c1",1200,1200);
  //c1->Divide(5,5);

  //int ch=10;
  for(int ch=0;ch<32;ch++)
    {
      std::cout<<"ch "<<ch<<" start!!"<<std::endl;
      if(ch==3||ch==4||ch==12||ch==13||ch==15||ch==29||ch==31)
        {
	  std::cout<<"ch "<< ch <<" is empty!!! "<<std::endl;
          continue;
        }
      else
	{
	  std::stringstream hist, uline, dline;
	  std::stringstream hi, dl, ul;
	  //hist<<"GeTdcTFAR["<<ch<<"]:GeAdc["<<ch<<"]>>h"<<ch<<"(8000,0,8000,1000,3500,4500)";
	  hist<<"GeTdcCRMR["<<ch<<"]:GeAdc["<<ch<<"]>>h"<<ch<<"(8000,0,8000,1000,3500,4500)";
	  dline<<par[ch][0]<<"*TMath::Exp("<<par[ch][1]<<"*GeAdc["<<ch<<"])+"<<par[ch][2]<<"+"<<par[ch][3]<<"*GeAdc["<<ch<<"]:GeAdc["<<ch<<"]>>d"<<ch<<"(8000,0,8000,1000,3500,4500)";
	  uline<<par[ch][4]<<"*TMath::Exp("<<par[ch][5]<<"*GeAdc["<<ch<<"])+"<<par[ch][6]<<"+"<<par[ch][7]<<"*GeAdc["<<ch<<"]:GeAdc["<<ch<<"]>>u"<<ch<<"(8000,0,8000,1000,3500,4500)";
	  hi<<"h"<<ch;
	  dl<<"d"<<ch;
	  ul<<"u"<<ch;
	  
	  //c1->cd(ch+1);
	  //std::cout<<"after cd"<<std::endl;
	  tr -> Draw ( hist.str().c_str() );
	  //std::cout<<hist.str().c_str()<<std::endl;
	  tr -> Draw ( dline.str().c_str() ,"" , "goff");
	  //std::cout<<dline.str().c_str()<<std::endl;
	  tr -> Draw ( uline.str().c_str() ,"" , "goff");
	  //std::cout<<uline.str().c_str()<<std::endl;
	  //std::cout<<"after draw"<<std::endl;
	  //TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hist.str().c_str()));
	  TH1D* d = dynamic_cast<TH1D*>(gDirectory->Get(dl.str().c_str()));
	  TH1D* u = dynamic_cast<TH1D*>(gDirectory->Get(ul.str().c_str()));
	  d->SetMarkerColor(2);
	  u->SetMarkerColor(2);
	  d->Draw("same");
	  u->Draw("same");
	  c1->Update();

	}
      std::cout<<"ch "<<ch<<" finished!!"<<std::endl;
      getchar();


      
      
    }

}
