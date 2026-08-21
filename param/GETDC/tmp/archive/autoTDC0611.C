#include <fstream>
#include <iostream>
#include <sstream>


void autoTDC0611()
//int main()
{
  std::fstream file;
  std::stringstream filename;
 
  file.open("./pre_auto0611.param",std::fstream::out | std::fstream::app );
  std::string run="8225";
  std::string date="0611";
  filename<<"~/directory/User_jaeyong/rootfile/HyperballJ0"<<run<<"_p18_"<<date<<".root";

  TFile* roofi = new TFile(filename.str().c_str(),"R");
  TTree* tr = dynamic_cast<TTree*>(gDirectory->Get("tree"));

  int nseg=6;
  int nch=32;  
  int r[7]={0,250,600,1500,3000,5500,8000};

  double fr=0;
  double fl=0,fh=0;
 
  double fr[7];
  double sl[7];
  double sh[7];


  int ch,j;
  //std::cin>>ch;
  //std::cin>>j;


  TCanvas *c1 = new TCanvas ("c1","c1",600,900);

  c1->Divide(2,3);
  //ch=27;
  //ch=10;
  for(ch=0;ch<nch;ch++)
    {
      std::cout<<"ch "<< ch <<" start!!! "<<std::endl;
      if(ch==3||ch==4||ch==12||ch==13||ch==15||ch==29||ch==31) 
	{
	  std::cout<<"ch "<< ch <<" is empty!!! "<<std::endl;
	  continue;
	}
      else
	{
	  for(j=0;j<nseg;j++)
	    {

	      std::stringstream draw, draw1, cut1, hist;
	      double mean[14]={0};
	      double sig[14]={0};  

	      //draw1<<"GeTdcTFAR["<<ch<<"]>>h"<<j<<"(100,3500,4500)";
	      draw1<<"GeTdcCRMR["<<ch<<"]>>h"<<j<<"(100,3500,4500)";

	      cut1<<"GeAdc["<<ch<<"]>"<<r[j]<<"&&GeAdc["<<ch<<"]<"<<r[j+1];
	      hist<< "h"<<j;

	      c1->cd(j+1);
	      TCut cut = cut1.str().c_str();

	      TF1 *f1 = new TF1("f1","gaus",0,8000);    
	    
	      tr->Draw( draw1.str().c_str(), cut );
  
	      TH1D* h = dynamic_cast<TH1D*>(gDirectory->Get(hist.str().c_str()));
	       
	      int binmax = h->GetMaximumBin();
	      double x = h->GetXaxis()->GetBinCenter(binmax);
      
	      std::cout<<"ch: "<<ch<<"  seg:"<<j;
	      if(j==0){
		fr=100;
		fl=5;
		fh=1;	
	      }
	      else if(j==1){
		fr=50;
		fl=4;
		fh=1;	
	      }
	      else if(j==2){
		fr=30;
		fl=4;
		fh=2;
	      }
	      else if(j==3){
		fr=30;
		fl=4;
		fh=2;
	      }
	      else{
		fr=40;
		fl=3;
		fh=3;
	      }
	      
	      //h->Fit("f1","","",x-fr,x+fr);
	      h->Fit("f1","Q","",x-1.5*fr,x+fr);
	      mean[j]=f1->GetParameter(1);
	      sig[j]=f1->GetParameter(2);
	      
	      double ymax = h->GetMaximum();
	      TLine* line1 = new TLine(mean[j]-fl*sig[j],0,mean[j]-fl*sig[j] ,ymax);
	      //TLine* line2 = new TLine(mean[j]+0.8*f*sig[j],0,mean[j]+0.8*f*sig[j] ,ymax);
	      TLine* line2 = new TLine(mean[j]+fh*sig[j],0,mean[j]+fh*sig[j] ,ymax);
	      line1->SetLineColor(3);
	      line2->SetLineColor(3);

	      line1->Draw("same");
	      line2->Draw("same");

	      c1->Update();	      
	      //std::cout<<""<<""<<""<<""<<std
	      //draw.str(std::string());
	      //draw1.str(std::string());
	      //cut1.str(std::string());
	      //hist.str(std::string());	      
	      /*
		if(file.is_open())
		{
		file<<10<<" "<<11<<""<<std::endl;
		}

	      */
	      //std::cout<<ch<<" "<<(r[j]+r[j+1])/2 <<" "<<mean[j]-3*sig[j]<<" "<<mean[j]+3*sig[j]<<std::endl;
	      file<<ch<<" "<<r[j]<<" "<<mean[j]-fl*sig[j]<<" "<<mean[j]+fh*sig[j]<<std::endl;
	      if(j==nseg-1)
	        file<<ch<<" "<<r[j+1]<<" "<<mean[j]-fl*sig[j]<<" "<<mean[j]+fh*sig[j]<<std::endl;
	     
	      //f1->Delete();
	      //h->Delete();
	      std::cout<<"     "<<"fitting done!"<<std::endl;	      
	    }
	}    
      std::cout<<"ch "<< ch <<" is finished!!! "<<std::endl;
      //c1->Clear();
      getchar();
    }
  

  filename.str(std::string());
  
  //file.close();
  //tr->Delete();
  //  roofi->Delete();
  //return 0;

}
