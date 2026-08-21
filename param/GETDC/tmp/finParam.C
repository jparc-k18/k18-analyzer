

void finParam()
{

  std::fstream prefile, outfile;
  std::stringstream prefilename, outfilename;
  std::string run;  
  std::cout<<"Input rum #:0xxxx"<<std::endl;
  std::cin>>run;
  getchar(); 
  int num = 0;
  char flag='n';

  do{
    if(flag=='y')
      std::cout<<"Input error!! Please select a function again!!"<<std::endl;

    std::cout<<"Select parameters you deal with: "<<std::endl;
    std::cout<<"1. CRM paramter"<<std::endl;
    std::cout<<"2. TFA paramter"<<std::endl;

    std::cin>>num;
    getchar();
    flag= 'y';
  }while( num<1 || num>3 );

  switch(num)
    {
    case 1:
      prefilename<<"preParam0"<<run<<".CRM.param";
      outfilename<<"finParam0"<<run<<".CRM.param";
      break;
    case 2:
      prefilename<<"preParam0"<<run<<".TFA.param";
      outfilename<<"finParam0"<<run<<".TFA.param";
      break;
    }

  prefile.open(prefilename.str().c_str(),std::fstream::in);

  int nch =32;
  int ch_flag[400], s[400];
 
  double x[400], y1[400], y2[400];
  std::cout<<"pre-parameter readout starts!!"<<std::endl;
  for(int i=0;i<7*nch;i++)
    {
      prefile>>ch_flag[i]>>x[i]>>y1[i]>>y2[i];
    }
  std::cout<<"pre-parameter readout finished!!"<<std::endl;

  TCanvas *c1 = new TCanvas ("c1","c1",600,600); 
  
  //int ch=23; 
  for( int ch=0;ch<32;ch++)
    { 
      std::cout<<"ch "<<ch<<" fitting.."<<std::endl;

      std::cout<<x[7+7*ch]<<" "<<y1[7+7*ch]<<" "<<y2[7+7*ch]<<std::endl;;
      std::cout<<x[8+7*ch]<<" "<<y1[7+7*ch]<<" "<<y2[8+7*ch]<<std::endl;;
      double p[8]={0};      
      //if(ch_flag[i]==ch)
      TF1 *fit1 = new TF1 ("fit1", "[0]*exp([1]*x)+[2]+[3]*x",0,7500);
      TF1 *fit2 = new TF1 ("fit2", "[0]*exp([1]*x)+[2]+[3]*x",0,7500);

      TGraph* gr1 = new TGraph(7 , x+7*ch , y1+7*ch);
      TGraph* gr2 = new TGraph(7 , x+7*ch , y2+7*ch);
      
      

      fit1->SetParameters(-300.,-0.005,4300.,0.002);
      gr1->Fit("fit1","","",0,8000);
      fit2->SetParameters(-150.,-0.02,4400.,0.002);
      gr2->Fit("fit2","","",0,8000);
      
      TMultiGraph* mg = new TMultiGraph;

      mg->Add(gr1);
      mg->Add(gr2);
      mg->Draw("a");
      //gr1->Draw("AP");
      //gr2->SetMarkerColor(2);
      //gr2->Draw("P");

      c1->Update();
      prefile.close();
      std::cout<<"ch "<< ch <<" fitting finished!!"<<std::endl;
      getchar();
      //fit1->Delete();
      //gr1->Delete();
      //gr2->Delete();
      double par[8]={0};
      fit1->GetParameters(par);
      fit2->GetParameters(par+4);
      par[0]= par[0]-100;
      par[1]=par[1]*0.8;
      par[4]=par[4]+30;

      //if(ch==3||ch==4||ch==11||ch==13||ch==15||ch==29||ch==31) //~0609
      if(ch==3||ch==4||ch==12||ch==13||ch==15||ch==29||ch==31) //0610~
	{
	  for(int k=0;k<8;k++)
	    par[k]=0;
	}
      outfile.open(outfilename.str().c_str(),std::fstream::out | std::fstream::app );  
      if(outfile.is_open())
	{
	  std::cout<<"Outfile open!!"<<std::endl;
	  //std::cout<<"par[0]: "<<par[0]<<std::endl;

	  outfile<<ch;

	  for(int k=0;k<8;k++){
	    outfile<<" "<<par[k];
	  }

	  outfile<<std::endl;
	  outfile.close();	
	}
      std::cout<<"ch "<< ch <<" parameters are written!!"<<std::endl;
    }
  if(num==1)
    std::cout<<"run0"<<run.c_str()<<"_CRM finishied!!"<<std::endl;
  else
    std::cout<<"run0"<<run.c_str()<<"_TFA finishied!!"<<std::endl;

}
