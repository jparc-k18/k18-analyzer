

void autoParam0611()
{
 
  

  TCanvas *c1 = new TCanvas ("c1","c1",600,600);


  //c1->Divide(2,1);

  std::fstream file, outfile;
  file.open("./pre_auto0611.param",std::fstream::in);

  int n =12;
  int ch_flag[400], s[400];
  int ch[25]={0,1,2,5,6,7,8,9,10,11,14,16,17,18,19,20,21,22,23,24,25,26,27,28,30};  
  //double slot[32];
 
  double x[400], y1[400], y2[400];
  std::cout<<"pre-parameter readout starts!!"<<std::endl;
  for(int i=0;i<7*25;i++)
    {
      file>>ch_flag[i]>>x[i]>>y1[i]>>y2[i];
    }
  std::cout<<"pre-parameter readout finished!!"<<std::endl;
  //int j=23;
  
  for( int j=0;j<25;j++)
    { 
      std::cout<<"ch "<<j<<" fitting.."<<std::endl;

      std::cout<<x[7+7*j]<<" "<<y1[7+7*j]<<" "<<y2[7+7*j]<<std::endl;;
      std::cout<<x[8+7*j]<<" "<<y1[7+7*j]<<" "<<y2[8+7*j]<<std::endl;;
      double p[8]={0};      
      //if(ch_flag[i]==j)
      TF1 *fit1 = new TF1 ("fit1", "[0]*exp([1]*x)+[2]+[3]*x",0,7500);
      TF1 *fit2 = new TF1 ("fit2", "[0]*exp([1]*x)+[2]+[3]*x",0,7500);

      TGraph* gr1 = new TGraph(7 , x+7*j , y1+7*j);
      TGraph* gr2 = new TGraph(7 , x+7*j , y2+7*j);
      
      

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
      file.close();
      std::cout<<"ch "<<ch[j]<<" fitting finished!!"<<std::endl;
      getchar();
      //fit1->Delete();
      //gr1->Delete();
      //gr2->Delete();
      double par[8]={0};
      fit1->GetParameters(par);
      fit2->GetParameters(par+4);

      outfile.open("./TDCauto0611.param",std::fstream::out | std::fstream::app );  
      if(outfile.is_open())
	{
	  std::cout<<"Outfile open!!"<<std::endl;
	  std::cout<<"par[0]: "<<par[0]<<std::endl;
	  //outfile<<"slot :"<<j<<std::endl;
	  outfile<<ch[j]<<" "<<par[0]<<" "<<par[1]<<" "<<par[2]<<" "<<par[3]/2<<" "<<par[4]<<" "<<par[5]<<" "<<par[6]<<" "<<par[7]<<std::endl;
	  outfile.close();	
	}
      std::cout<<"ch "<<ch[j]<<" parameters are written!!"<<std::endl;
    }
  

}
