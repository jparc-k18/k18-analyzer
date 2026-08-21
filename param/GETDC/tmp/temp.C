

void temp()
{



  std::stringstream rootfilename;
  TF1 *fit1 = new TF1 ("fit1", "[0]*exp([1]*x)+[2]+[3]*x",0,8000);


  //rootfilename<<"~/directory/User_jaeyong/rootfile/HyperballJ0"<<run<<"_p18_"<<date<<".root";
  rootfilename<<"~/directory/User_jaeyong/rootfile/HyperballJ0"<<"8148"<<"_p18_"<<"0608"<<".root";

  TFile* roofi = new TFile(rootfilename.str().c_str(),"R");
  TTree* tr = dynamic_cast<TTree*>(gDirectory->Get("tree"));

  TCanvas *c1 = new TCanvas ("c1","c1",900,900);

  std::stringstream draw1;
  int ch=0;
  draw1<<"GeTdcTFAR["<<ch<<"]:GeAdc["<<ch<<"]>>h(800,0,8000,500,0,5000)";

  tr -> Draw( draw1.str().c_str() );
  fit1->SetParameters(-400.,-0.02,4350.,0.0002);
  h->Draw("colz");
  //h->Fit("fit1","","",100,7500);
  //c1->Update();






}
