

void fillTdcParam()
{
  std::fstream oldfile, newfile, outfile;
  std::string oldfilename, newfilename, outfilename;

  //### define old parameter file ###
  std::cout<<"input old file name :"<<std::endl;
  std::cin>>oldfilename;
  std::cout<<oldfilename.c_str()<<std::endl;

    
  oldfile.open(oldfilename.c_str(),std::fstream::in);

  int ch_flag[32]={0};
  double oldpar[32][20];

  if(oldfile.is_open())
    {
      std::cout<<"Read old paramters ..."<<std::endl;
      for(int i=0;i<32;i++)
	{
	  oldfile>>ch_flag[i];
	  for(int j=0;j<20;j++)
	    {
	      oldfile>>oldpar[i][j];
	    }
	}
      oldfile.close();
      std::cout<<"Finish reading old paramters ..."<<std::endl;
    }
  
  //### define new parameter file ###
  std::cout<<"input new file name :"<<std::endl;
  
  std::cin>>newfilename;
  //std::cout<<newfilename.c_str()<<std::endl;
  
  //oldfile.open("GeTdc.param.e13.20150607_v2",std::fstream::in);
  newfile.open(newfilename.c_str(),std::fstream::in);

  int ch_flag[32]={0};
  double newpar[32][20];

  if(newfile.is_open())
    {
      std::cout<<"Read new paramters ..."<<std::endl;
      for(int i=0;i<32;i++)
	{
	  newfile>>ch_flag[i];
	  for(int j=0;j<8;j++)
	    {
	      newfile>>newpar[i][j];
	    }
	}
      newfile.close();
      std::cout<<"Finish reading new paramters ..."<<std::endl;
    }

  //### define new parameter file ###
  std::cout<<"input output file name :"<<std::endl;
  std::cin>>outfilename;
  outfile.open(outfilename.c_str(),std::fstream::out|std::fstream::app);

  
  int num = 0;
  char flag='n';
  
  do{
    if(flag=='y')
      std::cout<<"Input error!! Please select a function again!!"<<std::endl;
    
    std::cout<<"Select parameters you deal with: "<<std::endl;
    std::cout<<"1. CRM paramter"<<std::endl;
    std::cout<<"2. TFA paramter"<<std::endl;
  
    std::cin>>num;
    flag= 'y';
  }while( num<1 || num>3 );
  
  switch(num)
    {
    case 1: 
      if(outfile.is_open())
	{
	  std::cout<<"Write updated CRM parameters ..."<<std::endl;
	  for(int i=0; i<32;i++)
	    {
	      outfile<<ch_flag[i];
	      for(int j=0;j<8;j++)
		{
		  outfile<<" "<<newpar[i][j];
		}
	      for(int j=8;j<20;j++)
		{
		  outfile<<" "<<oldpar[i][j];
		}
	      outfile<<std::endl;
	    }      
	  outfile.close();
	  std::cout<<"Output file is produced ..."<<std::endl;
	}
      break;
    case 2:
      if(outfile.is_open())
	{
	  std::cout<<"Write updated TFA parameters ..."<<std::endl;
	  for(int i=0; i<32;i++)
	    {
	      outfile<<ch_flag[i];
	      for(int j=0;j<12;j++)
		{
		  outfile<<" "<<oldpar[i][j];
		}
	      for(int j=0;j<8;j++)
		{
		  outfile<<" "<<newpar[i][j];
		}
	      outfile<<std::endl;
	    }      
	  outfile.close();
	  std::cout<<"Output file is produced ..."<<std::endl;
	}
      break;
    }

}
