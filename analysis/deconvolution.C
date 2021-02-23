// this does 3W region analysis for epInelastic
// [process][ring][sector]
const int nMatrix=2;
const int nProcDef=10;
const int nProc=3;
double A[nProc][6][3],rate[nProc][6][3],sigma[nProc][6][3];
string procNm[nProcDef]={
  "moller",
  "epElastic",
  "epInelasticW1",
  "epInelasticW2",
  "epInelasticW3",
  "eAlElastic",
  "eAlQuasielastic",
  "eAlInelastic",
  "piminus",
  "neutralBknd"
};
string procNmComb[nProc]={
  "moller",
  "Bknd",
  "neutralBknd"
};

int verbose = 0;

const double neutralBkgndFactor = 1.;
const double neutralBkgndRate[6]={1e9,1e9,1e9,1e9,1e9,1e9};
//const double neutralBkgndRate[6]={5e7,8e7,11e7,7e7,33e7,3e7};

void analyzeOne(int ring, int sect);
void readSim(string fnm,int proc,int addProc);
void printAll();

void deconvolution(){

  for(int i=0;i<nProc;i++)
    for(int j=0;j<6;j++)
      for(int k=0;k<3;k++){
        A[i][j][k]=0;
        rate[i][j][k]=0;
        sigma[i][j][k]=0;
      }

  /// optimal positioning
  string fnms[nProcDef]={
    "histos/deconv_ee_offSet0_basicAnaV0.root",
    "histos/deconv_epE_offSet0_basicAnaV0.root",
    "histos/inelasticout_bkgAnaV4.root",
    "histos/inelasticout_bkgAnaV4.root",
    "histos/inelasticout_bkgAnaV4.root",
    "histos/deconv_eAlE_offSet0_basicAnaV0.root",
    "histos/deconv_eAlQ_offSet0_basicAnaV0.root",
    "histos/deconv_eAlI_offSet0_basicAnaV0.root",
    "histos/deconv_pi_offSet0_basicAnaV0.root",
    "byHand"
  };

  int addProc(0);
  for(int i=0;i<nProc;i++){
    readSim(fnms[i+addProc],i,addProc);
    if(i==1)
      while(addProc<7){
        addProc++;
        readSim(fnms[i+addProc],i,addProc);
      }
  }

  for(int i=0;i<6;i++)
    for(int j=0;j<3;j++)
      analyzeOne(i,j);

  printAll();
}

void analyzeOne(int ring, int sect){
  cout<<"\n\n\nprocessing Ring "<<ring<<" Sector "<<sect<<endl;

  TMatrixD F(nMatrix, nMatrix);
  for(int i=0;i<nMatrix;i++)
    for(int j=0;j<nMatrix;j++)
      F(i,j)=0;
  double B[nMatrix]={0};
  if(verbose){
    cout<<"\tAsym and rate\n";
    for(int i=0;i<nProc;i++)
      cout<<"\t\t"<<procNmComb[i]<<"\t"<<A[i][ring][sect]<<"\t"<<rate[i][ring][sect]<<endl;
  }

  const double polarization = 0.8;
  const double beamDays = 235 + 95 + 14;
  const double days2seconds = 24*60*60;

  for(int i=0;i<6;i++)
    for(int j=0;j<3;j++){
      double rateTot(0);
      //for(int k=0;k<nMatrix;k++)
      for(int k=0;k<nProc;k++)
        rateTot += rate[k][i][j];
      double sigmaAm = 1/sqrt(rateTot);

      double Ami(0),f[nMatrix];
      for(int k=0;k<nMatrix;k++){
	if(rate[k][i][j]==0 || A[k][ring][sect]==0){
	  f[k]=0;
	  if(verbose)
	    cout<<"\ti/j/k/rate\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<rate[k][i][j]<<"\t"
		<<A[k][ring]<<" "<<f[k]<<endl;
	}else{
	  f[k] = rate[k][i][j]/rateTot * A[k][i][j]/A[k][ring][sect];
	  if(verbose)
	    cout<<"\ti/j/k/rate**\t"<<i<<"\t"<<j<<"\t"<<k<<"\t"<<rate[k][i][j]<<"\t"<<f[k]<<endl;
	}
        Ami += f[k]*A[k][ring][sect];
	if(verbose)
	  cout<<"\t i/j/f/f*A/Ami\t"<<i<<"\t"<<j<<"\t"<<f[k]<<"\t"<<f[k]*A[k][ring][sect]
	      <<"\t"<<Ami<<endl;
      }

      for(int k1=0;k1<nMatrix;k1++){
        for(int k2=0;k2<nMatrix;k2++)
          F(k1,k2) += f[k1]*f[k2]/(sigmaAm*sigmaAm);
        B[k1] += Ami * f[k1]/(sigmaAm*sigmaAm);
      }
    }

  cout<<"\nBefore Inverstion\n\tDet: "<<F.Determinant()<<endl;
  cout<<"\tF matrix\n\t";
  for(int i=0;i<nMatrix;i++){
    for(int j=0;j<nMatrix;j++)
      cout<<"\t"<<F(i,j);
    cout<<"\n\t";
  }
  F.Invert();
  cout<<"\nAfter Inverstion\n\tDet: "<<F.Determinant()<<endl;
  cout<<"\tF matrix\n\t";
  for(int i=0;i<nMatrix;i++){
    for(int j=0;j<nMatrix;j++)
      cout<<"\t"<<F(i,j);
    cout<<"\n\t";
  }
  if(verbose){
    cout<<"B vector\n\t";
    for(int j=0;j<nMatrix;j++)
	cout<<"\t"<<B[j];
    cout<<endl;
    cin.ignore();
  }

  cout<<endl<<"\t\t\tAsymmetry extractor\n";
  for(int i=0;i<nMatrix;i++){
    double asym(0);
    for(int j=0;j<nMatrix;j++)
      asym += F(i,j) * B[j];
    cout<<"Asymmetry "<<procNmComb[i]<<"\t"<<asym<<endl;
  }

  cout<<endl<<endl<<"\t\t\toverall\nName\tAsymmetry\tuncert[ppb]\trelative uncer[ppb]\n";
  for(int i=0;i<nMatrix;i++){
    sigma[i][ring][sect] = sqrt( F(i,i) ) / ( 0.8 * sqrt(beamDays * days2seconds) ) * 1e9;
    cout<<procNmComb[i]<<"\t"<<A[i][ring][sect]<<"\t"<<sigma[i][ring][sect]<<"\t"<<sigma[i][ring][sect]/A[i][ring][sect]<<endl;
  }

  cout<<endl<<"\t\t\tin the selected bin\n";
  double totRate(0);
  for(int i=0;i<nProc;i++)
    totRate += rate[i][ring][sect];
  double stat = totRate / rate[0][ring][sect] /
    ( sqrt(totRate) * polarization * sqrt(beamDays * days2seconds)) *1e9;
  cout<<"stat moller\t"<< stat << "\t" << stat/33<<endl;

  double syst[nProc];
  double uncert[nProc]={};
  for(int i=1;i<nProc;i++){
    if(i<nMatrix)
      syst[i] = rate[i][ring][sect]/rate[0][ring][sect] * sigma[i][ring][sect];
    else
      syst[i] = rate[i][ring][sect]/rate[0][ring][sect] * fabs(A[i][ring][sect]) * uncert[i] ;
    cout<<"syst "<<procNmComb[i]<<"\t"<<syst[i]<<"\t"<<syst[i]/33<<endl;
  }

}

void printAll(){

  const double polarization = 0.8;
  const double beamDays = 235 + 95 + 14;
  const double days2seconds = 24*60*60;
  const char *sector[3] = {"C", "T", "O"};

  cout<<endl<<",";
  for(int k=0;k<nProc;k++)
    cout<<","<<procNmComb[k]<<",,,,";
  cout<<",,,"<<endl;

  cout<<"R,S";
  for(int k=0;k<nProc;k++)
    cout<<",A [ppb],f [% of rate],f*A [% of Am],dA [ppb],dA/A [% of A]";
  cout<<",Am [ppb],dAm [ppb],dAm/Am [%]"<<endl;

  for(int i=0;i<6;i++)
    for(int j=0;j<3;j++){
      double rateTot(0), fAsum(0);
      for(int k=0;k<nProc;k++){
        rateTot += rate[k][i][j];
        fAsum += rate[k][i][j] * A[k][i][j];
      }
      double sigmaAm = 1/sqrt(rateTot);

      double Ami2(0);
      cout<<i+1<<","<<sector[j];
      for(int k=0;k<nProc;k++){
        Ami2 += rate[k][i][j]/rateTot * A[k][i][j];
        cout<<","<<A[k][i][j]<<","<<rate[k][i][j]/rateTot*100.<<","<<rate[k][i][j]*A[k][i][j]/fAsum*100.
          <<","<<sigma[k][i][j]<<","<<(fabs(A[k][i][j])>1e-9?sigma[k][i][j]/fabs(A[k][i][j])*100.:0);
      }
      double dAmi2 = sigmaAm / ( 0.8 * sqrt(beamDays * days2seconds) ) * 1e9;
      cout<<","<<Ami2<<","<<dAmi2<<","<<dAmi2/fabs(Ami2)*100.<<endl;

    }
}

void readSim(string fnm,int proc,int addProc){

  if(procNm[proc+addProc] == "neutralBknd"){
    for(int i=0;i<6;i++)
      for(int j=0;j<3;j++){
        rate[proc][i][j] = neutralBkgndRate[i]/3*neutralBkgndFactor;
        A[proc][i][j] = 0;
      }
    return;
  }
  
  size_t inel = procNm[proc+addProc].find("epInelasticW");
  size_t Wbin = procNm[proc+addProc].size();
  if(verbose) cout<<"reading "<<fnm<<"\t"<<proc+addProc<<endl;
  TFile *fin=TFile::Open(fnm.c_str(),"READ");
  string hName="hRate";
  if(inel != string::npos)
    hName = Form("hRate_%s",procNm[proc+addProc].substr(Wbin-2,2).c_str());

  TH1D *hRate=(TH1D*)fin->Get(hName.c_str());

  const double rateFactor = 65./85;
  double gfFactor = 1;
  for(int i=0;i<6;i++)
    for(int j=0;j<3;j++){
      hName=Form("hAsym_R%d_S%d",i+1,j);
      if(inel != string::npos)
        hName=Form("hAsym_%s_R%d_S%d",procNm[proc+addProc].substr(Wbin-2,2).c_str(),i+1,j);

      if(verbose)
	cout<<"\tR/S\t"<<i<<"\t"<<j<<"\t"<<hName<<endl;

      double thisA;
      TH1D *hA=(TH1D*)fin->Get(hName.c_str());
      if(procNm[proc+addProc] == "piminus")
        thisA = -1660;
      else
        thisA = hA->GetMean()*gfFactor;
      
      double thisRate = hRate->GetBinContent(i*3+j+1) * rateFactor;
      if(verbose)
	cout<<"\tR/S\t"<<i<<"\t"<<j<<"\t"<<thisA<<"\t"<<thisRate<<endl;
      
      A[proc][i][j] = (rate[proc][i][j]*A[proc][i][j] + thisRate*thisA) / (rate[proc][i][j] + thisRate);
      rate[proc][i][j] += thisRate;
    }

  fin->Close();
  if(verbose)
    cin.ignore();
}


