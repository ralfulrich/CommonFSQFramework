#include <TStyle.h>
#include <tdrstyle.C>

//KAterina
double GetKaterina2013(int isec, int imod) 
 {

  double gain[16][14] = {{0.7510, 0.8700, 2.7370, 0.0000, 0.3630, 0.6430, 0.0000, 0.3100, 0.2120, 0.2740, 0.3030, 0.1690, 0.2650, 0.1550}, 
                    {0.6190, 0.6160, 1.8130, 0.8690, 0.1820, 0.6280, 0.0000, 0.5070, 0.1680, 0.2910, 0.3380, 0.1460, 0.2490, 0.1250}, 
                    {1.0700, 0.6510, 1.4250, 0.7660, 0.3040, 0.1930, 8.2170, 13.2900, 0.4650, 0.2350, 0.0000, 0.2950, 0.3430, 0.3510}, 
                    {0.5310, 0.3300, 0.8910, 0.8260, 0.1170, 0.3300, 0.0000, 0.0000, 0.0000, 0.6390, 0.0000, 0.3170, 0.0000, 0.4580}, 
                    {0.6120, 0.0000, 1.3410, 0.7020, 0.1560, 0.5690, 0.8360, 0.0000, 0.0000, 0.5230, 0.2360, 0.3290, 0.3990, 0.3420}, 
                    {1.3130, 0.4870, 1.4000, 0.6320, 0.1990, 0.7950, 1.2090, 0.0000, 0.5100, 0.7060, 0.2330, 0.2800, 0.4830, 0.4410}, 
                    {0.4160, 0.2820, 1.0430, 0.3130, 0.1140, 0.0860, 250.6690, 0.1950, 0.4200, 6.9160, 3.4790, 1.5110, 4.8590, 3.5340}, 
                    {0.3420, 0.2950, 1.1980, 1.4030, 0.2130, 1.0730, 0.0000, 0.2060, 0.6350, 27.2690, 9.4210, 3.3400, 3.4880, 1.0100}, 
                    {0.3030, 0.8460, 1.4120, 1.0000, 0.2180, 0.8830, 0.0000, 0.1320, 0.1950, 0.2490, 0.2250, 0.2270, 0.2990, 0.2780}, 
                    {0.9040, 1.4030, 2.6580, 1.1900, 0.2350, 1.5570, 0.0000, 0.3160, 0.1990, 0.3100, 0.1790, 0.2510, 0.2510, 0.2520}, 
                    {1.0160, 0.9930, 1.6950, 0.8870, 0.2850, 0.6230, 0.0000, 10.0790, 0.3730, 0.2440, 9.6350, 0.5240, 0.6990, 0.3790}, 
                    {1.1690, 1.1300, 2.1400, 1.3920, 0.2630, 1.2470, 0.0000, 0.0000, 0.5670, 0.3030, 99.3510, 0.3510, 0.1980, 0.3560}, 
                    {0.9160, 1.2700, 1.6430, 0.8070, 0.2310, 2.3020, 0.0000, 0.0000, 0.3230, 0.2910, 0.0000, 0.3430, 0.1280, 0.3080}, 
                    {0.6010, 0.9840, 2.1400, 0.8210, 0.1770, 1.0970, 0.0000, 0.0000, 0.2030, 0.2920, 16.6350, 0.3020, 0.3510, 0.3680}, 
                    {0.7590, 1.3650, 2.9620, 1.1740, 0.3800, 2.3370, 0.0000, 517.2540, 0.2690, 0.0000, 0.1940, 0.2740, 0.2800, 0.4100}, 
                    {0.7420, 0.9720, 2.4600, 0.9240, 0.2200, 0.1630, 3.9070, 0.1970, 0.2700, 0.2580, 0.1510, 0.1340, 0.2790, 0.2620}}; 

  return gain[isec][imod];

 }

bool IsBasChannel(int isec, int imod)
{
  bool bad = false;

  // int badChannelsSecMod[6][2] = {{2,10},{3,11},{12,12},{5,4},{3,8},{7,5}};//2015
  int badChannelsSecMod[8][2] = {{2,10},{3,11},{12,12},{5,4},{3,8},{7,5},{8,8},{4,11}} ;//2013
  for(int ipair=0; ipair<8; ipair++) {
    if( isec+1 == badChannelsSecMod[ipair][0] && 
        imod+1 == badChannelsSecMod[ipair][1])
      bad = true;
  }

  return bad;
}

void GetMuonSignalHist(TH1F* h[16][14],TFile* file)
{
  char name[254];

  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      sprintf(name,"data_PPMinBias_Run2013/MuonSignalSecCh_mod_%d_sec_%d",imod+1,isec+1);

      h[isec][imod] = (TH1F*)file->Get(name);
    }
  }
}

void GetMeanAndSigma(TH1F* h[16][14], 
                     double mean[16][14], 
                     double RMS[16][14], 
                     double Nent[16][14],
                     bool doEMcorrection = true)
{
  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      mean[isec][imod] = h[isec][imod]->GetMean();
      RMS[isec][imod] = h[isec][imod]->GetRMS();
      Nent[isec][imod] = h[isec][imod]->GetEntries();

      if(doEMcorrection && imod < 2) {
        mean[isec][imod] *= 2;
        RMS[isec][imod] *= 2;
      }
    }
  }
}

void PrintHistMap(TH2* h)
{
  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      double intval = h->GetBinContent(imod+1,isec+1);
      Printf("sec %d mod %d value %f",isec+1,imod+1,intval);
    }
  }
}

void WriteHistMapIntoFile(TH2* h,const char* fname)
{
  FILE* pf = fopen(fname,"w");

  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      double intval = h->GetBinContent(imod+1,isec+1);
      fprintf(pf,"sec %d mod %d value %f\n",isec+1,imod+1,intval);
    }
  }

  fclose(pf);
}

void draw_tool_2013()
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptStat("");
  gStyle->SetPaintTextFormat("4.2f");
  TGaxis::SetMaxDigits(3);
  setTDRStyle();
  tdrStyle->SetPadRightMargin(0.04);
  tdrStyle->SetPadLeftMargin(0.165);
  tdrStyle->SetOptFit(0);


  TFile * file = TFile::Open("d20Aug/plotsMuonselectioncuts_2_0000.root");

  TH1F* hMuonSignal[16][14];
  GetMuonSignalHist(hMuonSignal,file);

  TH1F* hMuonDist[16][14];
  

  double mean[16][14];
  double RMS[16][14];
  double Nent[16][14];
  double Kat2013[16][14];
  double sigma_Kat2013[16][14];


  GetMeanAndSigma(hMuonSignal,mean,RMS,Nent);//if you want to apply the doEMcorrection add false or true

  TH2F* hMuonSignalMap = new TH2F("hMuonSignalMap","hMuonSignalMap",
                                  14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonRMSMap = new TH2F("hMuonRMSMap","hMuonRMSMap",
                                14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonMultMap = new TH2F("hMuonMultMap","hMuonMultMap",
                                14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonSigmaMap = new TH2F("hMuonSigmaMap","hMuonSigmaMap",
                                14,0.5,14.5,16,0.5,16.5);

  TH2F* hMuonMeanMap2013 = new TH2F("hMuonMeanMap2013","hMuonMeanMap2013",
                                   14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonSigmaMap2013 = new TH2F("hMuonSigmaMap2013","hMuonSigmaMap2013",
                                   14,0.5,14.5,16,0.5,16.5);

  TH2F* hCalibMap2013 = new TH2F("hCalibMap2013","hCalibMap2013",
                                   14,0.5,14.5,16,0.5,16.5);
  TH2F* hSigmaCalibMap2013 = new TH2F("hSigmaCalibMap2013","hSigmaCalibMap2013",
                                   14,0.5,14.5,16,0.5,16.5);
  TH2F* hRelativeSigmaCalibMap2013 = new TH2F("hRelativeSigmaCalibMap2013","hRelativeSigmaCalibMap2013",
                                   14,0.5,14.5,16,0.5,16.5);
  TH1F* hRelativeSigmaDist = new TH1F("hRelativeSigmaDist","hRelativeSigmaDist",20,0,2);
  
  TH1F* hPullKatMuon2013 = new TH1F("hPullKatMuon2013","hPullKatMuon2013",30,-5,5);
  


  double refmean_pHV = 0;
  double refsigma_pHV = 0;
  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      hMuonSignalMap->Fill(imod+1,isec+1,mean[isec][imod]);
      hMuonRMSMap->Fill(imod+1,isec+1,RMS[isec][imod]);
      hMuonMultMap->Fill(imod+1,isec+1,Nent[isec][imod]);

      double sigma = 0;
      double newsigma = 0;
      if(Nent[isec][imod] != 0) {
        sigma = RMS[isec][imod]/sqrt(Nent[isec][imod]); // + ToyMC_RMS AND not Divide by N
        
      }

      double  sigma;
     

      hMuonSigmaMap->Fill(imod+1,isec+1,sigma);

      

      double mean_pHV = mean[isec][imod] ;/// inv_trpval;
      double sigma_pHV= mean_pHV * sqrt( (sigma/mean[isec][imod])*(sigma/mean[isec][imod])); 
                                        
      hMuonMeanMap2013->Fill(imod+1,isec+1,mean_pHV);
      hMuonSigmaMap2013->Fill(imod+1,isec+1,sigma_pHV);

      if(isec==8 && imod==3) {
        refmean_pHV = mean_pHV;
        refsigma_pHV = sigma_pHV;
      }
     
      Kat2013[isec][imod]= GetKaterina2013(8,3)/GetKaterina2013(isec, imod);
      sigma_Kat2013[isec][imod] =  Kat2013[isec][imod] *
        sqrt( GetKaterina2013(8,3)*(0.2)*GetKaterina2013(8,3) +
              GetKaterina2013(isec,imod)*(0.2)*GetKaterina2013(isec,imod) );




    }
  }

  if(refmean_pHV==0 || refsigma_pHV==0) {
    Printf("ERROR: Reference channel is 0!!!!!!!!");
    return;
  }


  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      double mean_pHV = hMuonMeanMap2013->GetBinContent(imod+1,isec+1);
      double sigma_pHV = hMuonSigmaMap2013->GetBinContent(imod+1,isec+1);

      if(mean_pHV==0 || sigma_pHV==0) continue;
      if(mean_pHV<0 || sigma_pHV<0) continue;

      // Printf("%d %d %d",isec,imod,int(IsBasChannel(isec,imod)));
      if( IsBasChannel(isec,imod) ) continue;

      double calib_pHV = refmean_pHV/mean_pHV;
      double sigma_calib_pHV = calib_pHV *
        sqrt( (refsigma_pHV/refmean_pHV)*(refsigma_pHV/refmean_pHV) + 
              (sigma_pHV/mean_pHV)*(sigma_pHV/mean_pHV) );

      hCalibMap2013->Fill(imod+1,isec+1,calib_pHV);
      hSigmaCalibMap2013->Fill(imod+1,isec+1,sigma_calib_pHV);
      hRelativeSigmaCalibMap2013->Fill(imod+1,isec+1,sigma_calib_pHV/calib_pHV);
      hRelativeSigmaDist->Fill(sigma_calib_pHV/calib_pHV);

      double pull_Katerina = (calib_pHV- Kat2013[isec][imod])/( sigma_calib_pHV*sigma_calib_pHV + sigma_Kat2013[isec][imod]*sigma_Kat2013[isec][imod] );
      hPullKatMuon2013->Fill(pull_Katerina);
    
    }
  }
 

  
  hCalibMap2013->GetXaxis()->SetTitle("Module");
  hCalibMap2013->GetYaxis()->SetTitle("Sector");
  hCalibMap2013->SetXTitle("Calibration constants Map2013");

  
  hPullKatMuon2013->SetXTitle("Pull Dist.MuonInterCalibration2013-2013 KAt");

  TCanvas * c3 = new TCanvas("c3","Muon Signal Physics HV",1000,500);
  c3->Divide(2,1); 
  c3->cd(1)->SetLogz(); hCalibMap2013->Draw("colz text");
  c3->cd(2); hPullKatMuon2013->Draw();
  hPullKatMuon2013->Fit("gaus");
  hPullKatMuon2013->SetXTitle("Pull Dist.MuonInterCalibration2013-2013 KAt");
  Printf("Mean Pull data-kat = %4.2f, RMS Pull data-kat = %4.2f",hPullKatMuon2013->GetMean(),hPullKatMuon2013->GetRMS());

  hRelativeSigmaCalibMap2013->SetTitle("Relative Err. on CalibConst. (Physics HV);Module;Sector");
  hRelativeSigmaCalibMap2013->GetXaxis()->SetTitle("Module");
  hRelativeSigmaCalibMap2013->GetYaxis()->SetTitle("Sector");
  hRelativeSigmaCalibMap2013->SetXTitle("RelativeSigmaCalibMap2013");
  hRelativeSigmaDist->SetXTitle("RelativeSigma");
  TCanvas * c31 = new TCanvas("c31","Relative Errors on Muon InterCalib. Physics HV",1000,500);
  c31->Divide(2,1);
  c31->cd(1); hRelativeSigmaCalibMap2013->Draw("colz text");
  c31->cd(2); hRelativeSigmaDist->Draw(); 

  WriteHistMapIntoFile(hCalibMap2013,"InterCalibValues2013.txt");
}
