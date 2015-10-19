#include <TStyle.h>
#include <tdrstyle.C>



double GetTrasposeValues(int isec, int imod)
{
  double gain[16][14] = {{17.8879, 13.2482, 14.0743, 13.3278, 2.12091, 2.12543, 3.21527, 2.98481, 2.11992, 2.06513, 2.10157, 4.18052, 2.10157, 2.10157},
                         {15.5813, 14.3688, 13.9962, 15.8174, 2.14671, 2.16158, 2.94851, 3.03287, 2.07367, 2.09709, 2.02886, 5.6682, 2.10157, 2.10157},
                         {17.1485, 14.2758, 13.6352, 13.7265, 2.19081, 2.16161, 3.12191, 3.12444, 2.11671, 2.14295, 1.68227, 2.62755, 2.90551, 2.25366}, 
                         {25.4937, 15.0353, 13.5396, 14.5396, 2.00951, 2.01052, 2.94567, 3.10935, 2.10436, 2.10581, 1.9966, 2.14238, 1.67304, 2.10112}, 
                         {22.6664, 14.0919, 14.2022, 13.6725, 2.14257, 2.13443, 3.07885, 3.15034, 2.12721, 2.111, 2.12661, 2.04572, 2.1345, 2.23077}, 
                         {13.8189, 13.0366, 12.7392, 14.4154, 2.14965, 2.13784, 2.99889, 3.00984, 2.09463, 2.09944, 2.10157, 0.507937, 2.24675, 2.14433},
                         {14.1375, 14.3402, 12.6179, 14.1134, 2.14431, 2.11088, 2.99728, 3.23413, 2.10955, 2.10675, 2.0835, 6.92593, 1.99355, 2.31884}, 
                         {13.8083, 13.8735, 13.0352, 13.8626, 2.00958, 2.13991, 3.18717, 3.15437, 2.09151, 2.11572, 2.21586, 2.13531, 2.51007, 2.03081},
                         {15.0498, 14.8564, 12.6001, 13.3532, 2.00998, 2.16025, 3.17426, 3.14376, 2.0985, 2.09938, 2.12139, 2.10157, 2.10157, 1.72961},
                         {13.4029, 14.6153, 13.4451, 13.2316, 2.12537, 2.11262, 3.21792, 3.06312, 2.16952, 2.08848, 2.10157, 1.15108, 2.10157, 1.96429}, 
                         {13.8487, 14.6547, 13.9845, 14.004, 2.16241, 1.98043, 3.01777, 3.02263, 2.12504, 2.16387, 1.9798, 2.0, 2.21538, 2.10157}, 
                         {13.5454, 14.1516, 11.5071, 14.1275, 2.08736, 2.11662, 3.24957, 3.10665, 2.13867, 2.04457, 2.22083, 2.04405, 1.98608, 1.98802}, 
                         {13.3621, 14.021, 13.4167, 11.0612, 2.15244, 2.09452, 3.13056, 3.04276, 2.04811, 2.10842, 2.125, 3.80645, 1.75349, 5.33108},
                         {13.4709, 14.0436, 13.4366, 14.2121, 2.12454, 2.12266, 2.99075, 3.13996, 2.06635, 2.07508, 2.31579, 2.23828, 2.27638, 2.10157},
                         {13.2229, 13.5576, 13.6674, 14.794, 2.00953, 2.02293, 3.12361, 3.20957, 2.08386, 2.10041, 0.850949, 2.31136, 4.57522, 2.21341},
                         {13.4984, 13.0597, 13.3461, 13.7108, 2.1543, 2.14047, 2.89888, 3.05443, 2.07617, 2.12138, 2.25668, 2.10459, 2.3318, 2.07519}};

  return gain[isec][imod];
}

double GetTrasposeErrors(int isec, int imod)
{
  double gain_error[16][14] = {{0.297684, 0.303952, 0.154135, 0.143911, 0.00635945, 0.00696305, 0.013834, 0.0321267, 0.0131778, 0.00933203, 0.050539, 0.202872, 0.050539, 0.050539},
                               {0.273256, 0.160228, 0.139015, 0.313545, 0.00771067, 0.00761597, 0.0624354, 0.0248926, 0.0202406, 0.00548391, 0.00676534, 0.234323, 0.050539, 0.050539},
                               {0.208304, 0.209137, 0.20524, 0.190448, 0.0232386, 0.00471184, 0.0277715, 0.0154876, 0.008732, 0.0134014, 0.030463, 0.0378628, 0.0626422, 0.0364017},
                               {0.542134, 0.393696, 0.132817, 0.242715, 0.0507952, 0.13058, 0.0520935, 0.0290168, 0.0121874, 0.0140919, 0.0130064, 0.0335709, 0.026466, 0.0332873},
                               {0.27, 0.160646, 0.16222, 0.195976, 0.00751524, 0.0067034, 0.0115408, 0.0185129, 0.00929156, 0.0129714, 0.0114801, 0.00506028, 1.0658, 0.606506},
                               {0.154156, 0.153823, 0.119356, 0.202521, 0.00716295, 0.0070185, 0.0120852, 0.0464483, 0.00874689, 0.0116634, 0.050539, 0.151539, 0.377385, 0.410735},
                               {0.184748, 0.160822, 0.393225, 0.178964, 0.00593741, 0.00470421, 0.017487, 0.0232038, 0.00917938, 0.0107265, 0.0160861, 0.356823, 0.0939621, 0.0386872},
                               {0.157523, 0.175195, 0.16801, 0.105009, 0.0290363, 0.00498382, 0.0212775, 0.022056, 0.0157854, 0.0104235, 0.194334, 0.206459, 0.261539, 0.0521384},
                               {0.151579, 0.143932, 0.135856, 0.171023, 0.253775, 0.0120092, 0.0237112, 0.0268974, 0.0142285, 0.0211197, 0.067731, 0.050539, 0.050539, 0.336312},
                               {0.163748, 0.161818, 0.181251, 0.253084, 0.00782319, 0.00626064, 0.0318485, 0.0249402, 0.00646081, 0.0111822, 0.050539, 0.189923, 0.050539, 0.653284},
                               {0.238548, 0.188285, 0.12227, 0.172699, 0.00682779, 0.00786426, 0.0209686, 0.0204322, 0.0101018, 0.0140415, 0.221228, 0.0962536, 0.365371, 0.050539},
                               {0.166892, 0.243934, 0.206652, 0.140395, 0.0051994, 0.00605701, 0.0256196, 0.0165193, 0.00994767, 0.0158936, 0.132698, 0.125007, 0.374697, 0.385424},
                               {0.115767, 0.173322, 0.126637, 0.12721, 0.00871158, 0.00421985, 0.0232185, 0.0220267, 0.00926294, 0.0144162, 0.0723647, 0.171277, 0.0925839, 0.208289},
                               {0.168528, 0.173265, 0.146743, 0.202867, 0.00571115, 0.00731292, 0.0177893, 0.0232495, 0.00932107, 0.00953935, 0.22871, 0.108933, 0.174756, 0.050539},
                               {0.183391, 0.137414, 0.124799, 0.168779, 0.0485851, 0.00782375, 0.0183694, 0.0220257, 0.00824439, 0.00902181, 0.0328126, 0.138288, 0.112013, 0.173584},
                               {0.257366, 0.167887, 0.272659, 0.184186, 0.00944811, 0.0108569, 0.0476646, 0.0224154, 0.00838952, 0.0083331, 0.13013, 0.0684296, 0.156254, 0.13299}};

  return gain_error[isec][imod];
}

double GetLED1800V(int isec, int imod)
{
  double gain[16][14]= {{75500.0, 1300000.0, 568000.0, 79100.0, 0.0, 362000.0, 330000.0, 512000.0, 2130000.0, 121000.0, 133000.0, 1760000.0, 104000.0, 1560000.0},
                  {72200.0, 530000.0, 344000.0, 0.0, 793000.0, 270000.0, 47400.0, 467000.0, 1140000.0, 0.0, 55200.0, 2460000.0, 611000.0, 525000.0}, 
                  {362000.0, 557000.0, 423000.0, 843000.0, 320000.0, 274000.0, 382000.0, 547000.0, 457000.0, 271000.0, 683000.0, 515000.0, 369000.0, 462000.0}, 
                  {579000.0, 522000.0, 377000.0, 746000.0, 293000.0, 248000.0, 297000.0, 463000.0, 1080000.0, 477000.0, 542000.0, 647000.0, 788000.0, 561000.0},
                  {354000.0, 94100.0, 326000.0, 0.0, 0.0, 307000.0, 486000.0, 456000.0, 325000.0, 0.0, 1570000.0, 393000.0, 365000.0, 638000.0},
                  {644000.0, 0.0, 0.0, 0.0, 310000.0, 282000.0, 226000.0, 486000.0, 463000.0, 420000.0, 0.0, 160000.0, 346000.0, 624000.0}, 
                  {535000.0, 486000.0, 414000.0, 703000.0, 285000.0, 83900.0, 397000.0, 661000.0, 832000.0, 311000.0, 423000.0, 748000.0, 309000.0, 320000.0}, 
                  {645000.0, 509000.0, 406000.0, 501000.0, 335000.0, 279000.0, 483000.0, 676000.0, 757000.0, 278000.0, 503000.0, 647000.0, 374000.0, 857000.0}, 
                  {0.0, 172000.0, 310000.0, 543000.0, 213000.0, 159000.0, 448000.0, 569000.0, 974000.0, 710000.0, 734000.0, 296000.0, 0.0, 403000.0}, 
                  {551000.0, 550000.0, 298000.0, 583000.0, 263000.0, 142000.0, 568000.0, 581000.0, 846000.0, 613000.0, 279000.0, 320000.0, 291000.0, 385000.0}, 
                  {327000.0, 515000.0, 303000.0, 506000.0, 160000.0, 140000.0, 506000.0, 532000.0, 756000.0, 920000.0, 392000.0, 324000.0, 288000.0, 304000.0}, 
                  {356000.0, 214000.0, 359000.0, 624000.0, 214000.0, 247000.0, 524000.0, 48600.0, 717000.0, 732000.0, 533000.0, 464000.0, 856000.0, 664000.0}, 
                  {483000.0, 564000.0, 0.0, 469000.0, 349000.0, 164000.0, 522000.0, 558000.0, 355000.0, 596000.0, 442000.0, 472000.0, 377000.0, 789000.0},
                  {573000.0, 431000.0, 293000.0, 401000.0, 346000.0, 180000.0, 435000.0, 0.0, 473000.0, 556000.0, 440000.0, 573000.0, 453000.0, 0.0},
                  {441000.0, 426000.0, 363000.0, 0.0, 345000.0, 122000.0, 459000.0, 436000.0, 219000.0, 663000.0, 314000.0, 631000.0, 517000.0, 363000.0},
                  {536000.0, 487000.0, 354000.0, 0.0, 375000.0, 212000.0, 274000.0, 483000.0, 504000.0, 806000.0, 844000.0, 825000.0, 506000.0, 552000.0}};

  return gain[isec][imod];
}

double GetLEDError1800V(int isec, int imod)
{
  double gain_error[16][14] = {{0.297684, 0.303952, 0.154135, 0.143911, 0.00635945, 0.00696305, 0.013834, 0.0321267, 0.0131778, 0.00933203, 0.050539, 0.202872, 0.050539, 0.050539}, 
                          {0.273256, 0.160228, 0.139015, 0.313545, 0.00771067, 0.00761597, 0.0624354, 0.0248926, 0.0202406, 0.00548391, 0.00676534, 0.234323, 0.050539, 0.050539}, 
                          {0.208304, 0.209137, 0.20524, 0.190448, 0.0232386, 0.00471184, 0.0277715, 0.0154876, 0.008732, 0.0134014, 0.030463, 0.0378628, 0.0626422, 0.0364017}, 
                          {0.542134, 0.393696, 0.132817, 0.242715, 0.0507952, 0.13058, 0.0520935, 0.0290168, 0.0121874, 0.0140919, 0.0130064, 0.0335709, 0.026466, 0.0332873}, 
                          {0.27, 0.160646, 0.16222, 0.195976, 0.00751524, 0.0067034, 0.0115408, 0.0185129, 0.00929156, 0.0129714, 0.0114801, 0.00506028, 1.0658, 0.606506}, 
                          {0.154156, 0.153823, 0.119356, 0.202521, 0.00716295, 0.0070185, 0.0120852, 0.0464483, 0.00874689, 0.0116634, 0.050539, 0.151539, 0.377385, 0.410735}, 
                          {0.184748, 0.160822, 0.393225, 0.178964, 0.00593741, 0.00470421, 0.017487, 0.0232038, 0.00917938, 0.0107265, 0.0160861, 0.356823, 0.0939621, 0.0386872}, 
                          {0.157523, 0.175195, 0.16801, 0.105009, 0.0290363, 0.00498382, 0.0212775, 0.022056, 0.0157854, 0.0104235, 0.194334, 0.206459, 0.261539, 0.0521384}, 
                          {0.151579, 0.143932, 0.135856, 0.171023, 0.253775, 0.0120092, 0.0237112, 0.0268974, 0.0142285, 0.0211197, 0.067731, 0.050539, 0.050539, 0.336312}, 
                          {0.163748, 0.161818, 0.181251, 0.253084, 0.00782319, 0.00626064, 0.0318485, 0.0249402, 0.00646081, 0.0111822, 0.050539, 0.189923, 0.050539, 0.653284}, 
                          {0.238548, 0.188285, 0.12227, 0.172699, 0.00682779, 0.00786426, 0.0209686, 0.0204322, 0.0101018, 0.0140415, 0.221228, 0.0962536, 0.365371, 0.050539}, 
                          {0.166892, 0.243934, 0.206652, 0.140395, 0.0051994, 0.00605701, 0.0256196, 0.0165193, 0.00994767, 0.0158936, 0.132698, 0.125007, 0.374697, 0.385424}, 
                          {0.115767, 0.173322, 0.126637, 0.12721, 0.00871158, 0.00421985, 0.0232185, 0.0220267, 0.00926294, 0.0144162, 0.0723647, 0.171277, 0.0925839, 0.208289}, 
                          {0.168528, 0.173265, 0.146743, 0.202867, 0.00571115, 0.00731292, 0.0177893, 0.0232495, 0.00932107, 0.00953935, 0.22871, 0.108933, 0.174756, 0.050539}, 
                          {0.183391, 0.137414, 0.124799, 0.168779, 0.0485851, 0.00782375, 0.0183694, 0.0220257, 0.00824439, 0.00902181, 0.0328126, 0.138288, 0.112013, 0.173584}, 
                          {0.257366, 0.167887, 0.272659, 0.184186, 0.00944811, 0.0108569, 0.0476646, 0.0224154, 0.00838952, 0.0083331, 0.13013, 0.0684296, 0.156254, 0.13299}};
  return gain_error[isec][imod];
}




bool IsBasChannel(int isec, int imod)
{
  bool bad = false;

  int badChannelsSecMod[6][2] = {{2,10},{3,11},{12,12},{5,4},{3,8},{7,5}};//2015

  for(int ipair=0; ipair<6; ipair++) {
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
      sprintf(name,"data_MinimumBias_Run2015A/MuonSignalSecCh_mod_%d_sec_%d",imod+1,isec+1);

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

void GetMeanDistributionHist(TH1F* h[16][14], TFile* file) 
{
  char name[254];

  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      sprintf(name,"data_MinimumBias_Run2015A/1DSignalMCCh_mod_%d_sec_%d",imod+1,isec+1);

      h[isec][imod] = (TH1F*)file->Get(name);
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

void draw_tool_2015(bool sigma_with_toymc = true)
{
  
  gStyle->SetOptStat(0);
  gStyle->SetOptStat("");
  gStyle->SetPaintTextFormat("4.2f");
  TGaxis::SetMaxDigits(3);
  setTDRStyle();
  tdrStyle->SetPadRightMargin(0.04);
  tdrStyle->SetPadLeftMargin(0.165);
  tdrStyle->SetOptFit(0);



  TFile * file = TFile::Open("input_for_toyMC.root");
  TFile * mcfile = TFile::Open("d24Sep/output_for_toyMC.root");

  TH1F* hMuonSignal[16][14];
  GetMuonSignalHist(hMuonSignal,file);

  TH1F* hMuonDist[16][14];
  GetMeanDistributionHist(hMuonDist,mcfile);

  double mean[16][14];
  double RMS[16][14];
  double Nent[16][14];

  double led[16][14];
  double sigma_led[16][14];

  GetMeanAndSigma(hMuonSignal,mean,RMS,Nent);//if you want to apply the doEMcorrection add false or true

  TH2F* hMuonSignalMap = new TH2F("hMuonSignalMap","hMuonSignalMap",
                                  14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonRMSMap = new TH2F("hMuonRMSMap","hMuonRMSMap",
                                14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonMultMap = new TH2F("hMuonMultMap","hMuonMultMap",
                                14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonSigmaMap = new TH2F("hMuonSigmaMap","hMuonSigmaMap",
                                14,0.5,14.5,16,0.5,16.5);

  TH2F* hMuonMeanPHVMap = new TH2F("hMuonMeanPHVMap","hMuonMeanPHVMap",
                                   14,0.5,14.5,16,0.5,16.5);
  TH2F* hMuonSigmaPHVMap = new TH2F("hMuonSigmaPHVMap","hMuonSigmaPHVMap",
                                   14,0.5,14.5,16,0.5,16.5);

  TH2F* hCalibPHVMap = new TH2F("hCalibPHVMap","hCalibPHVMap",
                                   14,0.5,14.5,16,0.5,16.5);
  TH2F* hSigmaCalibPHVMap = new TH2F("hSigmaCalibPHVMap","hSigmaCalibPHVMap",
                                   14,0.5,14.5,16,0.5,16.5);
  TH2F* hRelativeSigmaCalibPHVMap = new TH2F("hRelativeSigmaCalibPHVMap","hRelativeSigmaCalibPHVMap",
                                   14,0.5,14.5,16,0.5,16.5);
  TH1F* hRelativeSigmaDist = new TH1F("hRelativeSigmaDist","hRelativeSigmaDist",20,0,2);

  TH1F* hPullLEDMuon = new TH1F("hPullLEDMuon","hPullLEDMuon",30,-3,3);



  TH2F* hToyMC_Mean = new TH2F("hToyMC_Mean","hToyMC_Mean",
                                14,0.5,14.5,16,0.5,16.5);
  TH2F* hToyMC_RMS = new TH2F("hToyMC_RMS","hToyMC_RMS",
                                14,0.5,14.5,16,0.5,16.5);

  TH1F* hPullToyMCData = new TH1F("hPullToyMCData","hPullToyMCData",30,-3,3);


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
        sigma = RMS[isec][imod]/sqrt(Nent[isec][imod]); // + ToyMC_RMS AND not Divide by N zou did not add RMS from ToyMc
        newsigma = hMuonDist[isec][imod]->GetRMS();
      }

      double used_sigma = sigma;
      if(sigma_with_toymc) used_sigma = newsigma;

      hMuonSigmaMap->Fill(imod+1,isec+1,used_sigma);

      double inv_trpval = GetTrasposeValues(isec,imod);
      double inv_trperr = GetTrasposeErrors(isec,imod);

      double mean_pHV = mean[isec][imod] / inv_trpval;
      double sigma_pHV= mean_pHV * sqrt( (used_sigma/mean[isec][imod])*(used_sigma/mean[isec][imod]) + 
                                         (inv_trperr/inv_trpval)*(inv_trperr/inv_trpval) );
      hMuonMeanPHVMap->Fill(imod+1,isec+1,mean_pHV);
      hMuonSigmaPHVMap->Fill(imod+1,isec+1,sigma_pHV);

      if(isec==8 && imod==3) {
        refmean_pHV = mean_pHV;
        refsigma_pHV = sigma_pHV;
      }

      led[isec][imod] =  GetLED1800V(8,3)/GetLED1800V(isec,imod);
      sigma_led[isec][imod] = led[isec][imod] *
        sqrt( GetLEDError1800V(8,3)*GetLEDError1800V(8,3) +
              GetLEDError1800V(isec,imod)*GetLEDError1800V(isec,imod) );


      hToyMC_Mean->Fill(imod+1,isec+1,hMuonDist[isec][imod]->GetMean());
      hToyMC_RMS->Fill(imod+1,isec+1,hMuonDist[isec][imod]->GetRMS());
    }
  }

  if(refmean_pHV==0 || refsigma_pHV==0) {
    Printf("ERROR: Reference channel is 0!!!!!!!!");
    return;
  }


  for(int isec=0; isec<16; isec++) {
    for(int imod=0; imod<14; imod++) {
      double mean_pHV = hMuonMeanPHVMap->GetBinContent(imod+1,isec+1);
      double sigma_pHV = hMuonSigmaPHVMap->GetBinContent(imod+1,isec+1);

      if(mean_pHV==0 || sigma_pHV==0) continue;
      if(mean_pHV<0 || sigma_pHV<0) continue;

      // Printf("%d %d %d",isec,imod,int(IsBasChannel(isec,imod)));
      if( IsBasChannel(isec,imod) ) continue;

      double calib_pHV = refmean_pHV/mean_pHV;
      double sigma_calib_pHV = calib_pHV *
        sqrt( (refsigma_pHV/refmean_pHV)*(refsigma_pHV/refmean_pHV) + 
              (sigma_pHV/mean_pHV)*(sigma_pHV/mean_pHV) );

      hCalibPHVMap->Fill(imod+1,isec+1,calib_pHV);
      hSigmaCalibPHVMap->Fill(imod+1,isec+1,sigma_calib_pHV);
      hRelativeSigmaCalibPHVMap->Fill(imod+1,isec+1,sigma_calib_pHV/calib_pHV);
      hRelativeSigmaDist->Fill(sigma_calib_pHV/calib_pHV);

      double pull_val = (calib_pHV-led[isec][imod]) / ( sigma_calib_pHV*sigma_calib_pHV + sigma_led[isec][imod]*sigma_led[isec][imod] );

      hPullLEDMuon->Fill(pull_val);


      double pullToyMCData = (calib_pHV-hToyMC_Mean->GetBinContent(imod+1,isec+1)) / 
        (sigma_calib_pHV*sigma_calib_pHV + 
          hToyMC_RMS->GetBinContent(imod+1,isec+1)*hToyMC_RMS->GetBinContent(imod+1,isec+1) );

      hPullToyMCData->Fill(pullToyMCData);
    }
  }
 

  
  hCalibPHVMap->GetXaxis()->SetTitle("Module");
  hCalibPHVMap->GetYaxis()->SetTitle("Sector");
  hCalibPHVMap->SetXTitle("Calibration constants PHVMap");

  hSigmaCalibPHVMap->GetXaxis()->SetTitle("Module");
  hSigmaCalibPHVMap->GetYaxis()->SetTitle("Sector");
  hSigmaCalibPHVMap->SetXTitle("SigmaCalibration constants PHVMap");

  TCanvas * c3 = new TCanvas("c3","Muon Signal Physics HV",1000,500);
  c3->Divide(2,1); 
  c3->cd(1)->SetLogz(); hCalibPHVMap->Draw("colz text");
  c3->cd(2)->SetLogz(); hSigmaCalibPHVMap->Draw("colz text");

  hRelativeSigmaCalibPHVMap->SetTitle("Relative Err. on CalibConst. (Physics HV);Module;Sector");
  hRelativeSigmaCalibPHVMap->GetXaxis()->SetTitle("Module");
  hRelativeSigmaCalibPHVMap->GetYaxis()->SetTitle("Sector");
  hRelativeSigmaCalibPHVMap->SetXTitle("RelativeSigmaCalibPHVMap");
  hRelativeSigmaDist->SetXTitle("RelativeSigma");
  TCanvas * c31 = new TCanvas("c31","Relative Errors on Muon InterCalib. Physics HV",1000,500);
  c31->Divide(2,1);
  c31->cd(1); hRelativeSigmaCalibPHVMap->Draw("colz text");
  c31->cd(2); hRelativeSigmaDist->Draw();

  TCanvas * c4 = new TCanvas("c4","Pull Distribution MuonInterCalib. LED");
  c4->cd(); hPullLEDMuon->Draw();
  hPullLEDMuon->Fit("gaus");
  hPullLEDMuon->SetXTitle("Pull Distribution MuonInterCalib. LED");
  Printf("Mean Pull = %4.2f, RMS Pull = %4.2f",hPullLEDMuon->GetMean(),hPullLEDMuon->GetRMS());

  TCanvas * c5 = new TCanvas("c5","ToyMC Mean & RMS",1000,500);
  
  hToyMC_Mean->GetXaxis()->SetTitle("Module");
  hToyMC_Mean->GetYaxis()->SetTitle("Sector");
  hToyMC_Mean->SetXTitle(" Mean ToyMC");

  hToyMC_RMS->GetXaxis()->SetTitle("Module");
  hToyMC_RMS->GetYaxis()->SetTitle("Sector");
  hToyMC_RMS->SetXTitle("RMS ToyMC");
  c5->Divide(2,1);
  c5->cd(1);
  hToyMC_Mean->Draw("colz text");
  c5->cd(2);
  hToyMC_RMS->Draw("colz text");

  TCanvas * c6 = new TCanvas("c6","Pull Distribution ToyMC and Data");
  c6->cd(); hPullToyMCData->Draw();
  hPullToyMCData->Fit("gaus");
  hPullToyMCData->SetXTitle("Pull Distribution MuonInterCalib. ToyMC");
  Printf(" toy-data Mean Pull = %4.2f, RMS Pull = %4.2f",hPullToyMCData->GetMean(),hPullToyMCData->GetRMS());


  WriteHistMapIntoFile(hCalibPHVMap,"InterCalibValues_2015.txt");
}
