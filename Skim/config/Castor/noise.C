//void
//noise(TFile* f)
{
  gROOT->Reset();
  //TH1D* h1 = (TH1D*) f->Get("hottest_tower__BptxMinus");
  //TH1D* h2 = (TH1D*) f->Get("hottest_tower__BptxPlus");
  //TH1D* h3 = (TH1D*) f->Get("hottest_tower__NotBptxOR");

  gStyle->SetOptTitle(0);
  
  TCanvas* canvas = new TCanvas("tower", "tower", 3000, 1000);
  canvas->Divide(3);


  TH1D* h1 = (TH1D*) hottest_tower_BptxMinus;
  TH1D* h2 = (TH1D*) hottest_tower_BptxPlus;
  TH1D* h3 = (TH1D*) hottest_tower_NotBptxOR;
  
  
  h1->Scale(1./h1->Integral("width"));
  h1->SetLineWidth(3);
  h1->SetLineColor(kBlue);

  h2->Scale(1./h2->Integral("width"));
  h2->SetLineWidth(3);
  h2->SetLineColor(kRed);

  h3->Scale(1./h3->Integral("width"));
  h3->SetLineWidth(3);
  h3->SetLineColor(kGreen+2);  

  canvas->cd(1);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.13);
  
  //h1->SetTitle("Full tower energy");
  h1->SetXTitle("RecHit / fC");
  h1->SetYTitle("Probability / BX / fC");
  h1->Draw("hist");
  h2->Draw("hist,same");
  h3->Draw("hist,same");
  gPad->BuildLegend(0.2,0.5,0.9,0.7);
  


  TH1D* h1em = (TH1D*) hottest_tower_em_BptxMinus;
  TH1D* h2em = (TH1D*) hottest_tower_em_BptxPlus;
  TH1D* h3em = (TH1D*) hottest_tower_em_NotBptxOR;
  
  
  h1em->Scale(1./h1em->Integral("width"));
  h1em->SetLineWidth(3);
  h1em->SetLineColor(kBlue);

  h2em->Scale(1./h2em->Integral("width"));
  h2em->SetLineWidth(3);
  h2em->SetLineColor(kRed);

  h3em->Scale(1./h3em->Integral("width"));
  h3em->SetLineWidth(3);
  h3em->SetLineColor(kGreen+2);
  
  canvas->cd(2);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.13);

  //th1em->SetTitle("EM tower energy");
  h1em->SetXTitle("RecHit / fC");
  h1em->SetYTitle("Probability / BX / fC");
  h1em->Draw("hist");
  h2em->Draw("hist,same");
  h3em->Draw("hist,same");
  gPad->BuildLegend(0.2,0.5,0.9,0.7);
  



  TH1D* h1had = (TH1D*) hottest_tower_had_BptxMinus;
  TH1D* h2had = (TH1D*) hottest_tower_had_BptxPlus;
  TH1D* h3had = (TH1D*) hottest_tower_had_NotBptxOR;
  
  
  h1had->Scale(1./h1had->Integral("width"));
  h1had->SetLineWidth(3);
  h1had->SetLineColor(kBlue);

  h2had->Scale(1./h2had->Integral("width"));
  h2had->SetLineWidth(3);
  h2had->SetLineColor(kRed);

  h3had->Scale(1./h3had->Integral("width"));
  h3had->SetLineWidth(3);
  h3had->SetLineColor(kGreen+2);
  

  canvas->cd(3);
  gPad->SetLogy(1);
  gPad->SetLeftMargin(0.13);

  //h1had->SetTitle("Hadronic tower energy");
  h1had->SetXTitle("RecHit / fC");
  h1had->SetYTitle("Probability / BX / fC");
  h1had->Draw("hist");
  h2had->Draw("hist,same");
  h3had->Draw("hist,same");
  gPad->BuildLegend(0.2,0.5,0.9,0.7);
  
  
  
}
