void Print_PU(void){
	
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(kFALSE);

  
  TFile *f1 = new TFile("reweight_LP9Sep_NEW.root");
  TFile *f2 = new TFile("/cmshome/grippo/CMSSW_4_2_6/src/VHTauTau/TreeMaker/test/codice/batch_Analysis/5_10_2011/histo_file_WH.root");
  TFile *f3 = new TFile("/cmshome/grippo/CMSSW_4_2_6/src/VHTauTau/TreeMaker/test/codice/batch_Analysis/5_10_2011/histo_QCDmu.root");
  TFile *f4 = new TFile("/cmshome/grippo/CMSSW_4_2_6/src/VHTauTau/TreeMaker/test/codice/batch_Analysis/5_10_2011/histo_DYjets.root");
  
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,2);
  gStyle->SetOptStat(1111);

  f1->cd();
  c1->cd(1);  
  TH1F* h1 = (TH1F*)f1->Get("hNPU");
  h1->SetTitle("hNPU - WH115");
  h1->Draw();
  c1->cd(2);
  TH1F* h2 = (TH1F*)f1->Get("pileup");
  h2->SetTitle("pileup - data");
  h2->Draw();  
  c1->cd(3);
  TH1F* h3 = (TH1F*)f1->Get("plot_data_div_MC");
  h3->SetTitle("plot_data_div_MC");
  h3->Draw();

  c1->Print("reweighting.jpg");

  TCanvas *c2 = new TCanvas();
  c2->Divide(2,2);
  gStyle->SetOptStat(1111);

  f2->cd();
  c2->cd(1);  
  TH1F* h11 = (TH1F*)f2->Get("hNPU");
  h11->SetTitle("hNPU - WH115");
  h11->Draw();
  c2->cd(2);
  TH1F* h21 = (TH1F*)f2->Get("hnVertex");
  h21->SetTitle("hnVertex - WH115");
  h21->Draw();  
  c2->cd(3);
  TH1F* h31 = (TH1F*)f2->Get("hnVertex_PU");
  h31->SetTitle("hnVertex_PU - WH115");
  h31->Draw();

  c2->Print("WH_115.jpg");

  TCanvas *c3 = new TCanvas();
  c3->Divide(2,2);
  gStyle->SetOptStat(1111);

  f3->cd();
  c3->cd(1);  
  TH1F* h12 = (TH1F*)f3->Get("hNPU");
  h12->SetTitle("hNPU - QCDmu");
  h12->Draw();
  c3->cd(2);
  TH1F* h22 = (TH1F*)f3->Get("hnVertex");
  h22->SetTitle("hnVertex - QCDmu");
  h22->Draw();  
  c3->cd(3);
  TH1F* h32 = (TH1F*)f3->Get("hnVertex_PU");
  h32->SetTitle("hnVertex_PU - QCDmu");
  h32->Draw();

  c3->Print("QCDmu.jpg");

  TCanvas *c4 = new TCanvas();
  c4->Divide(2,2);
  gStyle->SetOptStat(1111);

  f4->cd();
  c4->cd(1);  
  TH1F* h13 = (TH1F*)f4->Get("hNPU");
  h13->SetTitle("hNPU - DYJets");
  h13->Draw();
  c4->cd(2);
  TH1F* h23 = (TH1F*)f4->Get("hnVertex");
  h23->SetTitle("hnVertex - DYJets");
  h23->Draw();  
  c4->cd(3);
  TH1F* h33 = (TH1F*)f4->Get("hnVertex_PU");
  h33->SetTitle("hnVertex_PU - DYJets");
  h33->Draw();

  c4->Print("DYJets.jpg");

  
}
