#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeBKGShape(TString folder_1, TString folder_2, TString sample){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
    //emb loose
    TFile *f_loose= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/Limits/auxiliaries/shapes/Italians/htt_tt.inputs-Hhh-8TeV_m_ttbb_kinfit_KinFitConvergedWithMassWindow.root");
    //MC medium
    TFile *f_medium= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/datacards_tautau_ZTT_MC_medium/htt_tt.inputs-Hhh-8TeV_m_ttbb_kinfit_KinFitConvergedWithMassWindow.root");
	

    //ZTT mc loose
    TH1D *histoZTT_medium = (TH1D*)f_medium->Get((folder_1+"/"+sample));
    //ZTT emb loose
    TH1D *histoZTT_loose = (TH1D*)f_loose->Get((folder_2+"/"+sample));
	
		
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

    histoZTT_medium->Scale(1/(histoZTT_medium->Integral()));
    histoZTT_medium->SetLineColor(kBlue);
    histoZTT_medium->SetMarkerColor(kBlue);
    histoZTT_medium->SetMarkerStyle(20);
	
    histoZTT_loose->Scale(1/(histoZTT_loose->Integral()));
    histoZTT_loose->SetLineColor(kRed);
    histoZTT_loose->SetMarkerColor(kRed);
    histoZTT_loose->SetMarkerStyle(20);

    histoZTT_medium->Draw("");
    histoZTT_loose->Draw("same");
    

    cout<<"--------------------------------------------------------------KS   "<<sample<<"--------------------------------------------------------------"<<endl;

//    TH1D* difference = (TH1D*)histoZTT_loose->Clone("difference");
//    difference->Add(histoZTT_medium,-1);
//    difference->Scale(1/difference->Integral());
    cout<<" Kolmogorov Test "<<histoZTT_medium->KolmogorovTest(histoZTT_loose,"")<<endl;
    histoZTT_medium->SetTitle("ZTT from Embedded");
    histoZTT_medium->GetXaxis()->SetTitle("M_{H} [GeV]");
    histoZTT_medium->GetYaxis()->SetTitleOffset(1.5);
    histoZTT_medium->GetYaxis()->SetTitle("N Events");
    histoZTT_medium->SetMaximum(0.7);
	
	
		
	
    TLegend* legend = new TLegend(0.65, 0.65, 0.99, 0.9);
	legend->SetFillColor(0);
    legend->SetTextSize(0.02);
	legend->SetEntrySeparation(0.05);
    legend->AddEntry(histoZTT_medium, " ZTT mc Loose CSV");
    legend->AddEntry(histoZTT_loose, " ZTT emb Loose CSV");
	
	legend->Draw();
	
    c1->SaveAs("./ZTT_MC_vs_Emb"+folder_1+".pdf");
}


void GetShape(){

    AnalyzeBKGShape("tauTau_2jetloose1tag", "tauTau_2jetloose1tag", "ZTT");
    AnalyzeBKGShape("tauTau_2jetloose2tag", "tauTau_2jetloose2tag", "ZTT");
		
}
