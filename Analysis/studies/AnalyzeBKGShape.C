#include <string>
#include <vector>
#include <TString.h>
#include <TStyle.h>
#include <iostream>


using namespace std;



void AnalyzeBKGShape(TString datacards_folder, TString folder, TString folder_loose, TString sample){
	gStyle->SetOptStat(000000);
	
	TH1::SetDefaultSumw2();
	
	TFile *f= TFile::Open("/Users/Tita/Desktop/analysis_HH_bbTauTau/src/HHbbTauTau/data/"+datacards_folder+"/htt_tt.inputs-Hhh-8TeV_m_ttbb_kinfit_KinFitConvergedWithMassWindow.root");
	

	
	TH1D *histo = (TH1D*)f->Get((folder+"/"+sample));
	TH1D *histoBLoose = (TH1D*)f->Get((folder_loose+"/"+sample));
	// TH1D *histoBLooseTauLoose = (TH1D*)fLooseBLooseIso->Get((folder+"/"+sample));
	
		
	TCanvas * c1 = new TCanvas("c1","c1", 800,800);

	histo->SetLineColor(kBlue);
	histo->SetMarkerColor(kBlue);
	histo->SetMarkerStyle(20);
	

	histoBLoose->SetLineColor(kRed);
	histoBLoose->SetMarkerColor(kRed);
	histoBLoose->SetMarkerStyle(20);

	
	
	// histoBLooseTauLoose->SetLineColor(kGreen);
// 	histoBLooseTauLoose->SetMarkerColor(kGreen);
// 	histoBLooseTauLoose->SetMarkerStyle(20);
	
    double error_histo = 0;
    double integral_histo = histo->IntegralAndError(1, histo->GetNbinsX(),error_histo);
    cout<<"--------------------------------------------------------------"<<sample<<"--------------------------------------------------------------"<<endl;
    cout<<"*************"<<folder<<"  "<<sample<<"*************"<<endl;
    cout<<" MC Shape "<<integral_histo<<"   ----  "<<error_histo<<endl;
    cout<<" MC Shape Loose B "<<histoBLoose->Integral()<<endl;
//     cout<<" MC Shape Loose B + Relax Tau Iso "<<histoBLooseTauLoose->Integral()<<endl;
    
    
//     histoBLooseTauLoose->Scale(1/(histoBLooseTauLoose->Integral()));
	TH1D* difference = (TH1D*)histoBLoose->Clone("difference");
	difference->Add(histo,-1);
	difference->Scale(1/difference->Integral());
    histoBLoose->Scale(1/(histoBLoose->Integral()));
    histo->Scale(1/(histo->Integral()));
	
	
    
    histo->Draw("");
    histoBLoose->Draw("same");
//     histoBLooseTauLoose->Draw("same");
    
    

    
    cout<<"--------------------------------------------------------------KS   "<<sample<<"--------------------------------------------------------------"<<endl;
    cout<<"*************"<<folder<<"  "<<sample<<"*************"<<endl;
    cout<<" Kolmogorov Test "<<histo->KolmogorovTest(difference,"")<<endl;
	histo->SetTitle("M_{H}");
	histo->GetXaxis()->SetTitle("M_{H} [GeV]");
	histo->GetYaxis()->SetTitleOffset(1.5);
	histo->GetYaxis()->SetTitle("Normalized Events");
	
	
	
		
	
	TLegend* legend = new TLegend(0.6, 0.65, 0.99, 0.9);
	legend->SetFillColor(0);
	legend->SetTextSize(0.03);
	legend->SetEntrySeparation(0.05);
	legend->AddEntry(histo, " Medium CSV, TauIso<1 ");
	legend->AddEntry(histoBLoose, " Loose CSV, TauIso<1 ");
	//legend->AddEntry(histoBLooseTauLoose, " No Btagging, TauIso<1 ");
	
	legend->Draw();
	
	c1->SaveAs("./plots_m_kinFit/"+sample+"_"+folder+".eps");
}


void GetShape(){


	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet1tag", "tauTau_2jetloose1tag", "ZLL");
	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet2tag", "tauTau_2jetloose2tag", "ZLL");
	
	
	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet1tag", "tauTau_2jetloose1tag", "QCD");
	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet2tag", "tauTau_2jetloose2tag", "QCD");
	
	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet1tag", "tauTau_2jetloose1tag", "VV");
	AnalyzeBKGShape("datacards_2014_12_18_tautau_Medium","tauTau_2jet2tag", "tauTau_2jetloose2tag", "VV");
	
	

		
}
