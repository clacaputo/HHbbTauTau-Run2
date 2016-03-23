#include <iostream>
#include <vector>

#include "TFile.h"
#include "TCut.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TEfficiency.h"

void BTagEff (const TString outputFile = "bTagEff.root"){
	
	TH1::SetDefaultSumw2(kTRUE);
	//TCanvas *c1 = new TCanvas ("c1","Canvas",800,400);
	//c1->Divide(2,1);
	
	TFile *fOut = new TFile(outputFile,"RECREATE");
	
	TFile *file = new TFile("data/76x_v2/TTBar_v2.root");
	
	TTree *tree = (TTree*) file->Get("sync");
	
	float PtBins[]  =  {20, 30, 40, 50, 70, 100, 150, 200, 300, 640} ;
    float EtaBins[] =  {0, 0.6, 1.2, 2.1, 2.4} ;
    int nPtBins  = sizeof(PtBins)/sizeof(float) - 1;
    int nEtaBins = sizeof(EtaBins)/sizeof(float) - 1;
	
	TH2F *h1_b = new TH2F("h1_b","NUM",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *h2_b = new TH2F("h2_b","DENUM",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *h1_c = new TH2F("h1_c","NUM",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *h2_c = new TH2F("h2_c","DENUM",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *h1_l = new TH2F("h1_l","NUM",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *h2_l = new TH2F("h2_l","DENUM",nPtBins,PtBins,nEtaBins,EtaBins);

	/* Muon Cuts */
	const muCuts         = TCut("pt_1 > 19 && TMath::Abs(eta_1) < 2.1 && iso_1 < 0.1");
	
	/* Tau Cuts */
	const tauCuts        = TCut("pt_2 > 20 && TMath::Abs(eta_2) < 2.3");
	const tauIsoCuts     = TCut("iso_2 > 0.2");
	const tauAgainstCuts = TCut("againstMuonTight3_2 && againstElectronVLooseMVA6_2");
	
    /* Common Cuts */
	const vetoCuts       = TCut("dilepton_veto");
	const chargeCut      = TCut("q_1*q_2 == -1");
	const csvCutMedium   = TCut("csv_jets > 0.80");
	const csvCutLoose    = TCut("csv_jets > 0.460");
	
	const lightCuts       = TCut(muCuts+tauCuts+tauAgainstCuts);
	const baselineCuts    = TCut(muCuts+tauCuts+tauIsoCuts+vetoCuts+chargeCut);
	const OS_AntiIso_cuts = TCut(muCuts+tauCuts+!tauIsoCuts+vetoCuts+chargeCut);

	const TCut cut1_b = TCut("TMath::Abs(partonFlavour_jets) == 5 "+lightCuts);
	const TCut cut2_b = TCut("TMath::Abs(partonFlavour_jets) == 5"+csvCutLoose+lightCuts);
	const TCut cut1_c = TCut("TMath::Abs(partonFlavour_jets) == 4 "+lightCuts);
	const TCut cut2_c = TCut("TMath::Abs(partonFlavour_jets) == 4"+csvCutLoose+lightCuts);
	const TCut cut1_l = TCut("TMath::Abs(partonFlavour_jets) != 5 && TMath::Abs(partonFlavour_jets) != 4 "+lightCuts);
	const TCut cut2_l = TCut("TMath::Abs(partonFlavour_jets) != 5 && TMath::Abs(partonFlavour_jets) != 4"+csvCutLoose+lightCuts);
	//tree->Draw("abs(eta_jets):pt_jets>>h1(10,0,500,4,0,2.4)",cut1);
	//tree->Draw("abs(eta_jets):pt_jets>>h2(10,0,500,4,0,2.4)",cut2);
	tree->Project("h1_b","abs(eta_jets):pt_jets",cut1_b);
	tree->Project("h2_b","abs(eta_jets):pt_jets",cut2_b);
	tree->Project("h1_c","abs(eta_jets):pt_jets",cut1_c);
	tree->Project("h2_c","abs(eta_jets):pt_jets",cut2_c);
	tree->Project("h1_l","abs(eta_jets):pt_jets",cut1_l);
	tree->Project("h2_l","abs(eta_jets):pt_jets",cut2_l);
	
	
//	TH2F *h1 = (TH2F*)gDirectory->Get("h1");
//	TH2F *h2 = (TH2F*)gDirectory->Get("h2");
	
// 	TH2F *all    = h1->Clone();
// 	TH2F *passed = h2->Clone();
// 	
// 	h2->Divide(h1);
	
	//h2->Draw("colz");

	TH2F *eff_b = new TH2F("eff_b","Efficiency BTag - b",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *eff_c = new TH2F("eff_c","Efficiency BTag - c",nPtBins,PtBins,nEtaBins,EtaBins);
	TH2F *eff_l = new TH2F("eff_l","Efficiency BTag - uds",nPtBins,PtBins,nEtaBins,EtaBins);
	TEfficiency *pEff_b,*pEff_c,*pEff_l;
// 	
 	if(TEfficiency::CheckConsistency((TH1)h2_b, (TH1)h1_b ))
             {
                pEff_b = new TEfficiency((TH1)h2_b, (TH1)h1_b );
                std::cout<<"b-Jets Confidence Intervall = "<<pEff_b->GetConfidenceLevel()<<std::endl;
                
                for (int gBin = 0; gBin < h2_b->GetSize(); gBin++)
                {
                    eff_b->SetBinContent(gBin, pEff_b->GetEfficiency(gBin));
                    
                    std::cout<<"Eff bin "<<gBin<<" = "<<pEff_b->GetEfficiency(gBin)<<" + "<<pEff_b->GetEfficiencyErrorUp(gBin)<<" - "<<pEff_b->GetEfficiencyErrorLow(gBin)<<std::endl;
                }
                //delete pEff_b; // ? could crash, stupid ROOT...
                pEff_b->Draw("AP");
             }
   if(TEfficiency::CheckConsistency((TH1)h2_c, (TH1)h1_c ))
             {
                pEff_c = new TEfficiency((TH1)h2_c, (TH1)h1_c );
                std::cout<<"c-Jets Confidence Intervall = "<<pEff_c->GetConfidenceLevel()<<std::endl;
                
                for (int gBin = 0; gBin < h2_c->GetSize(); gBin++)
                {
                    eff_c->SetBinContent(gBin, pEff_c->GetEfficiency(gBin));
                }
                //delete pEff_b; // ? could crash, stupid ROOT...
             }
    if(TEfficiency::CheckConsistency((TH1)h2_l, (TH1)h1_l ))
             {
                pEff_l = new TEfficiency((TH1)h2_l, (TH1)h1_l );
                
                std::cout<<"l-Jets Confidence Intervall = "<<pEff_l->GetConfidenceLevel()<<std::endl;
                
                for (int gBin = 0; gBin < h2_l->GetSize(); gBin++)
                {
                    eff_l->SetBinContent(gBin, pEff_l->GetEfficiency(gBin));

                }
                //delete pEff_b; // ? could crash, stupid ROOT...
             }         
// 	TCanvas *c1 = new TCanvas("c1","Eff",1200,400);
// 	c1->Divide(3,1);
// 	//pEff_b->Draw("colZTEXT");
// 	c1->cd(1);   c1->SetLogx(); eff_b->Draw("colZTEXT");
// 	c1->cd(2);   c1->SetLogx(); eff_c->Draw("colZTEXT");
// 	c1->cd(3);   c1->SetLogx(); eff_l->Draw("colZTEXT");
	
	
// 	std::cout<<"eff_b:  "<<eff_b->GetBinContent(45,1.5)<<"\n";
// 	std::cout<<"eff_b:  "<<eff_b->GetBinContent(250,1.5)<<"\n";
// 	std::cout<<"eff_b:  "<<eff_b->GetBinContent(21,1.5)<<"\n";
	
	fOut->cd();
	
	eff_b->Write(); eff_c->Write(); eff_l->Write();
	
	pEff_b->Write(); pEff_c->Write(); pEff_l->Write();
	
	fOut->Close();
}