#pragma once

#include <TChain.h>
#include <TH1.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

struct FlatTree {
    TChain* chain;

    std::vector<Float_t> *pt_BJets;
    std::vector<Float_t> *eta_BJets;
    std::vector<Float_t> *phi_BJets;
    std::vector<Float_t> *energy_BJets;
    std::vector<Float_t> *csv_BJets;

    Int_t nBjet;
    Float_t mvamet_phi, met, drLT, ptsv, etasv, phisv;
    Float_t pt1, pt2, eta1, eta2, phi1, phi2, e1, e2;
    Float_t byCombinedIsolationDeltaBetaCorrRaw3Hits1,byCombinedIsolationDeltaBetaCorrRaw3Hits2,pfRelIso1,mt1,mt2,weight;
    Int_t   Q1,Q2;
    Bool_t  againstMuonTight2,againstElectronLooseMVA2;


    FlatTree() : chain(new TChain("flatTree")), pt_BJets(0), eta_BJets(0), phi_BJets(0), energy_BJets(0), csv_BJets(0)
    {
        chain->SetBranchAddress( "pt_1", &pt1);
        chain->SetBranchAddress( "pt_2", &pt2);
        chain->SetBranchAddress( "eta_1", &eta1);
        chain->SetBranchAddress( "eta_2", &eta2);
        chain->SetBranchAddress( "phi_1", &phi1);
        chain->SetBranchAddress( "phi_2", &phi2);
        chain->SetBranchAddress( "energy_1", &e1);
        chain->SetBranchAddress( "energy_2", &e2);

        chain->SetBranchAddress( "pt_Bjets", &pt_BJets);
        chain->SetBranchAddress( "eta_Bjets", &eta_BJets);
        chain->SetBranchAddress( "phi_Bjets", &phi_BJets);
        chain->SetBranchAddress( "energy_Bjets", &energy_BJets);
        chain->SetBranchAddress( "csv_Bjets", &csv_BJets);
        chain->SetBranchAddress("njetspt20", &nBjet);
        chain->SetBranchAddress("mvametphi", &mvamet_phi);
        chain->SetBranchAddress("mvamet", &met);
        chain->SetBranchAddress("DeltaR_leptons", &drLT);
        chain->SetBranchAddress("pt_sv", &ptsv);
        chain->SetBranchAddress("eta_sv", &etasv);
        chain->SetBranchAddress("phi_sv", &phisv);

        chain->SetBranchAddress("againstMuonTight_2",&againstMuonTight2);
        chain->SetBranchAddress("againstElectronLooseMVA_2",&againstElectronLooseMVA2);
        chain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",&byCombinedIsolationDeltaBetaCorrRaw3Hits1);
        chain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits2);

        chain->SetBranchAddress("pfRelIso_1",&pfRelIso1);
        chain->SetBranchAddress("q_1",&Q1);
        chain->SetBranchAddress("q_2",&Q2);
        chain->SetBranchAddress("mt_1",&mt1);
        chain->SetBranchAddress("mt_2",&mt2);
        chain->SetBranchAddress("weight",&weight);
    }

    ~FlatTree() { delete chain; }

    Int_t Add(const char* name) { return chain->Add(name); }
    Long64_t GetEntriesFast() const { return chain->GetEntriesFast(); }
    Int_t GetEntry(Long64_t entry) { return chain->GetEntry(entry); }
};

inline TH1F* MakeTH1F(const std::string& name, const std::string& prefix, Int_t nbinsx,Double_t xlow,Double_t xup)
{
    const std::string full_name = name + "_" + prefix;
    return new TH1F(full_name.c_str(), full_name.c_str(), nbinsx, xlow, xup);
}

void TMVAtestAndTraining(TMVA::Factory *factory)
{
    TCut mycuts = "";
    TCut mycutb = "";

    factory->PrepareTrainingAndTestTree(mycuts, mycutb,"SplitMode=Random" );

    cout<<"*******************************************Call BDT***************************************"<<endl;
    // Adaptive Boost
//	factory->BookMethod( TMVA::Types::kBDT, "BDT",
//						"!H:!V:NTrees=850:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );


    factory->BookMethod( TMVA::Types::kBDT, "BDT",
                        "!H:!V:NTrees=850:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );


    cout<<"**************************************	Call Fisher******************************************"<<endl;

    factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                            "!H:!V:NTrees=50:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

    cout<<"********************************************************************************"<<endl;

    cout<<"**************************************	Call Decorrelation + AdaBoost******************************************"<<endl;
    // Decorrelation + Adaptive Boost
    factory->BookMethod( TMVA::Types::kBDT, "BDTD","!H:!V:NTrees=400:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );
    cout<<"********************************************************************************"<<endl;


    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------
}

