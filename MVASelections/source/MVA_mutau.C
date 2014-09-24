
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"

#include "../include/MVA_common.h"

void ApplySelection(TMVA::Factory *factory, FlatTree* tree, bool is_signal)
{
    Double_t sigWeight = 1.;
    Double_t bkgWeight = 1.;

    std::vector<Double_t> vars( 13 );

    const std::string prefix = is_signal ? "signal" : "bkg";
    TH1F* hTauPt = MakeTH1F("hTauPt",prefix,200,0,200);
    TH1F* hMuonPt = MakeTH1F("hMuonPt",prefix,200,0,200);
    TH1F* hPtb1 = MakeTH1F("hPtb1",prefix,300, 0, 300);
    TH1F* hptb2 = MakeTH1F("hPtb2",prefix,300, 0, 300);
    TH1F* hDRbb = MakeTH1F("hDRbb",prefix,100, 0, 10);
    TH1F* hDPhiBBMET = MakeTH1F("hDPhiBBMET",prefix,100, -3, 3);
    TH1F* hDRll = MakeTH1F("hDRll",prefix,100, 0, 10);
    TH1F* hPtHtt = MakeTH1F("hPtHtt",prefix,600, 0, 600);
    TH1F* hDRHBBTT= MakeTH1F("hDRHBBTT",prefix,100, 0, 10);
    TH1F* hPtHBB= MakeTH1F("hPtHBB",prefix,600, 0, 600);
    TH1F* hDeltaPhi_METTT= MakeTH1F("hDeltaPhi_METTT",prefix,100, -3, 3);
    TH1F* hPtH= MakeTH1F("hPtH",prefix,1000, 0, 1000);
    TH1F* hmT2= MakeTH1F("hmT2",prefix,600, 0, 600);

    for (Long64_t i=0; i<(tree->GetEntriesFast()); i++) {
        tree->GetEntry(i);
        if((tree->againstMuonTight2) && (tree->againstElectronLooseMVA2 ) && (tree->byCombinedIsolationDeltaBetaCorrRaw3Hits2<1.5) && (tree->pfRelIso1<0.1) && (tree->Q1*tree->Q2<0) && (tree->mt1<30)){
            //if(!(nBjet>1)) continue;
            vars[2]=0;
            vars[3]=0;

            std::vector<int> IndexNjetEta2p4, IndexBjetsMedium;
            for (unsigned int k =0; k<tree->nBjet; k++){
                if((TMath::Abs(tree->eta_BJets->at(k)) <2.4   )  ) {IndexNjetEta2p4.push_back(k);}
            }

            if ((IndexNjetEta2p4.size() > 1)) {


                for (unsigned int k =0; k<IndexNjetEta2p4.size() ; k++){
                    //cout<<i<<"pt bjets  "<<(*pt_BJets).at(IndexNjetEta2p4.at(k)) <<" csv  "<<(*csv_BJets).at(IndexNjetEta2p4.at(k))<<endl;
                    if(! ( tree->pt_BJets->at(IndexNjetEta2p4.at(k)) > 20  )) continue;
                    if(!(TMath::Abs(tree->eta_BJets->at(IndexNjetEta2p4.at(k))) <2.4   )  ) continue;
                    if(! (tree->csv_BJets->at(IndexNjetEta2p4.at(k)) > 0.689 )  ) continue;
                    IndexBjetsMedium.push_back(IndexNjetEta2p4.at(k));
                }

                bool BjetsCut =false;
                if(  IndexBjetsMedium.size() >1 )  BjetsCut = true;
                if ((BjetsCut == true)){

                    TLorentzVector b1, b2, BB, MET, l1, l2, TT, H;
                    b1.SetPtEtaPhiE(tree->pt_BJets->at(IndexBjetsMedium.at(0)), tree->eta_BJets->at(IndexBjetsMedium.at(0)), tree->phi_BJets->at(IndexBjetsMedium.at(0)) ,tree->energy_BJets->at(IndexBjetsMedium.at(0)) );
                    b2.SetPtEtaPhiE(tree->pt_BJets->at(IndexBjetsMedium.at(1)), tree->eta_BJets->at(IndexBjetsMedium.at(1)), tree->phi_BJets->at(IndexBjetsMedium.at(1)) ,tree->energy_BJets->at(IndexBjetsMedium.at(1)) );

                    l1.SetPtEtaPhiE( tree->pt1 , tree->eta1, tree->phi1, tree->e1 );
                    l2.SetPtEtaPhiE( tree->pt2 , tree->eta2, tree->phi2, tree->e2 );

                    MET.SetPtEtaPhiE(tree->met, 0., tree->mvamet_phi, tree->met);

                    BB = b1+b2;
                    TT=l2+l1;
                    H = BB+TT;



                    //////////////////////////////
                    vars[0] = tree->pt1;
                    vars[1] = tree->pt2;
                    vars[2] = tree->pt_BJets->at(IndexBjetsMedium.at(0));
                    vars[3] = tree->pt_BJets->at(IndexBjetsMedium.at(1));
                    vars[4] = b1.DeltaR(b2);
                    vars[5] = MET.DeltaPhi(BB);  //cout<<"  Delta phi MET*****************  "<<vars[5]<<endl;
                    vars[6] = tree->drLT;              // cout<<"  Delta leptoni *****************  "<<vars[6]<<endl;
                    vars[7]=  TT.Pt();            //cout<<"  PT SV *****************  "<<vars[7]<<endl;
                    vars[8] = TT.DeltaR(BB);
                    vars[9] = BB.Pt();
                    vars[10] = MET.DeltaPhi(TT);
                    vars[11] = H.Pt();
                    vars[12] = tree->mt2;

                    hMuonPt->Fill(vars[0]);
                    hTauPt->Fill(vars[1]);
                    hPtb1->Fill(vars[2]);
                    hptb2->Fill(vars[3]);
                    hDRbb->Fill(vars[4] );
                    hDPhiBBMET->Fill(vars[5] );
                    hDRll->Fill(vars[6] );
                    hPtHtt->Fill(vars[7] );
                    hDRHBBTT->Fill(vars[8] );
                    hPtHBB->Fill(vars[9] );
                    hDeltaPhi_METTT->Fill(vars[10] );
                    hPtH->Fill(vars[11] );
                    hmT2->Fill(vars[12] );

                    //cout<<" rnd "<<gRandom->Uniform(0, 1)<<endl;
                    double r = gRandom->Uniform(0, 1);


                    if (r<0.5)  {
                        //cout<<"Riempi Training ev="<<i<<endl;
                        if(is_signal)
                            factory->AddSignalTrainingEvent(vars, tree->weight*sigWeight );
                        else
                            factory->AddBackgroundTrainingEvent( vars, tree->weight*bkgWeight );
                    } else {
                        if(is_signal)
                            factory->AddSignalTestEvent(vars, tree->weight*sigWeight );
                        else
                            factory->AddBackgroundTestEvent(vars, tree->weight*bkgWeight );
                    }

                }
            }
        }
    }

    hMuonPt->Write();
    hTauPt->Write();
    hPtb1->Write();
    hptb2->Write();
    hDRbb->Write();
    hDPhiBBMET->Write();
    hDRll->Write();
    hPtHtt->Write();
    hDRHBBTT->Write();
    hPtHBB->Write();
    hDeltaPhi_METTT->Write();
    hPtH->Write();
    hmT2->Write();
}

void MVA_mutau(const TString& filePath)
{
	std::cout << "==> Start TMVAClassification" << std::endl;
	
	
	// --------------------------------------------------------------------------------------------------
	
    TString outfileName( "./out_mutau.root" );
	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

		
    TMVA::Factory *factory = new TMVA::Factory( "TMVA_muTau", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
    
    factory->AddVariable("pt_mu", 'F');
    factory->AddVariable("pt_tau", 'F');
    factory->AddVariable("pt_b1", 'F');
	factory->AddVariable("pt_b2", 'F');
	factory->AddVariable("DR_bb", 'F');
	factory->AddVariable("DPhi_BBMET", 'F');
	factory->AddVariable("DR_ll", 'F');
	factory->AddVariable("Pt_Htt", 'F');
	factory->AddVariable("DR_HBBHTT", 'F');
	factory->AddVariable("Pt_Hbb", 'F');
	factory->AddVariable("DeltaPhi_METTT", 'F');
	factory->AddVariable("PtH", 'F');
	factory->AddVariable("mT2", 'F');
    
    FlatTree* sigTree = new FlatTree();
    FlatTree* bkgTree = new FlatTree();

    sigTree->Add(filePath+"ggH_hh_bbtautau_*.root");

    bkgTree->Add(filePath+"tt_*.root");
			
	outputFile->cd();
	
    ApplySelection(factory, sigTree, true);
    cout<<"********************* Background *********************"<<endl;
    ApplySelection(factory, bkgTree, false);
	
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
	
	// Save the output
	outputFile->Close();
	
	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;
	
	delete factory;
    delete sigTree;
    delete bkgTree;
}
