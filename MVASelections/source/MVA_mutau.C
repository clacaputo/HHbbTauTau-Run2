
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

#include "/Users/rosmina/root/tmva/test/TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

GetMVA(){
	std::cout << "==> Start TMVAClassification" << std::endl;
	
	
	// --------------------------------------------------------------------------------------------------
	
	TString outfileName( "/Users/rosmina/Desktop/out_mutau.root" );
	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );
	
	
	
    TFile *inputSignal  = TFile::Open("/Users/rosmina/Desktop/FLAT_TREE/mutau/ggH_tot.root");
    TFile *inputBkg1    = TFile::Open("/Users/rosmina/Desktop/FLAT_TREE/mutau/tt_tot.root");
	
		
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
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
    
	TTree* sigTree = (TTree*)inputSignal->Get("flatTree");
    TTree* bkgTree = (TTree*)inputBkg1->Get("flatTree");

	
	//TTree* sigTree = sigTree1->CopyTree("(againstMuonTight_2) && (againstElectronLooseMVA_2 ) && (byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5) && (pfRelIso_1<0.1) && (q_1*q_2<0)");
	//TTree* bkgTree = bkgTree1->CopyTree("(againstMuonTight_2) && (againstElectronLooseMVA_2 ) && (byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5) && (pfRelIso_1<0.1) && (q_1*q_2<0)");
	
	
	Double_t sigWeight = 1.;
	Double_t bkgWeight = 1.;
	
    std::vector<Double_t> vars( 13 );
	std::vector<Float_t> *pt_BJets=0;
	std::vector<Float_t> *eta_BJets=0;
	std::vector<Float_t> *phi_BJets=0;
	std::vector<Float_t> *energy_BJets=0;
	std::vector<Float_t> *csv_BJets=0;
    Float_t  treevars[5];
    Double_t Weight, rndm;
	Int_t nBjet;
	float mvamet_phi;
	float met;
	float drLT;
	float  ptsv;
	float  etasv;
	float  phisv;
	float eta1, eta2, phi1, phi2, e1, e2,byCombinedIsolationDeltaBetaCorrRaw3Hits2,pfRelIso1,mt1,mt2,weight ;
	int   Q1,Q2;
	
	Bool_t  againstMuonTight2,againstElectronLooseMVA2;
	// Signal
	sigTree->SetBranchAddress( "pt_1", &(treevars[0]));
	sigTree->SetBranchAddress( "pt_2", &(treevars[1]));
	
	sigTree->SetBranchAddress( "eta_1", &eta1);
	sigTree->SetBranchAddress( "eta_2", &eta2);
	sigTree->SetBranchAddress( "phi_1", &phi1);
	sigTree->SetBranchAddress( "phi_2", &phi2);
	sigTree->SetBranchAddress( "energy_1", &e1);
	sigTree->SetBranchAddress( "energy_2", &e2);
	
	sigTree->SetBranchAddress( "pt_Bjets", &pt_BJets);
	sigTree->SetBranchAddress( "eta_Bjets", &eta_BJets);
	sigTree->SetBranchAddress( "phi_Bjets", &phi_BJets);
	sigTree->SetBranchAddress( "energy_Bjets", &energy_BJets);
	sigTree->SetBranchAddress( "csv_Bjets", &csv_BJets);
	sigTree->SetBranchAddress("njetspt20", &nBjet);
	sigTree->SetBranchAddress("mvametphi", &mvamet_phi);
	sigTree->SetBranchAddress("mvamet", &met);
	sigTree->SetBranchAddress("DeltaR_leptons", &drLT);
	sigTree->SetBranchAddress("pt_sv", &ptsv);
	sigTree->SetBranchAddress("eta_sv", &etasv);
	sigTree->SetBranchAddress("phi_sv", &phisv);
	
	sigTree->SetBranchAddress("againstMuonTight_2",&againstMuonTight2);
	sigTree->SetBranchAddress("againstElectronLooseMVA_2",&againstElectronLooseMVA2);
	sigTree->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits2);
	
	sigTree->SetBranchAddress("pfRelIso_1",&pfRelIso1);
	sigTree->SetBranchAddress("q_1",&Q1);
	sigTree->SetBranchAddress("q_2",&Q2);
	sigTree->SetBranchAddress("mt_1",&mt1);
	sigTree->SetBranchAddress("mt_2",&mt2);
	sigTree->SetBranchAddress("weight",&weight);
	
	outputFile->cd();
	
	
	TH1F* hTauPt = new TH1F("hTauPt","hTauPt",200,0,200);
	TH1F* hMuonPt = new TH1F("hMuonPt","hMuonPt",200,0,200);
	TH1F* hPtb1 = new TH1F("hPtb1","hPtb1",300, 0, 300);
	TH1F* hptb2 = new TH1F("hPtb2","hPtb2",300, 0, 300);
	TH1F* hDRbb = new TH1F("hDRbb","hDRbb",100, 0, 10);
	TH1F* hDPhiBBMET = new TH1F("hDPhiBBMET","hDPhiBBMET",100, -3, 3);
	TH1F* hDRll = new TH1F("hDRll","hDRll",100, 0, 10);
	TH1F* hPtHtt = new TH1F("hPtHtt","hPtHtt",600, 0, 600);
	TH1F* hDRHBBTT= new TH1F("hDRHBBTT","hDRHBBTT",100, 0, 10);
	TH1F* hPtHBB= new TH1F("hPtHBB","hPtHBB",600, 0, 600);
	TH1F* hDeltaPhi_METTT= new TH1F("hDeltaPhi_METTT","hDeltaPhi_METTT",100, -3, 3);
	TH1F* hPtH= new TH1F("hPtH","hPtH",1000, 0, 1000);
	TH1F* hmT2= new TH1F("hmT2","hmT2",600, 0, 600);
	
	for (UInt_t i=0; i<(sigTree->GetEntries()); i++) {
		sigTree->GetEntry(i);
		if((againstMuonTight2) && (againstElectronLooseMVA2 ) && (byCombinedIsolationDeltaBetaCorrRaw3Hits2<1.5) && (pfRelIso1<0.1) && (Q1*Q2<0) && (mt1<30)){
			//if(!(nBjet>1)) continue;
			vars[2]=0;
			vars[3]=0;
			
			std::vector<int> IndexNjetEta2p4, IndexBjetsMedium;
			for (unsigned int k =0; k<nBjet; k++){
				if((TMath::Abs((*eta_BJets).at(k)) <2.4   )  ) {IndexNjetEta2p4.push_back(k);}
			}
			
			if ((IndexNjetEta2p4.size() > 1)) {
				
				
				for (unsigned int k =0; k<IndexNjetEta2p4.size() ; k++){
					//cout<<i<<"pt bjets  "<<(*pt_BJets).at(IndexNjetEta2p4.at(k)) <<" csv  "<<(*csv_BJets).at(IndexNjetEta2p4.at(k))<<endl;
					if(! ( (*pt_BJets).at(IndexNjetEta2p4.at(k)) > 20  )) continue;
					if(!(TMath::Abs((*eta_BJets).at(IndexNjetEta2p4.at(k))) <2.4   )  ) continue;
					if(! ((*csv_BJets).at(IndexNjetEta2p4.at(k)) > 0.689 )  ) continue;  
					IndexBjetsMedium.push_back(IndexNjetEta2p4.at(k));
				}
				
				bool BjetsCut =false;
				if(  IndexBjetsMedium.size() >1 )  BjetsCut = true;
				if ((BjetsCut == true)){
					
					TLorentzVector b1, b2, BB, MET, l1, l2, TT, H;
					b1.SetPtEtaPhiE((*pt_BJets).at(IndexBjetsMedium.at(0)), (*eta_BJets).at(IndexBjetsMedium.at(0)), (*phi_BJets).at(IndexBjetsMedium.at(0)) ,(*energy_BJets).at(IndexBjetsMedium.at(0)) );    
					b2.SetPtEtaPhiE((*pt_BJets).at(IndexBjetsMedium.at(1)), (*eta_BJets).at(IndexBjetsMedium.at(1)), (*phi_BJets).at(IndexBjetsMedium.at(1)) ,(*energy_BJets).at(IndexBjetsMedium.at(1)) ); 
					
					l1.SetPtEtaPhiE( treevars[0] , eta1, phi1, e1 );
					l2.SetPtEtaPhiE( treevars[1] , eta2, phi2, e2 );
					
					MET.SetPtEtaPhiE(met, 0., mvamet_phi, met);
					
					BB = b1+b2;
					TT=l2+l1;
					H = BB+TT;
					
					
					
					//////////////////////////////
					for (UInt_t ivar=0; ivar<2; ivar++) vars[ivar] = treevars[ivar];
					vars[2] = (*pt_BJets).at(IndexBjetsMedium.at(0));
					vars[3] = (*pt_BJets).at(IndexBjetsMedium.at(1));
					vars[4] = b1.DeltaR(b2);
					vars[5] = MET.DeltaPhi(BB);  //cout<<"  Delta phi MET*****************  "<<vars[5]<<endl;
					vars[6] = drLT;              // cout<<"  Delta leptoni *****************  "<<vars[6]<<endl;
					vars[7]=  TT.Pt();            //cout<<"  PT SV *****************  "<<vars[7]<<endl;
					vars[8] = TT.DeltaR(BB);
					vars[9] = BB.Pt();
					vars[10] = MET.DeltaPhi(TT);
					vars[11] = H.Pt();
					vars[12] = mt2;
					
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
					
					TRandom3 rand1;
					//cout<<" rnd "<<gRandom->Uniform(0, 1)<<endl;
					double r = gRandom->Uniform(0, 1);
					
					
					if (r<0.5)  {
						//cout<<"Riempi Training ev="<<i<<endl;
						factory->AddSignalTrainingEvent(vars, weight*sigWeight );} 
						
					else {
						//cout<<"Riempi Test ev="<<i<<endl;
						factory->AddSignalTestEvent(vars, weight*sigWeight );
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
	
	
	cout<<"********************* Background *********************"<<endl;
	
    bkgTree->SetBranchAddress( "pt_1", &(treevars[0]) );
    bkgTree->SetBranchAddress( "pt_2", &(treevars[1]) );
	
	bkgTree->SetBranchAddress( "eta_1", &eta1);
	bkgTree->SetBranchAddress( "eta_2", &eta2);
	bkgTree->SetBranchAddress( "phi_1", &phi1);
	bkgTree->SetBranchAddress( "phi_2", &phi2);
	bkgTree->SetBranchAddress( "energy_1", &e1);
	bkgTree->SetBranchAddress( "energy_2", &e2);
	
	
	
    bkgTree->SetBranchAddress( "pt_Bjets", &pt_BJets );
	bkgTree->SetBranchAddress( "eta_Bjets", &eta_BJets);
	bkgTree->SetBranchAddress( "phi_Bjets", &phi_BJets);
	bkgTree->SetBranchAddress( "energy_Bjets", &energy_BJets);
	bkgTree->SetBranchAddress( "csv_Bjets", &csv_BJets);
	bkgTree->SetBranchAddress("njetspt20", &nBjet);
	bkgTree->SetBranchAddress("mvametphi", &mvamet_phi);
	bkgTree->SetBranchAddress("mvamet", &met);
	bkgTree->SetBranchAddress("DeltaR_leptons", &drLT);
	
	bkgTree->SetBranchAddress("pt_sv", &ptsv);
	bkgTree->SetBranchAddress("eta_sv", &etasv);
	bkgTree->SetBranchAddress("phi_sv", &phisv);
	
	bkgTree->SetBranchAddress("againstMuonTight_2",&againstMuonTight2);
	bkgTree->SetBranchAddress("againstElectronLooseMVA_2",&againstElectronLooseMVA2);
	bkgTree->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits2);
	
	bkgTree->SetBranchAddress("pfRelIso_1",&pfRelIso1);
	bkgTree->SetBranchAddress("q_1",&Q1);
	bkgTree->SetBranchAddress("q_2",&Q2);
    bkgTree->SetBranchAddress("mt_1",&mt1);
	bkgTree->SetBranchAddress("mt_2",&mt2);
	bkgTree->SetBranchAddress("weight",&weight);
	
	
	for (UInt_t i=0; i<bkgTree->GetEntries(); i++) {
		//for (UInt_t i=0; i<2000; i++) {
        bkgTree->GetEntry(i);
		//		 if(!(nBjet>1)) continue;
        //cout<<tmp<<"\r";
		if((againstMuonTight2) && (againstElectronLooseMVA2 ) && (byCombinedIsolationDeltaBetaCorrRaw3Hits2<1.5) && (pfRelIso1<0.1) && (Q1*Q2<0) &&  (mt1<30)){
			
			vars[2]=0; vars[3]=0;
			
			std::vector<int> IndexNjetEta2p4, IndexBjetsMedium;
			for (unsigned int k =0; k<nBjet; k++){
				if((TMath::Abs((*eta_BJets).at(k)) <2.4   )  ) {IndexNjetEta2p4.push_back(k);}
			}
			
			if ((IndexNjetEta2p4.size() > 1)) {
				
						
				for (unsigned int k =0; k<IndexNjetEta2p4.size() ; k++){
					if(! ( (*pt_BJets).at(IndexNjetEta2p4.at(k)) > 20  )) continue;
					if(!(TMath::Abs((*eta_BJets).at(IndexNjetEta2p4.at(k))) <2.4   )  ) continue;
					if(! ((*csv_BJets).at(IndexNjetEta2p4.at(k)) > 0.689 )  ) continue;  
					IndexBjetsMedium.push_back(IndexNjetEta2p4.at(k));
				}
				
				bool BjetsCut =false;
				if(  IndexBjetsMedium.size() >1 )  BjetsCut = true;
				if ((BjetsCut == true)) {
					
					
				TLorentzVector b1, b2, BB, MET, l1, l2,TT, H;
				b1.SetPtEtaPhiE((*pt_BJets).at(IndexBjetsMedium.at(0)), (*eta_BJets).at(IndexBjetsMedium.at(0)), (*phi_BJets).at(IndexBjetsMedium.at(0)) ,(*energy_BJets).at(IndexBjetsMedium.at(0)) );    
				b2.SetPtEtaPhiE((*pt_BJets).at(IndexBjetsMedium.at(1)), (*eta_BJets).at(IndexBjetsMedium.at(1)), (*phi_BJets).at(IndexBjetsMedium.at(1)) ,(*energy_BJets).at(IndexBjetsMedium.at(1)) ); 
				
				
				l1.SetPtEtaPhiE( treevars[0] , eta1, phi1, e1 );
				l2.SetPtEtaPhiE( treevars[1] , eta2, phi2, e2 );
				
				
				MET.SetPtEtaPhiE(met, 0., mvamet_phi, met);
				BB = b1+b2;
				TT = l1+l2;
				H= TT+BB;
				////////////////////////////////////////////////
				for (UInt_t ivar=0; ivar<2; ivar++) vars[ivar] = treevars[ivar];
				vars[2] = (*pt_BJets).at(IndexBjetsMedium.at(0));
				vars[3] = (*pt_BJets).at(IndexBjetsMedium.at(1));
				vars[4] = b1.DeltaR(b2);
				vars[5] = MET.DeltaPhi(BB); //cout<<"  Delta phi MET*****************  "<<vars[5]<<endl;
				vars[6] = drLT;              //cout<<"  Delta leptoni *****************  "<<vars[6]<<endl;
				vars[7]=  TT.Pt();            //cout<<"  PT SV *****************  "<<vars[7]<<endl;
				vars[8] = TT.DeltaR(BB);
				vars[9] = BB.Pt();
				vars[10] = MET.DeltaPhi(TT); 
				vars[11] = H.Pt();
				vars[12] =mt2;
				
				//cout<<"pt_1 = "<<vars[0]<<endl;
				//cout<<"eta_1 = "<<vars[1]<<endl;
				//cout <<i<<"# bjet"<<nBjet<<endl; 
				
				//cout <<i<<" DR Bjet "<<vars[4]<<endl;	
					TRandom3 rand1;
					double r = gRandom->Uniform(0, 1);
					if (r<0.5)  {	factory->AddBackgroundTrainingEvent( vars, weight*bkgWeight );}
					else{factory->AddBackgroundTestEvent(vars, weight*bkgWeight );}
								
				}
				
			}
		}
    }
	
	
	// Set individual event weights (the variables must exist in the original TTree)
	//    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
	//    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
	//factory->SetBackgroundWeightExpression( "weight" );
	
	// Apply additional cuts on the signal and background samples (can be different)
	// TCut mycuts = "againstMuonTight_2>0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5 && pfRelIso_1<0.1 && q_1*q_2<0"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
	// TCut mycutb = "againstMuonTight_2>0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5 && pfRelIso_1<0.1 && q_1*q_2<0"; // for example: TCut mycutb = "abs(var1)<0.5";
	
	TCut mycuts = "";
	TCut mycutb = "";
	
    factory->PrepareTrainingAndTestTree(mycuts, mycutb,"SplitMode=Random" );
	// Tell the factory how to use the training and testing events
	//
	// If no numbers of events are given, half of the events in the tree are used
	// for training, and the other half for testing:
	//    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
	// To also specify the number of testing events, use:
	//    factory->PrepareTrainingAndTestTree( mycut,
	//                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
	//factory->PrepareTrainingAndTestTree( mycuts, mycutb,
	//                                   "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );
	
	// ---- Book MVA methods
	//
	// Please lookup the various method configuration options in the corresponding cxx files, eg:
	// src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
	// it is possible to preset ranges in the option string in which the cut optimisation should be done:
	// "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable
	
	
	// Cut optimisation
    
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
	
	
	
	
	// For an example of the category classifier usage, see: TMVAClassificationCategory
	
	// --------------------------------------------------------------------------------------------------
	
	// ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events
	
	// factory->OptimizeAllMethods("SigEffAt001","Scan");
	// factory->OptimizeAllMethods("ROCIntegral","FitGA");
	
	// --------------------------------------------------------------------------------------------------
	
	// ---- Now you can tell the factory to train, test, and evaluate the MVAs
	
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
	// Launch the GUI for the root macros
	if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
