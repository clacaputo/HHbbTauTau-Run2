#define VHAnalyzer_WH_New_cxx
#include "VHAnalyzer_WH_New.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <TLorentzVector.h>
#include <TMath.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>
using namespace std;


void VHAnalyzer_WH_New::Loop()
{
	
   if (fChain == 0) return;
	TFile* file=new TFile("/lustre/cms/store/user/rvenditti/PostLS2/FullSim/ggH2Tau_8TeV/iso.root","RECREATE");
	TTree* t1=new TTree("Tree","Iso_mu");
	float iso_mu, pt_mu;
	//	t1->Branch("isoMu",&iso_mu,"isoMu/F");
	t1->Branch("ptMu",&pt_mu,"ptMu/F");
	
	int puWeight=1;
    //reweighting pileup
    
   Long64_t nentries = fChain->GetEntriesFast();
   
    Int_t loose=0;
    /////////////////CICLO SU TUTTI GLI EVENTI
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //for (Long64_t jentry=0; jentry<3000; jentry++) { //for sui primi 100 eventi

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
    
	   ///////////////////////// SELEZIONATE LE PARTICELLE ED OSSERVATA LA FISICA DEL CANALE IN ANALISI FACCIO HISTO PRELIMINARI
       initial_entries->Fill(1);
       hnMuon->Fill(nMuon);
       hnTau->Fill(nTau);
       hnElectron->Fill(nElectron);
       hnMET->Fill(nMET);
       //hnMET_PU->Fill(nMET,puWeight);
       hMET->Fill(MET_met[0]);
       //hMET_PU->Fill(MET_met[0],puWeight);
	//cout<<nTau<<endl;
	/////////////////////////TAGLI VARI
	/// VERTEX CUT
	std::vector<int> indexVertex;
	Double_t sumPtmax_squared = 0;
	Int_t idmax = -1;
	Int_t NTVmax = -1;
	Int_t lmax = -1;
	Double_t zpv = -10;
	int nGoodVertices = 0;
	
	bool VertexCuts=false;
	
	// std::cout<<"Evento: "<<jentry<<std::endl;
	
	/*	for(Int_t l=0; l<nVertex; l++){
		VertexX->Fill(Vertex_x[l]);
		VertexY->Fill(Vertex_y[l]);
		VertexZ->Fill(Vertex_z[l]);
		VertexXY->Fill(Vertex_x[l],Vertex_y[l]);
		double z = Vertex_z[l];
		double rho = TMath::Sqrt(TMath::Power(Vertex_x[l],2)+TMath::Power(Vertex_y[l],2));
		if(Vertex_ndf[l] > 4 && fabs(z) < 24 && fabs(rho)<2)
		{
			++nGoodVertices;
			if (nGoodVertices==1)
				zpv = z;
		}
	}
	
	if (nGoodVertices!=0){VertexCuts =true;}
	
	if (VertexCuts==false) continue;
	std::cout<<jentry<<"good events for vertex selection "<<std::endl;
	entries_after_vertex->Fill(1);
	//entries_after_vertex_PU->Fill(1,puWeight);
	hnVertex->Fill(nVertex);
	//hnVertex_PU->Fill(nVertex, puWeight);
	*/
	  TLorentzVector MET;
	  MET.SetPtEtaPhiM(mvaMET_met[0], 0., mvaMET_metphi[0], 0.);	
	
	//////////////////////////// MUON CUT /////////////////////////////////////////
	std::vector<int> indexMuon;
	std::vector<int> indexMuonPt;
	std::vector<int> indexMuonAccept;
	std::vector<int> indexMuonID;
	std::vector<int> indexMuonISO;
	TLorentzVector mu; mu.SetPtEtaPhiE(1,0,0,0);
	Int_t muon_index=-1;
	TLorentzVector mu_tmp; mu_tmp.SetPtEtaPhiE(1,0,0,0);
	bool MuonCuts=false;
	bool MuonPt=false;
	Double_t muPt_max = 0;
	Int_t imax = -1;
	
	Double_t massa_trasversa=0;
	bool MtCut = true;
	Float_t dB_err = -10;
	Float_t dB = 10;
	std::vector<int>        indexMuonTrackGlob;
	std::vector<int>        indexMuonPixHits;
	std::vector<int>        indexMuonTrackHits;
	std::vector<int>        indexTrackD0;
	std::vector<int>        indexDistZ;
	std::vector<int>        indexNMatch;
	std::vector<int>        indexMuonChi2;
	std::vector<int>        indexMuonTrackdB;
	std::vector<int>        indexNChambers;
	std::vector<int>        indexNMatchStat;
	std::vector<int> relIso01, relIso02, relIso03, relIso04, relIso05, relIso06, relIso07, relIso08, relIso09, relIso1, relIso12, relIso14,relIso13,relIso16,relIso18;
	//std::cout<<" nMuon "<<nMuon<<std::endl;

	std::vector<int> indexNJet; bool OneJet=false;
	for (Int_t p=0; p<nJet; p++ ){
	  if (!( Jet_pt[p] > 30.0 )) continue;
	  if (!( TMath::Abs(Jet_eta[p]) < 4.7 )) continue;
		  indexNJet.push_back(p);
	  } //fine loop jet
	  if(indexNJet.size()==1){OneJet=true;}
	  if(OneJet==false)continue;
	  hOneJet->Fill(1);


	for (Int_t i=0 ; i<nMuon ; i++){
		mu_tmp.SetPtEtaPhiE(Muon_pt[i], Muon_eta[i], Muon_phi[i], Muon_energy[i]);
		indexMuonID.push_back(i);
		//		pt_mu=Muon_pt[i];
		//		t1->Fill();
		
		if(!(Muon_pt[i]>20))continue;
		indexMuonPt.push_back(i);
		if (!(TMath::Abs(Muon_eta[i]) < 2.1)) continue;
		indexMuonAccept.push_back(i);
		if (!((Muon_isTrackerMuon[i]==1)&&(Muon_isGlobalMuonPromptTight[i]==1) && (Muon_isPFMuon[i]==1)))continue;
		//if (!((Muon_isTrackerMuon[i]==1)&&(Muon_isGlobalMuonPromptTight[i]==1)))continue;
		indexMuonTrackGlob.push_back(i);
		if(!(Muon_nChambers[i]>=2))continue;
		indexNChambers.push_back(i);
		if(!(Muon_nMatches[i]>=2))continue;
		if(!(Muon_nMatchedStations[i]>=2))continue;
		indexNMatchStat.push_back(i);
		if(!(Muon_trackerLayersWithMeasurement[i]>5))continue;//5
		indexMuonTrackHits.push_back(i);
		if (!(Muon_pixHits[i]>=1)) continue;
		indexMuonPixHits.push_back(i);
		if(!(Muon_globalChi2[i]<10.0))continue;
		indexMuonChi2.push_back(i);
		if(!(TMath::Abs(Muon_dB[i])<0.045))continue;
		indexMuonTrackdB.push_back(i);
		//          if(!(TMath::Abs(Muon_vtxDistZ[i])<0.2))continue;
		//	if(!(TMath::Abs(Muon_vz[i]-zpv)<0.1))continue;
		indexDistZ.push_back(i);
		//if (!(Muon_pfRelIsoDB04v2[i] < 0.1)) continue;
		hrelIso->Fill(Muon_pfRelIso[i],puWeight);
		//iso_mu=Muon_pfRelIso[i];
		//		t1->Fill();

		if (!(Muon_pfRelIso[i] < 0.1)) continue;
		indexMuonISO.push_back(i);
		indexMuon.push_back(i);
		if(Muon_pt[i] > muPt_max){
		  muPt_max = Muon_pt[i];
		  imax = i;
		}
		hPtMuon->Fill(Muon_pt[i],puWeight);
		hEtaMuon->Fill(Muon_eta[i],puWeight);
		hPtEle->Fill(Muon_phi[i],puWeight);
	}
	
	// end for sui mu
	if (indexMuonAccept.size()!=0){
		sizeMuonAccept->Fill(1,puWeight);}
	if(indexMuonPt.size()!=0){
	  sizeMuonPt->Fill(1,puWeight);}
	if (indexMuonISO.size()!=0){
	  sizeMuonISO->Fill(1,puWeight);}
	
	if(indexMuonTrackGlob.size()>0){
	  sizeMuonTrackGlob->Fill(1,puWeight);
	}
	if(indexMuonPixHits.size()>0){
	  sizeMuonPixHits->Fill(1,puWeight);
	}
	if(indexMuonTrackHits.size()>0){
	  sizeMuonTrackHits->Fill(1,puWeight);
	}
	if(indexMuonChi2.size()>0){
	  sizeMuonChi2->Fill(1,puWeight);
	}
	
	if (indexMuonTrackdB.size()>0){
	  sizeMuondB->Fill(1,puWeight);
	}
	if(indexDistZ.size()>0){
	  sizeMuonDistZ->Fill(1,puWeight);
	}
	
	if(indexNChambers.size()>0){
	  sizeMuonChambers->Fill(1,puWeight);
	}
	
	if (indexMuonID.size()!=0){
	  sizeMuonID->Fill(1,puWeight);}
	
	if(indexNMatchStat.size()!=0){
	  sizeMuonMatchStat->Fill(1,puWeight);}
	if(indexMuon.size()!=0){
	  //sizeMuon->Fill(indexMuon.size(),puWeight);}
	  sizeMuon->Fill(1,puWeight);}
	
	
	if (indexMuon.size()>0){MuonCuts=true;}
	
	if(MuonCuts==false)continue;
	
	
	muon_index=indexMuon[0];
	mu.SetPtEtaPhiE(Muon_pt[muon_index], Muon_eta[muon_index], Muon_phi[muon_index],Muon_energy[muon_index] );
	//	  mu.SetPtEtaPhiM(Muon_pt[imax], Muon_eta[imax], Muon_phi[imax], 0.105);
	dB_err = Muon_dB[imax]/Muon_edB[imax];
	dB = Muon_dB[imax];
	hMuonPtReco->Fill(mu.Pt(),puWeight);
	hMuonEtaReco->Fill(mu.Eta(),puWeight);
	
	TLorentzVector mu1;
	mu1=mu;
	
	hnVertex2->Fill(nVertex);
	hnVertex_PU2->Fill(nVertex, puWeight);
	tagli[1]->Fill(1,puWeight);
	
	//	std::cout<<jentry<<" after muon  "<<endl;
	  
	  
	  //////////////////////////////////////// TAU CUT
	  
	  std::vector<int> indexTauAccept1;
	  std::vector<int> indexTauPt1;
	  std::vector<int> indexTauDeltaZ1;
	  std::vector<int> indexTauID1;
	  std::vector<int> indexTauISO1;
	  std::vector<int> indexTauAgainstLepton1;
	  std::vector<int> indexTauAccept2;
	  std::vector<int> indexTauPt2;
	  std::vector<int> indexTauDeltaZ2;
	  std::vector<int> indexTauID2;
	  std::vector<int> indexTauISO2;
	  std::vector<int> indexTauAgainstLepton2;
	  std::vector<int> indexTau1DR;
	  std::vector<int> indexTau1;
	  std::vector<int> indexTau2;
	  std::vector<int> indexTau2DRMu;
	  std::vector<int> indexTau2DRTau1;
	  std::vector<int> indexTau1AgainstMuon;
	  std::vector<int> indexTau2AgMu;
	  std::vector<int> indexTau1Fake;
	  std::vector<int> indexTau1Good;
	  std::vector<int> indexTau2Fake;
	  std::vector<int> indexTau2Good;
	  
	  std::vector<int> indexTauOS;
	  std::vector<int> indexTauSS;
	  float fattore_tau = 0;
	  std::vector<float> correzione;
	  float fattore_tau1 = 0;
	  std::vector<float> correzione1;
	  
	  Int_t tau1_indexprova=0; Int_t  tau2_indexprova=0; Int_t tau1_index=-1; Int_t tau2_index=-1;
	  TLorentzVector tau1; tau1.SetPtEtaPhiM(1,0,0,0);
	  TLorentzVector tau2; tau2.SetPtEtaPhiM(1,0,0,0);
	  Int_t chargeTauTau=-1;
	  bool TauCuts=false;
	  bool TauCuts1Fake=false;
	  bool TauCuts2Fake=false;
	  bool TaulooseCuts=false;
	  bool Tau1Pt=false;
	  bool Tau2Pt=false;
	  bool TauChargeCut=false;
	  TLorentzVector tau_tmp; tau_tmp.SetPtEtaPhiM(1,0,0,0);
	  Double_t tauPt_max = 0;
	  Int_t jmax = -1;
	  
	  bool isTauOS=false;
	  bool TauCuts1=false;
	  bool TauCuts2=false;
	  std::vector<int> indexEleIso;
	  
	  //std::cout<<" nTau "<<nTau<<std::endl;
	  
	  for(Int_t j=0; j<nTau; j++){
	    hTau_jetEta->Fill(Tau_jetEta[j],puWeight);
	    tau_tmp.SetPtEtaPhiM(Tau_pt[j], Tau_eta[j], Tau_phi[j], Tau_mass[j]);
	    indexTau1DR.push_back(j);
	    if (tau_tmp.DeltaR(mu)<=  0.5) continue;
	    indexTauDeltaZ1.push_back(j);
	    if(!(Tau_pt[j]>20)) continue;
	    indexTauPt1.push_back(j);
	    if (!(TMath::Abs(Tau_eta[j]) < 2.3)) continue;
	    indexTauAccept1.push_back(j);
	    if (!(Tau_decayModeFinding[j] == 1)) continue;
	    indexTauID1.push_back(j);
	    hTauIso_Bef->Fill(Tau_looseIsolation[j]); 
	    hTauIso_Aft->Fill((Tau_ptSumPFChargedHadronsIsoCone[j]+Tau_ptSumPFNeutralHadronsIsoCone[j]+Tau_ptSumPhotonsIsoCone[j])/(Tau_leadTrkPt[j])); 
	    if(!(Tau_byLooseIsolationDeltaBetaCorr[j] == 1)) continue;
	    //if (!(Tau_byLooseCombinedIsolationDeltaBetaCorr[j] == 1)) continue;

	    indexTauISO1.push_back(j);
	    if(! (Tau_againstMuonTight[j] == 1)) continue;
	    indexTau1AgainstMuon.push_back(j);
	    if (!(Tau_againstElectronLoose[j] == 1))continue;
	    indexTauAgainstLepton1.push_back(j);
	    //  if(!(TMath::Abs(Tau_zvertex[j]-zpv)<0.2)) continue;
	    indexTau1.push_back(j);
	    ZCut->Fill(tau_tmp.Pt());
		  
		  // std::cout<<" evento "<<jentry<<" tau 1 reco after selection pt "<<tau_tmp.Pt()<<" tau gen eta  "<<tau_tmp.Eta()<<" tau gen phi  "<<tau_tmp.Phi()<<" tau gen en  "<<tau_tmp.E()<<std::endl;
		  
	  } //for tau
	  
	  
	  if (indexTauAccept1.size()>0){
		  sizeTauAccept1->Fill(1,puWeight);
	  }
	  if(indexTau1DR.size()>0){
		  sizeTau1DR->Fill(1,puWeight);
		  
	  }
	  
	  if (indexTauPt1.size()>0){
		  sizeTauPt1->Fill(1,puWeight);
	  }
	  
	  if (indexTauISO1.size()>0){
		  sizeTauISO1->Fill(1,puWeight);
	  }
	  if (indexTauDeltaZ1.size()>0){
		  sizeTauDeltaZ1->Fill(1,puWeight);
	  }
	  if (indexTauID1.size()>0){
		  sizeTauID1->Fill(1,puWeight);
	  }
	  
	
	  if (indexTau1AgainstMuon.size()>0){
		  sizeTau1AgMu->Fill(1,puWeight);}
	  
	  if (indexTauAgainstLepton1.size()>0){
		  sizeTauAgainstLepton1->Fill(1,puWeight);
	  }
	  if (indexTau1.size()!=0){
		  nTau_baseSelection1->Fill(1,puWeight);}
	  
	  if (indexTau1.size()>0){TauCuts1 = true;}
	  
	  if(TauCuts1 ==false)continue;
	  
	  tagli[2]->Fill(1,puWeight);
	  
	  std::vector<int> indexZMu;
	  
	  indexZMu.push_back(1);
	  if(indexZMu.size()>0){sizeMuTau1Z->Fill(1,puWeight);}
	  
	  
	  hnVertex_PU3->Fill(nVertex, puWeight);
	  
	  
	  TLorentzVector  Z;
	  Z.SetPtEtaPhiM(1.,0.,0.,0.);
	  
	  Double_t m_MuMet=TMath::Sqrt(2*(mu1.Pt())*(MET.Pt())*(1-TMath::Cos(mu1.DeltaPhi(MET))));
	  MuTau1_AfterTau1->Fill(Z.M(),puWeight);
	  MuMET_AfterTau1->Fill(m_MuMet,puWeight);
	  Tau1Pt_AfterTau1->Fill(tau1.Pt(),puWeight);
	  MuPt_AfterTau1->Fill(mu1.Pt(),puWeight);  
	  
	  //	  cout<<jentry<<" numero muoni buoni: "<<indexMuon.size()<<"  numero tau:  "<<indexTau1.size()<<endl;
	  
	  
	  
	  std::vector<int>	  extramuon;
	  std::vector<int>  extrab;
	  std::vector<int> indexMtmumet;
	  std::vector<int> indexMuoncoppiaSegnale;
	  std::vector<int> indexTau1coppiaSegnale;
	  std::vector<float> pesoSegnale;
	  std::vector<int> indexCoppiaSegnale;
	  TLorentzVector  mu_charge, tau_charge;
	  mu_charge.SetPtEtaPhiM(1.,0.,0.,0.);
	  tau_charge.SetPtEtaPhiM(1.,0.,0.,0.);
	  
	  for(int t=0;t<indexMuon.size();t++){
		  mu_charge.SetPtEtaPhiM(Muon_pt[indexMuon[t]], Muon_eta[indexMuon[t]], Muon_phi[indexMuon[t]],Muon_energy[indexMuon[t]] );
		  for(int k=0;k<indexTau1.size();k++){
			  
			  tau_charge.SetPtEtaPhiM(Tau_pt[indexTau1[k]], Tau_eta[indexTau1[k]], Tau_phi[indexTau1[k]], Tau_energy[indexTau1[k]]);
			  
			  //   cout<<jentry<<"-----t1   "<<Tau_pt[indexTau1[t]]<<endl;
			  // cout<<jentry<<"-----t2  "<<Tau_pt[indexTau2[k]]<<endl;
			  
			  
			  if(!(mu_charge.DeltaR(tau_charge)>0.5)) continue;
			  if ((Tau_charge[indexTau1[k]]+Muon_charge[indexMuon[t]])==0 )isTauOS=true;
			  if(!(isTauOS==true)) continue;
			  double LT=mu_charge.Pt()+tau_charge.Pt();
			  indexMuoncoppiaSegnale.push_back(indexMuon[t]);
			  indexTau1coppiaSegnale.push_back(indexTau1[k]);
			  indexCoppiaSegnale.push_back(1);
		  }}
	  if(indexCoppiaSegnale.size()==0) continue;
	  
	  tagli[3]->Fill(1,puWeight);
	  //////////////////////////MT MUMET CUT///////////////////////////////////
	  double Mtmumet=TMath::Sqrt(2*(mu.Pt())*(MET.Pt())*(1-TMath::Cos(mu.DeltaPhi(MET))));
	  MtMuMetVsMet_afterCharge->Fill(Mtmumet,(MET_met[0]));
	  //  cout<<jentry<<" mt mumet  "<<Mtmumet<<endl;
	  if(!(Mtmumet<30))continue;
	  indexMtmumet.push_back(1);
	  sizeMtmumet->Fill(1,puWeight);
	  // cout<<jentry<<" AFTER MTMUMET "<<endl;
	  /////////////////////////zero goodMu/////////////////////////////////////////////////
	  TLorentzVector muDY1,muDY2,ZmumuCand; std::vector<int> indexZmumuCand;
	  for (Int_t i=0 ; i<nMuon ; i++){
	    if (!(TMath::Abs(Muon_eta[i]) < 2.4)) continue;
	    // if(!(TMath::Abs(Muon_vz[i]-zpv)<0.2))continue;                                                                                                                                        
	    if (!(Muon_pt[i] > 15)) continue;
	    if (!((Muon_isTrackerMuon[i]==1)&&(Muon_isGlobalMuonPromptTight[i]==1)&& (Muon_isPFMuon[i]==1)))continue;
	    if (!(Muon_pfRelIso[i] < 0.3)) continue;
	    if(!(TMath::Abs(Muon_dB[i])<0.045))continue;
	    indexZmumuCand.push_back(i);
          }
	  
	  bool isZmumu=false;
	  
	  for (Int_t i=0 ; i<indexZmumuCand.size() ; i++){
	    muDY1.SetPtEtaPhiE(Muon_pt[indexZmumuCand[i]], Muon_eta[indexZmumuCand[i]], Muon_phi[indexZmumuCand[i]], Muon_energy[indexZmumuCand[i]]);
	    for (Int_t j=0 ; j<indexZmumuCand.size() ; j++){
	      muDY2.SetPtEtaPhiE(Muon_pt[indexZmumuCand[j]], Muon_eta[indexZmumuCand[j]], Muon_phi[indexZmumuCand[j]], Muon_energy[indexZmumuCand[j]]);
	      if( (muDY1.DeltaR(muDY2) > 0.15)  && (Muon_charge[indexZmumuCand[i]]+Muon_charge[indexZmumuCand[j]] ==0   )  ){
		isZmumu=true;
		ZmumuCand=muDY1+muDY1;
		cout<<jentry<<"  "<<ZmumuCand.M()<<endl;
	      }
	    }
	  }
	  if (!(isZmumu==false)) continue;
	  hZmumuVeto->Fill(1);


	  std::vector<int> indexExtraMuon;
	  TLorentzVector muschifo;
	  muschifo.SetPtEtaPhiM(1.,0.,0.,0.);
	  for (Int_t i=0 ; i<nMuon ; i++){
		  muschifo.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], 0.105);
		  if (!(muschifo.DeltaR(mu)>0.5) ) continue;
		  if (!(TMath::Abs(Muon_eta[i]) < 2.4)) continue;
		  //  if(!(TMath::Abs(Muon_vz[i]-Tau_zvertex[indexTauLead[t]])<0.14))continue;
		  // if(!(TMath::Abs(Muon_vz[i]-zpv)<0.2))continue;
		  if (!(Muon_pt[i] > 10)) continue;
		  if (!((Muon_isTrackerMuon[i]==1)&&(Muon_isGlobalMuonPromptTight[i]==1)&& (Muon_isPFMuon[i]==1)))continue;
		  if (!(Muon_pfRelIso[i] < 0.3)) continue;
		  if(!(Muon_trackerLayersWithMeasurement[i]>5))continue;
		  if(!(Muon_globalChi2[i]<10.0))continue;
		  if(!(TMath::Abs(Muon_dB[i])<0.045))continue;
		 

		  indexExtraMuon.push_back(i);
	  }
	  if(!(indexExtraMuon.size()==0))continue;
	  extramuon.push_back(1);
	  NoExtraMuon->Fill(1,puWeight);
	  //cout<<jentry<<" AFTER EXTRA MU  "<<endl;
	  //////////////////////////////////jetCuts/////////////////////////////
	  std::vector<int> indexJet;
	  std::vector<int> indexJetMatch;
	  bool JetCuts=false;
	  bool JetCutsMatch=false;
	  TLorentzVector jet; jet.SetPtEtaPhiE(1,0,0,0);
	  
	  for (Int_t p=0; p<nJet; p++ ){
		  if (!( Jet_pt[p] > 30.0 )) continue;
		  if (!( TMath::Abs(Jet_eta[p]) < 2.4 )) continue;
		  if(!(Jet_combinedSecondaryVertexBTag[p]>0.244))continue;
		  bdiscrAleternative->Fill(Jet_simpleSecondaryVertexHighEffBTag[p],puWeight);
		  hPtJet_baseSelection->Fill(Jet_pt[p],puWeight);
		  indexJet.push_back(p);
		  
	  } //fine loop jet
	  nJet_baseSelection->Fill(indexJet.size(),puWeight);
	  if(indexJet.size()==0){JetCuts=true;}
	  
	  if(JetCuts==false)continue;
	  extrab.push_back(1);
	  sizeJetMatch->Fill(1,puWeight);
	  //	  cout<<jentry<<" AFTER EXTRA bjet  "<<endl;
	   ////////////////////////////////////////////ZERO ELETTRONI//////////////////////////////////////
		
	  std::vector<int> indexEle, extraele;
	  TLorentzVector el_prova;
	  el_prova.SetPtEtaPhiM(1,0,0,0);
	  bool ElectronCuts=false;
	  bool eleMin20=false;
	  bool eleGreat20=false;
	  
               
               for (Int_t h=0; h<nElectron; h++ ){
                   
                   TLorentzVector el_prova;
                   el_prova.SetPtEtaPhiE(1,0,0,0);
                   el_prova.SetPtEtaPhiM(Electron_pt[h], Electron_eta[h], Electron_phi[h], Electron_energy[h]);
                   
                   if (!(Electron_pt[h] > 10.0) ) continue;
                   // indexElePt.push_back(h);
                   
                   if (!(TMath::Abs(Electron_eta[h]) < 2.5) ) continue;
                   if (!(TMath::Abs(Electron_eta[h]) < 1.4442 || TMath::Abs(Electron_eta[h]) > 1.566)) continue; //crackVeto
		   if(!(TMath::Abs(Electron_vz[h]-zpv)<0.2))continue;
                   
		   if(!(Electron_missingHits[h] < 1)) continue;
		   if(!(Electron_hasMatchedConv[h] < 1)) continue;
		   if (!(TMath::Abs(Electron_dB[h]) < 0.045) ) continue;
                   
                   bool ele1_l20=false;
                   bool ele2_l20=false;
                   bool ele3_l20=false;
                   bool ele_l20=false;
                   //		cout<<h<< " of  "<<nElectron<<" idMVA "<<Electron_mvaPOGNonTrig[h]<<endl;
                   // replaced Electron_idMVA with Electron_mvaPOGNonTrig
                   if(Electron_pt[h] < 20 && TMath::Abs(Electron_scEta[h])<0.8 && Electron_mvaPOGNonTrig[h] > 0.925) {
                       ele1_l20=true;}
                   
                   if(Electron_pt[h] < 20 && (TMath::Abs(Electron_scEta[h])<1.479)  && TMath::Abs(Electron_scEta[h])>0.8 && (Electron_mvaPOGNonTrig[h] >  0.915) ) {
                       ele2_l20=true;}
                   
                   if(Electron_pt[h] < 20 && TMath::Abs(Electron_scEta[h])>1.479 &&   (Electron_mvaPOGNonTrig[h] >  0.965) ) {
                       ele3_l20=true;}
                   
                   if((ele1_l20 || ele2_l20 || ele3_l20)==true) ele_l20=true;
                   
                   bool ele1_g20=false;
                   bool ele2_g20=false;
                   bool ele3_g20=false;
                   bool ele_g20=false;
                   
                   if(Electron_pt[h] > 20 && TMath::Abs(Electron_scEta[h])<0.8 && (Electron_mvaPOGNonTrig[h] > 0.905 )) {
                       ele1_g20=true;}
                   
                   if(Electron_pt[h] > 20 && (TMath::Abs(Electron_scEta[h])<1.479 ) && (TMath::Abs(Electron_scEta[h])>0.8 ) && (Electron_mvaPOGNonTrig[h] >  0.955) ) {
                       ele2_g20=true;}
                   
                   if(Electron_pt[h] > 20 && TMath::Abs(Electron_scEta[h])>1.479 && (Electron_mvaPOGNonTrig[h] >   0.975) ) {
                       ele3_g20=true;
                   }
                   
                   if((ele1_g20 || ele2_g20 || ele3_g20)==true) ele_g20=true;
                   
                   bool eleFinal=false;
                   if(ele_l20 || ele_g20) eleFinal=true;
                   
                   if(eleFinal==false)continue;
                   indexEle.push_back(h);
               }
               
               
               if (!(indexEle.size()==0))continue;
               extraele.push_back(1);
               sizeZeroEle->Fill(1,puWeight);
	       cout<<jentry<<" AFTER EXTRA ELE  "<<endl;
	       //////////////////////////////////////////////////////////////////////////////////////////////
	      
	       for(int i=0; i<indexMuon.size(); i++){
		  pt_mu = Muon_pt[indexMuon[i]];
		 cout<<" muon pt "<<Muon_pt[indexMuon[i]]<<endl;
		 t1->Fill();         
	       }


	  
	  TLorentzVector MH;
	  for(Int_t c=0; c<indexCoppiaSegnale.size(); c++){
		  tau_charge.SetPtEtaPhiE(Tau_pt[indexTau1coppiaSegnale[c]], Tau_eta[indexTau1coppiaSegnale[c]], Tau_phi[indexTau1coppiaSegnale[c]], Tau_energy[indexTau1coppiaSegnale[c]]);
		  mu_charge.SetPtEtaPhiE(Muon_pt[indexMuoncoppiaSegnale[c]], Muon_eta[indexMuoncoppiaSegnale[c]], Muon_phi[indexMuoncoppiaSegnale[c]], Muon_energy[indexMuoncoppiaSegnale[c]]);
		  if(!(tau1.DeltaR(tau2)>0.5)) continue;
		  
		  cout<<jentry<<"-----final mu  "<<Tau_pt[indexMuoncoppiaSegnale[c]]<<endl;
		  cout<<jentry<<"-----final tau "<<Tau_pt[indexTau1coppiaSegnale[c]]<<endl;
		  
		  MH = mu_charge + tau_charge;

		  // std::cout<<jentry<<" peso segnale  "<<pesoSegnale.at(c)<<" peso  "<<peso<<" puWeight  "<<puWeight<<endl;;
		  hTau1Pt[7]->Fill(tau_charge.Pt(), puWeight);
		  hTau1Eta[7]->Fill(tau_charge.Eta(),puWeight);
		  hMuon1Pt[7]->Fill(mu_charge.Pt(),puWeight);
		  hMuon1Eta[7]->Fill(mu_charge.Eta(),puWeight);
		  massainv_TauTau[7]->Fill(MH.M(),puWeight);
		  tagli[7]->Fill(1,puWeight);


	  }




       // stampa qualcosa per vedere se gira...
       if (jentry % 10000 == 0 ){
           std::cout << "jentry = " << jentry << std::endl;}
} //parentesi del ciclo su tutti gli eventi
    
  cout<<" Initial  "<<initial_entries->Integral()<<endl;
  cout<<" Vertex cut  "<<entries_after_vertex->Integral()<<endl;
  cout<<" OneJet  "<<hOneJet->Integral()<<endl;
  cout<<" Muon cut  "<<tagli[1]->Integral()<<endl;	
  cout<<" Tau cut  "<<tagli[2]->Integral()<<endl;	
  cout<<" Pair Cut  "<<tagli[3]->Integral()<<endl;	
  cout<<" Mt(mu, MET) cut "<<sizeMtmumet->Integral()<<endl;	
  cout<<" Zmumu cut "<<hZmumuVeto->Integral()<<endl;	
  cout<<" Extra Muon Veto "<<NoExtraMuon->Integral()<<endl;
  cout<<" BJet Veto "<<sizeJetMatch->Integral()<<endl;
  cout<<" Extra ele veto "<<sizeZeroEle->Integral()<<endl;
	file->cd();
	//	t1->Write();
	file->Close();

    
} //parentesi del void

void VHAnalyzer_WH_New::bookHistograms() {
    if (!bookedHistos) {
        // Open Output File
        std::string root_file = fileName;
        outputFile = new TFile(root_file.c_str(), "RECREATE");
        outputFile->cd();
        TH1::SetDefaultSumw2();

        //////////////////////////////////////// H per conteggi /////////////////////////////////////////////
        mupt_gen= new TH1F(" mupt_gen"," mupt_gen" ,100, 0, 200);
        mass_MuMet= new TH1F("mass_MuMet","mass_MuMet" ,100, 0, 200);
        mass_MuTau1= new TH1F("mass_MuTau1","mass_MuTau1", 100, 0, 200);
        mass_TauTau_SS= new TH1F("mass_TauTau_SS","mass_TauTau_SS", 100, 0, 200);
        mass_MuMet_PU= new TH1F("mass_MuMet_PU","mass_MuMet_PU" ,100, 0, 200);
        mass_MuTau1_PU= new TH1F("mass_MuTau1_PU","mass_MuTau1_PU", 100, 0, 200);
        mass_TauTau_SS_PU= new TH1F("mass_TauTau_SS_PU","mass_TauTau_SS_PU", 100, 0, 200);
        MtMuMetVsMet_afterCharge= new TH2F("MtMuMetVsMet_afterCharge","MtMuMetVsMet_afterCharge", 200,0,200,200,0,200);
        MtMuMetVsMet_afterZ= new TH2F("MtMuMetVsMet_afterZ","MtMuMetVsMet_afterZ", 200,0,200,200,0,200);
        MtMuTauVSMET_afterCharge= new TH2F("MtMuTauVSMET_afterCharge","MtMuTauVSMET_afterCharge", 200,0,200,200,0,200);
        MtMuTauVSMET_afterZ= new TH2F("MtMuTauVSMET_afterZ","MtMuTauVSMET_afterZ", 200,0,200,200,0,200);
        cosTauLeadVSMtMuMet_afterCharge= new TH2F("cosTauLeadVSMtMuMet_afterCharge","cosTauLeadVSMtMuMet_afterCharge", 10,-1.1,1.1,200,0,200);
        cosTauSubVSMtMuMet_afterCharge= new TH2F("cosTauSubVSMtMuMet_afterCharge","cosTauSubVSMtMuMet_afterCharge", 10,-1.1,1.1,200,0,200);
        cosTauLeadVSMtMuMet_afterZ= new TH2F("cosTauLeadVSMtMuMet_afterZ","cosTauLeadVSMtMuMet_afterZ", 10,-1.1,1.1,200,0,200);
        cosTauSubVSMtMuMet_afterZ= new TH2F("cosTauSubVSMtMuMet_afterZ","cosTauSubVSMtMuMet_afterZ", 10,-1.1,1.1,200,0,200);
        
        cosTauLeadVSMet_afterCharge= new TH2F("cosTauLeadVSMet_afterCharge","cosTauLeadVSMet_afterCharge", 10,-1.1,1.1,200,0,200);
        cosTauSubVSMet_afterCharge= new TH2F("cosTauSubVSMet_afterCharge","cosTauSubVSMet_afterCharge", 10,-1.1,1.1,200,0,200);
        cosTauSubVSMet_afterZ=new TH2F("cosTauSubVSMet_afterZ","cosTauSubVSMet_afterZ", 10,-1.1,1.1,200,0,200);
        cosTauLeadVSMet_afterZ=new TH2F("cosTauLeadVSMet_afterZ","cosTauLeadVSMet_afterZ", 10,-1.1,1.1,200,0,200);
        
        hNPU = new TH1F("hNPU", "Number of PileUp", 36, -0.5, 35.5);
        hNPU_PU = new TH1F("hNPU_PU", "Number of PileUp ripesato", 36, -0.5, 35.5);
        hnVertex = new TH1F("hnVertex", "Number of Vertex", 80, 40, 120);
        hnVertex2 = new TH1F("hnVertex2", "Number of Vertex", 80, 40, 120);
        hnVertex_PU = new TH1F("hnVertex_PU", "Number of Vertex reweighted", 80, 40, 120);
        hnVertex_PU2 = new TH1F("hnVertex_PU2", "Number of Vertex reweighted", 80, 40, 120);
        hnVertex3 = new TH1F("hnVertex3", "Number of Vertex", 80, 40, 120);
        hnVertex_PU3 = new TH1F("hnVertex_PU3", "Number of Vertex reweighted", 80, 40, 120);
        hnMuon = new TH1F("nMuon","number of muon", 60, -0.5, 10);
        hnTau = new TH1F("nTau","number ef tau", 6, 97, 102);
        hnElectron = new TH1F("nElectron","number of ele", 31, -0.5, 10);
        hnMET_PU = new TH1F("nMET_PU","number of MET PU", 31, -0.5, 5);
        hnMET = new TH1F("nMET","number of MET", 31, -0.5, 5);
        hMET_PU = new TH1F("hMET_PU","histogram of MET PU", 50, 0, 150);
        hMET = new TH1F("hMET","histogram of MET", 50, 0, 150);
        hEtaMuon = new TH1F("EtaMuon","eta of muon", 50, -5, 5);
        hrelIso = new TH1F("RelIso","relIso of muon",100,0,20);
        hPtMuon = new TH1F("PtMuon","pt of muon candidate", 200, 0, 200);
        sizeMuonPt = new TH1F("sizeMuonPt", "Number of Muon after Pt cut", 31, -0.5, 30.5);
        sizeMuon = new TH1F("sizeMuon", "Number of Muon after cuts", 31, -0.5, 30.5);
        sizeMuonISO = new TH1F("sizeMuonISO", "Number of Muon after ISO cuts", 31, -0.5, 30.5);
        sizeMuonID = new TH1F("sizeMuonID", "Number of Muon after ID cuts", 31, -0.5, 30.5);
        sizeMuonAccept = new TH1F("sizeMuonAccept", "Number of Muon after Accept cuts", 31, -0.5, 30.5);
        hPtEle = new TH1F("hMuPhi","Phi of Mu Candidate", 70, 0, 3.5);
        nEle_baseSelection = new TH1F("nEle_baseSelection","number of electrons after baseline selection", 10, 0., 10);
        nTau_baseSelection1 = new TH1F("nTau_baseSelection1","number of tau1 after baseline selection", 10, 0., 10);
        sizeTauISO1 = new TH1F("sizeTauISO1","number of tau1 after ISO selection", 10, -0.5, 10.5);
        sizeTauID1 = new TH1F("sizeTauID1","number of tau1 after ID selection", 10,-0.5, 10.5);
        sizeTauDeltaZ1 = new TH1F("sizeTauDeltaZ1","number of tau1 after delta z selection", 10,-0.5, 10.5);
        sizeTauAccept1 = new TH1F("sizeTauAccept1","number of tau1 after Accept selection", 10, 0., 10);
        sizeTauPt1 = new TH1F("sizeTauPt1","number of tau1 after pt selection", 10, 0., 10);
        sizeTauAgainstLepton1 = new TH1F("sizeTauAgainstLepton1","number of tau1 after AgainstLepton selection", 10, 0., 10);
        nTau_baseSelection2 = new TH1F("nTau_baseSelection2","number of tau2 after baseline selection", 10, 0., 10);
        sizeTauISO2 = new TH1F("sizeTauISO2","number of tau2 after ISO selection", 10, -0.5, 10.5);
        sizeTauID2 = new TH1F("sizeTauID2","number of tau2 after ID selection", 10,-0.5, 10.5);
        sizeTauDeltaZ2 = new TH1F("sizeTauDeltaZ2","number of tau2 after delta z selection", 10,-0.5, 10.5);
        sizeTauAccept2 = new TH1F("sizeTauAccept2","number of tau2 after Accept selection", 10, 0., 10);
        sizeTauPt2 = new TH1F("sizeTauPt2","number of tau2 after pt selection", 10, 0., 10);
        sizeTauAgainstLepton2 = new TH1F("sizeTauAgainstLepton2","number of tau2 after AgainstLepton selection", 10, 0., 10);
        hTau_jetEta = new TH1F ("hTau_jetEta","eta of jet tau candidate", 50, -5, 5);
        hTau_deltaZ1 = new TH1F ("hTau_deltaZ2","delta Z of tau vertex and primary vertex ", 11, -0.1, 1);
        hTau_deltaZ2 = new TH1F ("hTau_deltaZ1","delta Z of tau vertex and primary vertex ", 11, -0.1, 1);
        hMuonEtaReco = new TH1F("hMuonEtaReco ","eta of muon ", 50, -5, 5);
        hMuonPtReco = new TH1F("hMuonPtReco ","pt of muon ", 100, 0, 200);
        hPtJet = new TH1F("PtJet","pt of jets", 100, 0, 200);
        hJettrackCountingHighEffBTag = new TH1F("hJet_trackCountingHighEffBTag","Jet_trackCountingHighEffBTag", 65,0,13);
        hPtJet_baseSelection = new TH1F("PtJet_baseSelection","pt of jets after baseline selection", 100, 0, 200);
        nJet_baseSelection = new TH1F("nJet_baseSelection","number of jets after baseline selection", 10, 0., 10);
        sizeJetMatch = new TH1F("sizeJetsMatch", "Number of event that has at least two good jets", 31, -0.5, 30.5);
        hWH_Pt = new TH1F("hWH_Pt","pt of mu + tau1+ tau2", 100, 0, 200);
        hTau1Pt_ottimizzazione = new TH1F("hTau1Pt_ottimizzazione","pt tau1", 100, 0, 100);
        hTau2Pt_ottimizzazione = new TH1F("hTau2Pt_ottimizzazione","pt tau2", 100, 0, 100);
        sizeTauZDist= new TH1F("sizeTauZDist","sizeTauZDist", 10, 0., 10);
        VertexX = new TH1F("VertexX","VertexX", 500, 0.1, 0.5);
        VertexY= new TH1F("VertexY","VertexY", 500, 0.1, 0.5);
        VertexZ= new TH1F("VertexZ","VertexZ", 100, 0, 30);
        VertexXY= new TH2F("VertexXY","VertexXY", 500, 0.1, 0.5, 500, 0.1, 0.5);
        MuTauVSHPt= new TH2F("MuTauVSHPt","MuTauVSHPt", 200, 0, 200, 200, 0, 200);
        sizeZeroEle = new TH1F("sizeZeroEle","sizeZeroEle", 10, 0., 10);
        ZCut = new TH1F("ZCut","ZCut", 50, 0., 100); ////////////////////<---------------------
        NoExtraMuon= new TH1F("NoExtraMuon","NoExtraMuon", 10, 0., 10);
        sizeMuTau1Z= new TH1F("sizeMuTau1Z","sizeMuTau1Z", 10, 0., 10);
        bdiscrAleternative= new TH1F("bdiscrAlternative","simpleSecondaryVertexHighEffBTag", 10, -5, 5);
        bdiscrAleternativeAfter= new TH1F("bdiscrAlternativeAfter","simpleSecondaryVertexHighEffBTag", 10, -5, 5);
        sizeMuonMatchStat= new TH1F("sizeMuonMatchStat","sizeMuonMatchStat", 10, 0., 10.);
        sizeTau1DR= new TH1F("sizeTau1DR","sizeTau1DR", 10,-0.5, 10.5);
        sizeTau1AgMu= new TH1F("sizeTau1AgMu","sizeTau1AgMu", 10,-0.5, 10.5);
        sizeTau2AgMu= new TH1F("sizeTau2AgMu","sizeTau2AgMu", 10,-0.5, 10.5);
        sizeMuonChambers= new TH1F("sizeMuonChambers","sizeMuonChambers", 10,-0.5, 10.5);
        sizeTau2DRTau1=new TH1F("sizeTau2DRTau1","sizeTau2DRTau1", 10,-0.5, 10.5);
        sizeTau2DRMu=new TH1F("sizeTau2DRMu","sizeTau2DRMu", 10,-0.5, 10.5);
        VertexSS=new TH1F("VertexSS","VertexSS", 10,-0.5, 10.5);
        ZeroMuSS=new TH1F("ZeroMuSS","ZeroMuSS", 10,-0.5, 10.5);
        zeroBJetSS=new TH1F("ZeroBJetSS","ZeroBJetSS", 10,-0.5, 10.5);
        ZeroEleSS=new TH1F("ZeroEleSS","ZeroEleSS", 10,-0.5, 10.5);
        AfterSS=new TH1F("AfterSS","AfterSS", 10,-0.5, 10.5);
        mass_TauTau_AfterSS_PU = new TH1F("mass_TauTau_AfterSS_PU","mass_TauTau_AfterSS_PU", 200, 0, 100);
        MuTau1_AfterTau1= new TH1F("MuTau1_AfterTau1"," MuTau1_AfterTau1", 200, 0, 100);
        MuMET_AfterTau1= new TH1F("MuMET_AfterTau1"," MuMET_AfterTau1", 200, 0, 100);
        Tau1Pt_AfterTau1= new TH1F("Tau1Pt_AfterTau1"," Tau1Pt_AfterTau1", 200, 0, 100);
        MuPt_AfterTau1= new TH1F("MuPt_AfterTau1"," MuPt_AfterTau1", 200, 0, 100);
        
        //////////////////////////////////////// H per conteggi /////////////////////////////////////////////
        H_tot = new TH1F("H_tot","H tot",3,0,3);
        W_tot = new TH1F("W_tot","W tot",3,0,3);
        Z_tot = new TH1F("Z_tot","Z tot",3,0,3);
        t_tot = new TH1F("t_tot","t tot",3,0,3);
        muonStatus1 = new TH1F("muonStatus1","muonStatus1",3,0,3);
        WplusMuon = new TH1F("WplusMuon","WplusMuon",3,0,3);
        muonW = new TH1F("muonW","",3,0,3);
        muonWstatus3 = new TH1F("muonWstatus3","",3,0,3);
        tautauEvent = new TH1F("tautauEvent","",3,0,3);
        controlloTauTau = new TH1F("controlloTauTau","",3,0,3);
        FHEvent= new TH1F("FHEvent","",3,0,3);
        GTFHEvent=new TH1F("GTFHEvent","",3,0,3);
        elejetEvent=new TH1F("elejetEvent","",3,0,3);
        GTelejetEvent=new TH1F("GTelejetEvent","",3,0,3);
        mujetEvent=new TH1F("mujetEvent","",3,0,3);
        GTmujetEvent=new TH1F("GTmujetEvent","",3,0,3);
        FLEvent=new TH1F("FLEvent","",3,0,3);
        GTFLEvent=new TH1F("GTFLEvent","",3,0,3);
        
        initial_entries=new TH1F("initial_entries","",3,0,3);;
        entries_after_vertex=new TH1F("entries_after_vertex","",3,0,3);;
        entries_after_vertex_PU=new TH1F("entries_after_vertex_PU","",3,0,3);;
        sizeAfterTrigger= new TH1F("sizeAfterTrigger", "sizeAfterTrigger", 3, 0, 3);
        sizeAfterTrigger_PU = new TH1F("sizeAfterTrigger_PU", "sizeAfterTrigger_PU", 3, 0, 3);
        
        MuonMatchEvent=new TH1F("MuonMatchEvent","",3,0,3);
        TauMatchEvent=new TH1F("TauMatchEvent","",3,0,3);
        JetMatchEvent=new TH1F("JetMatchEvent","",3,0,3);
        
        for(Int_t i=0; i<8; i++){
            char title[12];
            //itoa (i,title,10);
            sprintf(title,"tagli_%d",i);
            tagli[i] = new TH1F(title,"tagli",3,0,3);
            char title2[15];
            //itoa (i,title,10);
            sprintf(title2,"tagli_%d_PU",i);
            tagli_PU[i] = new TH1F(title2,"tagli_PU",3,0,3);
            
        }
        
        
        
        mass_Mu1Mu2= new TH1F("mass_Mu1Mu2","mass_Mu1Mu2", 100, 0, 200);
        mass_Mu1Tau1= new TH1F("mass_Mu1Tau1","mass_Mu1Tau1", 100, 0, 200);
        mass_Mu2Tau1= new TH1F("mass_Mu2Tau1","mass_Mu2Tau1", 100, 0, 200);
        mass_Mu1Tau2= new TH1F("mass_Mu1Tau2","mass_Mu1Tau2", 100, 0, 200);
        mass_Mu2Tau2= new TH1F("mass_Mu2Tau2","mass_Mu2Tau2", 100, 0, 200);
        
        mass_Mu1Mu2_jtotau1= new TH1F("mass_Mu1Mu2_jtotau1","mass_Mu1Mu2", 100, 0, 200);
        mass_Mu1Tau1_jtotau1= new TH1F("mass_Mu1Tau1_jtotau1","mass_Mu1Tau1", 100, 0, 200);
        mass_Mu2Tau1_jtotau1= new TH1F("mass_Mu2Tau1_jtotau1","mass_Mu2Tau1", 100, 0, 200);
        mass_Mu1Tau2_jtotau1= new TH1F("mass_Mu1Tau2_jtotau1","mass_Mu1Tau2", 100, 0, 200);
        mass_Mu2Tau2_jtotau1= new TH1F("mass_Mu2Tau2_jtotau1","mass_Mu2Tau2", 100, 0, 200);
        
        hTau1Eta_jtotau1 =new TH1F("hTau1Eta_jtotau1" ,"eta of tau ", 12, -3, 3);
        hTau1Pt_jtotau1 = new TH1F("hTau1Pt_jtotau1" ,"pt of tau ", 100, 0, 200);
        hTau2Eta_jtotau1 =new TH1F("hTau2Eta_jtotau1" ,"eta of tau ", 12, -3, 3);
        hTau2Pt_jtotau1 = new TH1F("hTau2Pt_jtotau1" ,"pt of tau ", 100, 0, 200);
        hMuon1Pt_jtotau1 = new TH1F("hMuon1Pt_jtotau1","pt of muon1 ", 100, 0, 200);
        hMuon1Eta_jtotau1 = new TH1F("hMuon1Eta_jtotau1" ,"eta of muon1 ", 12, -3, 3);
        hMuon2Pt_jtotau1 = new TH1F("hMuon2Pt_jtotau1" ,"pt of muon2 ", 100, 0, 200);
        hMuon2Eta_jtotau1 = new TH1F("hMuon2Eta_jtotau1" ,"eta of muon2 ", 12, -3, 3);
        massainv_TauTau_jtotau1 = new TH1F("massainv_TauTau_jtotau1","visible mass tau muon", 50, 0, 200);
        tagli_7_jtotau1 = new TH1F("tagli_7_jtotau1","tagli7",3,0,3);
        
        mass_Mu1Mu2_jtotau2= new TH1F("mass_Mu1Mu2_jtotau2","mass_Mu1Mu2", 100, 0, 200);
        mass_Mu1Tau1_jtotau2= new TH1F("mass_Mu1Tau1_jtotau2","mass_Mu1Tau1", 100, 0, 200);
        mass_Mu2Tau1_jtotau2= new TH1F("mass_Mu2Tau1_jtotau2","mass_Mu2Tau1", 100, 0, 200);
        mass_Mu1Tau2_jtotau2= new TH1F("mass_Mu1Tau2_jtotau2","mass_Mu1Tau2", 100, 0, 200);
        mass_Mu2Tau2_jtotau2= new TH1F("mass_Mu2Tau2_jtotau2","mass_Mu2Tau2", 100, 0, 200);
        hTau1Eta_jtotau2 =new TH1F("hTau1Eta_jtotau2" ,"eta of tau ", 12, -3, 3);
        hTau1Pt_jtotau2 = new TH1F("hTau1Pt_jtotau2" ,"pt of tau ", 100, 0, 200);
        hTau2Eta_jtotau2 =new TH1F("hTau2Eta_jtotau2" ,"eta of tau ", 12, -3, 3);
        hTau2Pt_jtotau2 = new TH1F("hTau2Pt_jtotau2" ,"pt of tau ", 100, 0, 200);
        hMuon1Pt_jtotau2 = new TH1F("hMuon1Pt_jtotau2","pt of muon1 ", 100, 0, 200);
        hMuon1Eta_jtotau2 = new TH1F("hMuon1Eta_jtotau2" ,"eta of muon1 ", 12, -3, 3);
        hMuon2Pt_jtotau2 = new TH1F("hMuon2Pt_jtotau2" ,"pt of muon2 ", 100, 0, 200);
        hMuon2Eta_jtotau2 = new TH1F("hMuon2Eta_jtotau2" ,"eta of muon2 ", 12, -3, 3);
        massainv_TauTau_jtotau2 = new TH1F("massainv_TauTau_jtotau2","visible mass tau muon", 50, 0, 200);
        tagli_7_jtotau2 = new TH1F("tagli_7_jtotau2","tagli7",3,0,3);
        
        mass_Mu1Mu2_DC= new TH1F("mass_Mu1Mu2_DC","mass_Mu1Mu2", 100, 0, 200);
        mass_Mu1Tau1_DC= new TH1F("mass_Mu1Tau1_DC","mass_Mu1Tau1", 100, 0, 200);
        mass_Mu2Tau1_DC= new TH1F("mass_Mu2Tau1_DC","mass_Mu2Tau1", 100, 0, 200);
        mass_Mu1Tau2_DC= new TH1F("mass_Mu1Tau2_DC","mass_Mu1Tau2", 100, 0, 200);
        mass_Mu2Tau2_DC= new TH1F("mass_Mu2Tau2_DC","mass_Mu2Tau2", 100, 0, 200);
        hTau1Eta_DC =new TH1F("hTau1Eta_DC" ,"eta of tau ", 12, -3, 3);
        hTau1Pt_DC = new TH1F("hTau1Pt_DC" ,"pt of tau ", 100, 0, 200);
        hTau2Eta_DC =new TH1F("hTau2Eta_DC" ,"eta of tau ", 12, -3, 3);
        hTau2Pt_DC = new TH1F("hTau2Pt_DC" ,"pt of tau ", 100, 0, 200);
        hMuon1Pt_DC = new TH1F("hMuon1Pt_DC","pt of muon1 ", 100, 0, 200);
        hMuon1Eta_DC = new TH1F("hMuon1Eta_DC" ,"eta of muon1 ", 12, -3, 3);
        hMuon2Pt_DC = new TH1F("hMuon2Pt_DC" ,"pt of muon2 ", 100, 0, 200);
        hMuon2Eta_DC = new TH1F("hMuon2Eta_DC" ,"eta of muon2 ", 12, -3, 3);
        massainv_TauTau_DC = new TH1F("massainv_TauTau_DC","visible mass tau muon", 50, 0, 200);
        tagli_7_DC = new TH1F("tagli_7_DC","tagli7",3,0,3);
        
        
        
        for(Int_t i=0; i<8; i++){
            char titleptTau1[25];
            char titleetaTau1[25];
            char titleptTau2[25];
            char titleetaTau2[25];
            char titleptMu1[25];
            char titleetaMu1[25];
            char titleptMu2[25];
            char titleetaMu2[25];
            char titleminv[25];
            char titledB_err[25];
            char titledB[25];
            sprintf(titleetaTau1,"hTau1Eta_taglio%d",i);
            sprintf(titleptTau1,"hTau1Pt_taglio%d",i);
            sprintf(titleetaTau2,"hTau2Eta_taglio%d",i);
            sprintf(titleptTau2,"hTau2Pt_taglio%d",i);
            sprintf(titleetaMu1,"hMuon1Eta_taglio%d",i);
            sprintf(titleptMu1,"hMuon1Pt_taglio%d",i);
            sprintf(titleetaMu2,"hMuon2Eta_taglio%d",i);
            sprintf(titleptMu2,"hMuon2Pt_taglio%d",i);
            sprintf(titleminv,"massainv_TauTau_taglio%d",i);
            sprintf(titledB_err,"dB_edB_taglio%d",i);
            sprintf(titledB,"MuondB_taglio%d",i);
            
            char titleptTau1_PU[25];
            char titleetaTau1_PU[25];
            char titleptTau2_PU[25];
            char titleetaTau2_PU[25];
            char titleptMu_PU[25];
            char titleetaMu_PU[25];
            char titleminv_PU[25];
            char titledB_err_PU[30];
            char titledB_PU[25];
            sprintf(titleetaTau1_PU,"hTau1Eta_taglio%d_PU",i);
            sprintf(titleptTau1_PU,"hTau1Pt_taglio%d_PU",i);
            sprintf(titleetaTau2_PU,"hTau2Eta_taglio%d_PU",i);
            sprintf(titleptTau2_PU,"hTau2Pt_taglio%d_PU",i);
            sprintf(titleetaMu_PU,"hMuonEta_taglio%d_PU",i);
            sprintf(titleptMu_PU,"hMuonPt_taglio%d_PU",i);
            sprintf(titleminv_PU,"massainv_TauTau_taglio%d_PU",i);
            sprintf(titledB_err_PU,"dB_edB_taglio%d_PU",i);
            sprintf(titledB_PU,"MuondB_taglio%d_PU",i);
            
            
            hTau1Eta[i]=new TH1F(titleetaTau1 ,"eta of tau ", 12, -3, 3);
            hTau1Pt[i]= new TH1F(titleptTau1 ,"pt of tau ", 100, 0, 200);
            hTau2Eta[i]=new TH1F(titleetaTau2 ,"eta of tau ", 12, -3, 3);
            hTau2Pt[i]= new TH1F(titleptTau2 ,"pt of tau ", 100, 0, 200);
            hMuon1Pt[i]= new TH1F(titleptMu1 ,"pt of muon1 ", 100, 0, 200);
            hMuon1Eta[i]= new TH1F(titleetaMu1 ,"eta of muon1 ", 12, -3, 3);
            hMuon2Pt[i]= new TH1F(titleptMu2 ,"pt of muon2 ", 100, 0, 200);
            hMuon2Eta[i]= new TH1F(titleetaMu2 ,"eta of muon2 ", 12, -3, 3);
            massainv_TauTau[i] = new TH1F(titleminv,"visible mass tau muon", 200, 0, 400);
            hdB_err[i]=new TH1F(titledB_err," dB/dB_err per mu",100,-3,10);
            muondB[i] = new TH1F(titledB,"db of muon ",100,-0.1,0.1);
            
            hTau1Eta_PU[i]=new TH1F(titleetaTau1_PU ,"eta of tau ", 12, -3, 3);
            hTau1Pt_PU[i]= new TH1F(titleptTau1_PU ,"pt of tau ", 100, 0, 200);
            hTau2Eta_PU[i]=new TH1F(titleetaTau2_PU ,"eta of tau ", 12, -3, 3);
            hTau2Pt_PU[i]= new TH1F(titleptTau2_PU ,"pt of tau ", 100, 0, 200);
            hMuonPt_PU[i]= new TH1F(titleptMu_PU ,"pt of muon ", 100, 0, 200);
            hMuonEta_PU[i]= new TH1F(titleetaMu_PU ,"eta of muon ", 12, -3, 3);
            massainv_TauTau_PU[i] = new TH1F(titleminv_PU,"visible mass tau muon", 200, 0, 400);
            hdB_err_PU[i]=new TH1F(titledB_err_PU," dB/dB_err per mu",100,-3,10);
            muondB_PU[i] = new TH1F(titledB_PU,"db of muon ",100,-0.1,0.1);
            
            
        }
        hmassa_trasversa = new TH1F("hmassa_trasversa", "massa trasversa",12,0,120);
        deltaPhi_muMET = new TH1F("deltaPhi_muMET","DeltaPhi tra mu e MET", 20, -1.5, 1.5);
        hPt_ratio = new TH1F("pt_ratio","ratio of pt",100,-0.1,1.1);
        hPt_ratio_PU = new TH1F("pt_ratio_PU","ratio of pt",100,-0.1,1.1);
        hPt_ratio_Final = new TH1F("pt_ratio_Final","ratio of pt",100,-0.1,1.1);
        hPt_ratio_Final_PU = new TH1F("pt_ratio_Final_PU","ratio of pt",100,-0.1,1.1);
        
        sizeMuonTrackGlob= new TH1F("sizeMuonTrackGlob", "sizeMuonTrackGlob", 31, -0.5, 30.5);
        sizeMuonPixHits= new TH1F("sizeMuonPixHits", "sizeMuonPixHits", 31, -0.5, 30.5);
        sizeMuonTrackHits= new TH1F("sizeMuonTrackHits", "sizeMuonTrackHits", 31, -0.5, 30.5);
        sizeMuonChi2= new TH1F("sizeMuonChi2", "Number of Muon after chi2", 31, -0.5, 30.5);
        sizeMuondB= new TH1F("sizeMuondB", "Number of Muon after dB cuts", 31, -0.5, 30.5);
        sizeMuonDistZ= new TH1F("sizeMuonDistZ", "Number of Muon after distZ cuts", 31, -0.5, 30.5);
        sizeMuonNMatch= new TH1F("sizeMuonNMatch", "Number of Muon after NMatch", 31, -0.5, 30.5);
        sizeMtmumet= new TH1F("sizeMtmumet", "sizeMtmumet", 31, -0.5, 30.5);
		hTauPairCuts= new TH1F("hTauPairCuts", "hTauPairCuts", 31, -0.5, 30.5);
		MaxPtGreat25= new TH1F("MaxPtGreat25", "MaxPtGreat25", 31, -0.5, 30.5);
		sizeTauSS= new TH1F("sizeTauSS", "sizeTauSS", 31, -0.5, 30.5);
		sizeTauOS= new TH1F("sizeTauOS", "sizeTauOS", 31, -0.5, 30.5);
		hBDT_Fake= new TH1F("hBDT_Fake", "hBDT_Fake", 10, -1, 1);
		hMLP_WZ= new TH1F("hMLP_WZ", "hMLP_WZ", 30, -3, 3);

		hMuonIsoOpt=new TH1F ("hMuonIsoOpt","hMuonIsoOpt",20,0,20);
		hTauIso_Aft=new TH1F ("hTauIso_Aft","hTauIso_Aft",100,0,20);
		hTauIso_Bef=new TH1F ("hTauIso_Bef","hTauIso_Bef",100,0,20);
		hOneJet=new TH1F ("hOneJet","hOneJet",10,0,10);
		hZmumuVeto=new TH1F ("hZmumuVeto","hZmumuVeto",10,0,10);
    }
}


void VHAnalyzer_WH_New::saveHistograms() {
    if (outputFile) {
        //f_puweights->Close();
        std::string root_file = fileName;
        outputFile = new TFile(root_file.c_str(), "RECREATE");
        outputFile->cd();
        
        hNPU->Write();
        hNPU_PU->Write();
        hnVertex->Write();
        hnVertex2->Write();
        hnVertex_PU->Write();
        hnVertex_PU2->Write();
	VertexX->Write();
	VertexY->Write();
	VertexZ->Write();
	VertexXY->Write();
	hOneJet->Write();
        hnElectron->Write();
        hnTau->Write();
        hnMuon->Write();
        hnMET->Write();
        hnMET_PU->Write();
        hMET->Write();
        hMET_PU->Write();
        hPtMuon->Write();
        hEtaMuon->Write();
        hrelIso->Write();
        sizeMuon->Write();
        sizeMuonAccept->Write();
        sizeMuonISO->Write();
        sizeMuonTrackGlob->Write();
        sizeMuonPixHits->Write();
        sizeMuonTrackHits->Write();
        sizeMuonChi2->Write();
        sizeMuondB->Write();
        sizeMuonDistZ->Write();
        sizeMuonNMatch->Write();
        sizeMuonChambers->Write();
        sizeMuonID->Write();
        sizeMuonMatchStat->Write();
        sizeTau1DR->Write();
        sizeTauAgainstLepton1->Write();
        sizeTauISO1->Write();
        sizeTau1AgMu->Write();
        sizeTauDeltaZ1->Write();
        sizeTauID1->Write();
        sizeTauAccept1->Write();
        sizeTauPt1->Write();
        sizeTauAgainstLepton2->Write();
        sizeTauISO2->Write();
        sizeTauDeltaZ2->Write();
        sizeTauID2->Write();
        sizeTauAccept2->Write();
        sizeTauPt2->Write();
        sizeTau2DRTau1->Write();
        sizeTau2DRMu->Write();
        sizeTau2AgMu->Write();
        hMuonPtReco->Write();
        hMuonEtaReco->Write();
        hPtEle->Write();
        nEle_baseSelection->Write();
        nTau_baseSelection1->Write();
        nTau_baseSelection2->Write();
        hTau_jetEta->Write();
        hTau_deltaZ1->Write();
        hTau_deltaZ2->Write();
        
        hWH_Pt->Write();
	    
        MuTauVSHPt->Write();
        H_tot->Write();
        W_tot->Write();
        Z_tot->Write();
        t_tot->Write();
        muonStatus1->Write();
        WplusMuon->Write();
        muonW->Write();
        muonWstatus3->Write();
        tautauEvent->Write();
        controlloTauTau->Write();
        FHEvent->Write();
        GTFHEvent->Write();
        elejetEvent->Write();
        GTelejetEvent->Write();
        mujetEvent->Write();
        GTmujetEvent->Write();
        FLEvent->Write();
        GTFLEvent->Write();
        initial_entries->Write();
        entries_after_vertex->Write();
        entries_after_vertex_PU->Write();
        sizeAfterTrigger->Write();
        sizeAfterTrigger_PU->Write();
        
        sizeMuonPt->Write();
        
        
        
        sizeTauZDist->Write();
        sizeMuTau1Z->Write();
        sizeZeroEle->Write();
        sizeJetMatch->Write();
        ZCut->Write();
        NoExtraMuon->Write();
        
        mass_MuMet->Write();
        mass_MuTau1->Write();
        mass_TauTau_SS->Write();
        mass_MuMet_PU->Write();
        mass_MuTau1_PU->Write();
        mass_TauTau_SS_PU->Write();
        
        
        MuTau1_AfterTau1->Write();
        MuMET_AfterTau1->Write();
        Tau1Pt_AfterTau1->Write();
        MuPt_AfterTau1->Write();
        
        
        hTau1Pt[7]->Write();
        hTau1Eta[7]->Write();
        hTau2Pt[7]->Write();
        hTau2Eta[7]->Write();
        hMuon1Pt[7]->Write();
        hMuon2Pt[7]->Write();
        hMuon1Eta[7]->Write();
        hMuon2Eta[7]->Write();
        massainv_TauTau[7]->Write();
        mass_Mu1Mu2->Write();
        mass_Mu1Tau1->Write();
        mass_Mu2Tau1->Write();
        mass_Mu1Tau2->Write();
        mass_Mu2Tau2->Write();
        tagli[1]->Write();
	tagli[2]->Write();
	tagli[3]->Write();
        tagli[7]->Write();
        tagli[6]->Write();
        tagli[4]->Write();
        tagli[0]->Write();
        
        MtMuMetVsMet_afterCharge->Write();
        
        
        MtMuMetVsMet_afterZ->Write();
        MtMuTauVSMET_afterCharge->Write();
        MtMuTauVSMET_afterZ->Write();
        cosTauLeadVSMtMuMet_afterCharge->Write();
        cosTauSubVSMtMuMet_afterCharge->Write();
        cosTauLeadVSMtMuMet_afterZ->Write();
        cosTauSubVSMtMuMet_afterZ->Write();
        
        cosTauLeadVSMet_afterCharge->Write();
        cosTauSubVSMet_afterCharge->Write();
        cosTauSubVSMet_afterZ->Write();
        cosTauLeadVSMet_afterZ->Write();
        mupt_gen->Write();
        sizeMtmumet->Write();
        hTauPairCuts->Write();
	hMuonIsoOpt->Write();
		MaxPtGreat25->Write();
		sizeTauOS->Write();
		sizeTauSS->Write();
		hBDT_Fake->Write();
		hMLP_WZ->Write();
		hZmumuVeto->Write();
		hTauIso_Aft->Write();
		hTauIso_Bef->Write();
        //delete outputFile;    
	}
}
