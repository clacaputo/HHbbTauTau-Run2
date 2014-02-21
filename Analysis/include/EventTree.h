//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 21 14:20:54 2014 by ROOT version 5.32/00
// from TTree vhtree/VH Analysis Tree
// found on file: /data/mc/analysis/Tree_RadionM500.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <TClonesArray.h>
#include <TObject.h>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxEvent = 1;
const Int_t kMaxVertex = 100;
const Int_t kMaxElectron = 20;
const Int_t kMaxGenParticle = 2000;
const Int_t kMaxJet = 100;
const Int_t kMaxMET = 1;
const Int_t kMaxcorrMET = 1;
const Int_t kMaxmvaMET = 1;
const Int_t kMaxMuon = 20;
const Int_t kMaxTau = 300;

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Event_;
   UInt_t          Event_fUniqueID[kMaxEvent];   //[Event_]
   UInt_t          Event_fBits[kMaxEvent];   //[Event_]
   UInt_t          Event_run[kMaxEvent];   //[Event_]
   UInt_t          Event_event[kMaxEvent];   //[Event_]
   UInt_t          Event_lumis[kMaxEvent];   //[Event_]
   UInt_t          Event_bunch[kMaxEvent];   //[Event_]
   UInt_t          Event_orbit[kMaxEvent];   //[Event_]
   Double_t        Event_time[kMaxEvent];   //[Event_]
   Bool_t          Event_isdata[kMaxEvent];   //[Event_]
   Bool_t          Event_isPhysDeclared[kMaxEvent];   //[Event_]
   Bool_t          Event_isBPTX0[kMaxEvent];   //[Event_]
   Bool_t          Event_isBSCMinBias[kMaxEvent];   //[Event_]
   Bool_t          Event_isBSCBeamHalo[kMaxEvent];   //[Event_]
   Bool_t          Event_isPrimaryVertex[kMaxEvent];   //[Event_]
   Bool_t          Event_isBeamScraping[kMaxEvent];   //[Event_]
   Bool_t          Event_passHBHENoiseFilter[kMaxEvent];   //[Event_]
   vector<int>     Event_nPU[kMaxEvent];
   vector<int>     Event_bunchCrossing[kMaxEvent];
   vector<int>     Event_trueNInt[kMaxEvent];
   vector<int>     *nPU;
   vector<int>     *bunchCrossing;
   vector<int>     *trueNInt;
   Int_t           Vertex_;
   UInt_t          Vertex_fUniqueID[kMaxVertex];   //[Vertex_]
   UInt_t          Vertex_fBits[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_x[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_y[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_z[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_xErr[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_yErr[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_zErr[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_rho[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_chi2[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_ndf[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_ntracks[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_ntracksw05[kMaxVertex];   //[Vertex_]
   Bool_t          Vertex_isfake[kMaxVertex];   //[Vertex_]
   Bool_t          Vertex_isvalid[kMaxVertex];   //[Vertex_]
   Double_t        Vertex_sumPt[kMaxVertex];   //[Vertex_]
   Int_t           Vertex_selbit[kMaxVertex];   //[Vertex_]
   Int_t           nVertex;
   Int_t           Electron_;
   UInt_t          Electron_fUniqueID[kMaxElectron];   //[Electron_]
   UInt_t          Electron_fBits[kMaxElectron];   //[Electron_]
   Double_t        Electron_eta[kMaxElectron];   //[Electron_]
   Double_t        Electron_phi[kMaxElectron];   //[Electron_]
   Double_t        Electron_pt[kMaxElectron];   //[Electron_]
   Bool_t          Electron_ecalDriven[kMaxElectron];   //[Electron_]
   Bool_t          Electron_hasGsfTrack[kMaxElectron];   //[Electron_]
   Double_t        Electron_trackPt[kMaxElectron];   //[Electron_]
   Double_t        Electron_trackPtError[kMaxElectron];   //[Electron_]
   Double_t        Electron_energy[kMaxElectron];   //[Electron_]
   Double_t        Electron_caloEnergy[kMaxElectron];   //[Electron_]
   Double_t        Electron_caloEnergyError[kMaxElectron];   //[Electron_]
   Int_t           Electron_charge[kMaxElectron];   //[Electron_]
   Int_t           Electron_pixHits[kMaxElectron];   //[Electron_]
   Int_t           Electron_trkHits[kMaxElectron];   //[Electron_]
   Int_t           Electron_nValidHits[kMaxElectron];   //[Electron_]
   Float_t         Electron_trkD0[kMaxElectron];   //[Electron_]
   Float_t         Electron_trkD0Error[kMaxElectron];   //[Electron_]
   Double_t        Electron_hoe[kMaxElectron];   //[Electron_]
   Double_t        Electron_hoeDepth1[kMaxElectron];   //[Electron_]
   Double_t        Electron_eop[kMaxElectron];   //[Electron_]
   Double_t        Electron_sigmaEtaEta[kMaxElectron];   //[Electron_]
   Double_t        Electron_sigmaIEtaIEta[kMaxElectron];   //[Electron_]
   Double_t        Electron_deltaPhiTrkSC[kMaxElectron];   //[Electron_]
   Double_t        Electron_deltaEtaTrkSC[kMaxElectron];   //[Electron_]
   Int_t           Electron_classif[kMaxElectron];   //[Electron_]
   Double_t        Electron_e1x5overe5x5[kMaxElectron];   //[Electron_]
   Double_t        Electron_e2x5overe5x5[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoEcal03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoHcal03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoTrk03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoEcal04[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoHcal04[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoTrk04[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoRel03[kMaxElectron];   //[Electron_]
   Double_t        Electron_isoRel04[kMaxElectron];   //[Electron_]
   Double_t        Electron_vx[kMaxElectron];   //[Electron_]
   Double_t        Electron_vy[kMaxElectron];   //[Electron_]
   Double_t        Electron_vz[kMaxElectron];   //[Electron_]
   Double_t        Electron_scEn[kMaxElectron];   //[Electron_]
   Double_t        Electron_scEta[kMaxElectron];   //[Electron_]
   Double_t        Electron_scPhi[kMaxElectron];   //[Electron_]
   Double_t        Electron_scET[kMaxElectron];   //[Electron_]
   Double_t        Electron_scRawEnergy[kMaxElectron];   //[Electron_]
   Double_t        Electron_vtxDist3D[kMaxElectron];   //[Electron_]
   Int_t           Electron_vtxIndex[kMaxElectron];   //[Electron_]
   Double_t        Electron_vtxDistZ[kMaxElectron];   //[Electron_]
   Double_t        Electron_relIso[kMaxElectron];   //[Electron_]
   Double_t        Electron_pfRelIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_chargedHadronIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_neutralHadronIso[kMaxElectron];   //[Electron_]
   Float_t         Electron_photonIso[kMaxElectron];   //[Electron_]
   Int_t           Electron_missingHits[kMaxElectron];   //[Electron_]
   Double_t        Electron_dB[kMaxElectron];   //[Electron_]
   Double_t        Electron_edB[kMaxElectron];   //[Electron_]
   Double_t        Electron_dB3d[kMaxElectron];   //[Electron_]
   Double_t        Electron_edB3d[kMaxElectron];   //[Electron_]
   Double_t        Electron_scE1E9[kMaxElectron];   //[Electron_]
   Double_t        Electron_scS4S1[kMaxElectron];   //[Electron_]
   Double_t        Electron_sckOutOfTime[kMaxElectron];   //[Electron_]
   Double_t        Electron_scEcalIso[kMaxElectron];   //[Electron_]
   Double_t        Electron_scHEEPEcalIso[kMaxElectron];   //[Electron_]
   Double_t        Electron_scHEEPTrkIso[kMaxElectron];   //[Electron_]
   Int_t           Electron_nBrems[kMaxElectron];   //[Electron_]
   Float_t         Electron_fbrem[kMaxElectron];   //[Electron_]
   Float_t         Electron_dist_vec[kMaxElectron];   //[Electron_]
   Float_t         Electron_dCotTheta[kMaxElectron];   //[Electron_]
   Float_t         Electron_hasMatchedConv[kMaxElectron];   //[Electron_]
   Float_t         Electron_mva[kMaxElectron];   //[Electron_]
   Float_t         Electron_mvaPOGTrig[kMaxElectron];   //[Electron_]
   Float_t         Electron_mvaPOGNonTrig[kMaxElectron];   //[Electron_]
   Bool_t          Electron_mvaPreselection[kMaxElectron];   //[Electron_]
   Bool_t          Electron_isTriggerElectron[kMaxElectron];   //[Electron_]
   Float_t         Electron_isoMVA[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIso03v1[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIso03v2[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB03v1[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB03v2[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB03v3[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIso04v1[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIso04v2[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB04v1[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB04v2[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB04v3[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIso03[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIso04[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB03[kMaxElectron];   //[Electron_]
   Float_t         Electron_pfRelIsoDB04[kMaxElectron];   //[Electron_]
   Int_t           Electron_selbit[kMaxElectron];   //[Electron_]
   Int_t           Electron_fidFlag[kMaxElectron];   //[Electron_]
   Int_t           nElectron;
   Int_t           GenParticle_;
   UInt_t          GenParticle_fUniqueID[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_fBits[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_index[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_eta[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_phi[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_p[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_px[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_py[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_pz[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_pt[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_energy[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_pdgId[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_vx[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_vy[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_vz[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_status[kMaxGenParticle];   //[GenParticle_]
   Double_t        GenParticle_charge[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_numDaught[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_numMother[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_motherIndex_1[kMaxGenParticle];   //[GenParticle_]
   UInt_t          GenParticle_motherIndex_2[kMaxGenParticle];   //[GenParticle_]
   Int_t           nGenParticle;
   Int_t           Jet_;
   UInt_t          Jet_fUniqueID[kMaxJet];   //[Jet_]
   UInt_t          Jet_fBits[kMaxJet];   //[Jet_]
   Double_t        Jet_eta[kMaxJet];   //[Jet_]
   Double_t        Jet_phi[kMaxJet];   //[Jet_]
   Double_t        Jet_pt[kMaxJet];   //[Jet_]
   Double_t        Jet_pt_raw[kMaxJet];   //[Jet_]
   Double_t        Jet_energy[kMaxJet];   //[Jet_]
   Double_t        Jet_energy_raw[kMaxJet];   //[Jet_]
   Double_t        Jet_jecUnc[kMaxJet];   //[Jet_]
   Double_t        Jet_resJEC[kMaxJet];   //[Jet_]
   Int_t           Jet_partonFlavour[kMaxJet];   //[Jet_]
   Float_t         Jet_puIdMVA[kMaxJet];   //[Jet_]
   Int_t           Jet_puIdFlag[kMaxJet];   //[Jet_]
   Int_t           Jet_puIdBits[kMaxJet];   //[Jet_]
   Double_t        Jet_chargedEmEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_chargedHadronEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_chargedMuEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_electronEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_muonEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_neutralEmEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_neutralHadronEnergyFraction[kMaxJet];   //[Jet_]
   Double_t        Jet_photonEnergyFraction[kMaxJet];   //[Jet_]
   Int_t           Jet_chargedHadronMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_chargedMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_electronMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_muonMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_neutralHadronMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_neutralMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_photonMultiplicity[kMaxJet];   //[Jet_]
   Int_t           Jet_nConstituents[kMaxJet];   //[Jet_]
   Double_t        Jet_trackCountingHighEffBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_trackCountingHighPurBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_simpleSecondaryVertexHighEffBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_simpleSecondaryVertexHighPurBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_jetProbabilityBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_jetBProbabilityBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_combinedSecondaryVertexBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_combinedSecondaryVertexMVABTag[kMaxJet];   //[Jet_]
   Double_t        Jet_combinedInclusiveSecondaryVertexBTag[kMaxJet];   //[Jet_]
   Double_t        Jet_combinedMVABTag[kMaxJet];   //[Jet_]
   Int_t           Jet_passLooseID[kMaxJet];   //[Jet_]
   Int_t           Jet_passTightID[kMaxJet];   //[Jet_]
   Int_t           Jet_selbit[kMaxJet];   //[Jet_]
   Int_t           nJet;
   Int_t           MET_;
   UInt_t          MET_fUniqueID[kMaxMET];   //[MET_]
   UInt_t          MET_fBits[kMaxMET];   //[MET_]
   Double_t        MET_met[kMaxMET];   //[MET_]
   Double_t        MET_metphi[kMaxMET];   //[MET_]
   Double_t        MET_sumet[kMaxMET];   //[MET_]
   Double_t        MET_metuncorr[kMaxMET];   //[MET_]
   Double_t        MET_metphiuncorr[kMaxMET];   //[MET_]
   Double_t        MET_sumetuncorr[kMaxMET];   //[MET_]
   Int_t           nMET;
   Int_t           corrMET_;
   UInt_t          corrMET_fUniqueID[kMaxcorrMET];   //[corrMET_]
   UInt_t          corrMET_fBits[kMaxcorrMET];   //[corrMET_]
   Double_t        corrMET_met[kMaxcorrMET];   //[corrMET_]
   Double_t        corrMET_metphi[kMaxcorrMET];   //[corrMET_]
   Double_t        corrMET_sumet[kMaxcorrMET];   //[corrMET_]
   Double_t        corrMET_metuncorr[kMaxcorrMET];   //[corrMET_]
   Double_t        corrMET_metphiuncorr[kMaxcorrMET];   //[corrMET_]
   Double_t        corrMET_sumetuncorr[kMaxcorrMET];   //[corrMET_]
   Int_t           corrnMET;
   Int_t           mvaMET_;
   UInt_t          mvaMET_fUniqueID[kMaxmvaMET];   //[mvaMET_]
   UInt_t          mvaMET_fBits[kMaxmvaMET];   //[mvaMET_]
   Double_t        mvaMET_met[kMaxmvaMET];   //[mvaMET_]
   Double_t        mvaMET_metphi[kMaxmvaMET];   //[mvaMET_]
   Double_t        mvaMET_sumet[kMaxmvaMET];   //[mvaMET_]
   Double_t        mvaMET_metuncorr[kMaxmvaMET];   //[mvaMET_]
   Double_t        mvaMET_metphiuncorr[kMaxmvaMET];   //[mvaMET_]
   Double_t        mvaMET_sumetuncorr[kMaxmvaMET];   //[mvaMET_]
   Int_t           mvanMET;
   Int_t           Muon_;
   UInt_t          Muon_fUniqueID[kMaxMuon];   //[Muon_]
   UInt_t          Muon_fBits[kMaxMuon];   //[Muon_]
   Bool_t          Muon_isTrackerMuon[kMaxMuon];   //[Muon_]
   Bool_t          Muon_isPFMuon[kMaxMuon];   //[Muon_]
   Double_t        Muon_eta[kMaxMuon];   //[Muon_]
   Double_t        Muon_phi[kMaxMuon];   //[Muon_]
   Double_t        Muon_pt[kMaxMuon];   //[Muon_]
   Double_t        Muon_ptError[kMaxMuon];   //[Muon_]
   Double_t        Muon_p[kMaxMuon];   //[Muon_]
   Double_t        Muon_energy[kMaxMuon];   //[Muon_]
   Int_t           Muon_charge[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkD0[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkD0Error[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkDz[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkDzError[kMaxMuon];   //[Muon_]
   Double_t        Muon_globalChi2[kMaxMuon];   //[Muon_]
   Double_t        Muon_trkIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_ecalIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_hcalIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_hoIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_relIso[kMaxMuon];   //[Muon_]
   Int_t           Muon_passID[kMaxMuon];   //[Muon_]
   Double_t        Muon_vtxDist3D[kMaxMuon];   //[Muon_]
   Int_t           Muon_vtxIndex[kMaxMuon];   //[Muon_]
   Double_t        Muon_vtxDistZ[kMaxMuon];   //[Muon_]
   Int_t           Muon_pixHits[kMaxMuon];   //[Muon_]
   Int_t           Muon_trkHits[kMaxMuon];   //[Muon_]
   Int_t           Muon_muoHits[kMaxMuon];   //[Muon_]
   Int_t           Muon_matches[kMaxMuon];   //[Muon_]
   Int_t           Muon_trackerLayersWithMeasurement[kMaxMuon];   //[Muon_]
   Double_t        Muon_pfRelIso[kMaxMuon];   //[Muon_]
   Double_t        Muon_vx[kMaxMuon];   //[Muon_]
   Double_t        Muon_vy[kMaxMuon];   //[Muon_]
   Double_t        Muon_vz[kMaxMuon];   //[Muon_]
   Double_t        Muon_dB[kMaxMuon];   //[Muon_]
   Double_t        Muon_edB[kMaxMuon];   //[Muon_]
   Double_t        Muon_dB3d[kMaxMuon];   //[Muon_]
   Double_t        Muon_edB3d[kMaxMuon];   //[Muon_]
   Bool_t          Muon_isGlobalMuonPromptTight[kMaxMuon];   //[Muon_]
   Bool_t          Muon_isAllArbitrated[kMaxMuon];   //[Muon_]
   Int_t           Muon_nChambers[kMaxMuon];   //[Muon_]
   Int_t           Muon_nMatches[kMaxMuon];   //[Muon_]
   Int_t           Muon_nMatchedStations[kMaxMuon];   //[Muon_]
   UInt_t          Muon_stationMask[kMaxMuon];   //[Muon_]
   UInt_t          Muon_stationGapMaskDistance[kMaxMuon];   //[Muon_]
   UInt_t          Muon_stationGapMaskPull[kMaxMuon];   //[Muon_]
   Int_t           Muon_muonID[kMaxMuon];   //[Muon_]
   Float_t         Muon_idMVA[kMaxMuon];   //[Muon_]
   Float_t         Muon_isoRingsMVA[kMaxMuon];   //[Muon_]
   Float_t         Muon_isoRingsRadMVA[kMaxMuon];   //[Muon_]
   Float_t         Muon_idIsoCombMVA[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIso03v1[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIso03v2[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIsoDB03v1[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIsoDB03v2[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIso04v1[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIso04v2[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIsoDB04v1[kMaxMuon];   //[Muon_]
   Float_t         Muon_pfRelIsoDB04v2[kMaxMuon];   //[Muon_]
   Int_t           Muon_selbit[kMaxMuon];   //[Muon_]
   Int_t           nMuon;
   Int_t           Tau_;
   UInt_t          Tau_fUniqueID[kMaxTau];   //[Tau_]
   UInt_t          Tau_fBits[kMaxTau];   //[Tau_]
   Double_t        Tau_eta[kMaxTau];   //[Tau_]
   Double_t        Tau_phi[kMaxTau];   //[Tau_]
   Double_t        Tau_pt[kMaxTau];   //[Tau_]
   Double_t        Tau_energy[kMaxTau];   //[Tau_]
   Int_t           Tau_charge[kMaxTau];   //[Tau_]
   Double_t        Tau_mass[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkPt[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkPtError[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkEta[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkPhi[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkCharge[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkD0[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkD0Error[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkDz[kMaxTau];   //[Tau_]
   Double_t        Tau_leadTrkDzError[kMaxTau];   //[Tau_]
   Int_t           Tau_vtxIndex[kMaxTau];   //[Tau_]
   Double_t        Tau_vtxDxy[kMaxTau];   //[Tau_]
   Double_t        Tau_vtxDz[kMaxTau];   //[Tau_]
   Double_t        Tau_leadChargedParticlePt[kMaxTau];   //[Tau_]
   Double_t        Tau_leadNeutralParticlePt[kMaxTau];   //[Tau_]
   Double_t        Tau_leadParticlePt[kMaxTau];   //[Tau_]
   Int_t           Tau_numChargedHadronsSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numNeutralHadronsSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numPhotonsSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numParticlesSignalCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numChargedHadronsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numNeutralHadronsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numPhotonsIsoCone[kMaxTau];   //[Tau_]
   Int_t           Tau_numParticlesIsoCone[kMaxTau];   //[Tau_]
   Double_t        Tau_ptSumPFChargedHadronsIsoCone[kMaxTau];   //[Tau_]
   Double_t        Tau_ptSumPFNeutralHadronsIsoCone[kMaxTau];   //[Tau_]
   Double_t        Tau_ptSumPhotonsIsoCone[kMaxTau];   //[Tau_]
   Double_t        Tau_sigChHadCandPt[kMaxTau][40];   //[Tau_]
   Double_t        Tau_sigChHadCandEta[kMaxTau][40];   //[Tau_]
   Double_t        Tau_sigChHadCandPhi[kMaxTau][40];   //[Tau_]
   Double_t        Tau_sigNeHadCandPt[kMaxTau][20];   //[Tau_]
   Double_t        Tau_sigNeHadCandEta[kMaxTau][20];   //[Tau_]
   Double_t        Tau_sigNeHadCandPhi[kMaxTau][20];   //[Tau_]
   Double_t        Tau_sigGammaCandPt[kMaxTau][20];   //[Tau_]
   Double_t        Tau_sigGammaCandEta[kMaxTau][20];   //[Tau_]
   Double_t        Tau_sigGammaCandPhi[kMaxTau][20];   //[Tau_]
   Double_t        Tau_isoChHadCandPt[kMaxTau][40];   //[Tau_]
   Double_t        Tau_isoChHadCandEta[kMaxTau][40];   //[Tau_]
   Double_t        Tau_isoChHadCandPhi[kMaxTau][40];   //[Tau_]
   Double_t        Tau_isoNeHadCandPt[kMaxTau][20];   //[Tau_]
   Double_t        Tau_isoNeHadCandEta[kMaxTau][20];   //[Tau_]
   Double_t        Tau_isoNeHadCandPhi[kMaxTau][20];   //[Tau_]
   Double_t        Tau_isoGammaCandPt[kMaxTau][20];   //[Tau_]
   Double_t        Tau_isoGammaCandEta[kMaxTau][20];   //[Tau_]
   Double_t        Tau_isoGammaCandPhi[kMaxTau][20];   //[Tau_]
   Float_t         Tau_decayModeFinding[kMaxTau];   //[Tau_]
   Float_t         Tau_decayModeFindingNewDMs[kMaxTau];   //[Tau_]
   Float_t         Tau_decayModeFindingOldDMs[kMaxTau];   //[Tau_]
   Float_t         Tau_chargedIsoPtSum[kMaxTau];   //[Tau_]
   Float_t         Tau_neutralIsoPtSum[kMaxTau];   //[Tau_]
   Float_t         Tau_puCorrPtSum[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPt3Prong_1track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPt3Prong_2track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPt3Prong_3track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandEta3Prong_1track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandEta3Prong_2track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandEta3Prong_3track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPhi3Prong_1track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPhi3Prong_2track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPhi3Prong_3track[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPt1Prong[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandEta1Prong[kMaxTau];   //[Tau_]
   Double_t        Tau_ChHadCandPhi1Prong[kMaxTau];   //[Tau_]
   Float_t         Tau_CombinedIsolationDeltaBetaCorrRaw3Hits[kMaxTau];   //[Tau_]
   Float_t         Tau_CombinedIsolationDeltaBetaCorrRaw[kMaxTau];   //[Tau_]
   Float_t         Tau_byIsolationMVA3newDMwLTraw[kMaxTau];   //[Tau_]
   Float_t         Tau_byIsolationMVA3newDMwoLTraw[kMaxTau];   //[Tau_]
   Float_t         Tau_byIsolationMVA3oldDMwLTraw[kMaxTau];   //[Tau_]
   Float_t         Tau_byIsolationMVA3oldDMwoLTraw[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseCombinedIsolationDeltaBetaCorr[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumCombinedIsolationDeltaBetaCorr[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightCombinedIsolationDeltaBetaCorr[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightCombinedIsolationDeltaBetaCorr3Hits[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseIsolationMVA3newDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseIsolationMVA3newDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseIsolationMva3oldDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseIsolationMVA3oldDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumIsolationMVA3newDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumIsolationMVA3newDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumIsolationMva3oldDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumIsolationMVA3oldDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightIsolationMVA3newDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightIsolationMVA3newDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightIsolationMva3oldDMwLT[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightIsolationMVA3oldDMwoLT[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonLoose[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonLoose2[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonLoose3[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonLooseMVA[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonTight[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonTight2[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonTight3[kMaxTau];   //[Tau_]
   Float_t         Tau_againstMuonTightMVA[kMaxTau];   //[Tau_]
   Float_t         Tau_againstElectronLoose[kMaxTau];   //[Tau_]
   Float_t         Tau_againstElectronMedium[kMaxTau];   //[Tau_]
   Float_t         Tau_againstElectronTight[kMaxTau];   //[Tau_]
   Float_t         Tau_pfElectronMVA[kMaxTau];   //[Tau_]
   Float_t         Tau_againstElectronMVALooseMVA5[kMaxTau];   //[Tau_]
   Float_t         Tau_againstElectronMVAMediumMVA5[kMaxTau];   //[Tau_]
   Float_t         Tau_againstElectronMVATightMVA5[kMaxTau];   //[Tau_]
   Float_t         Tau_byVLooseCombinedIsolationDeltaBetaCorr[kMaxTau];   //[Tau_]
   Float_t         Tau_byLooseIsolationMVA[kMaxTau];   //[Tau_]
   Float_t         Tau_byMediumIsolationMVA[kMaxTau];   //[Tau_]
   Float_t         Tau_byTightIsolationMVA[kMaxTau];   //[Tau_]
   Double_t        Tau_jetPt[kMaxTau];   //[Tau_]
   Double_t        Tau_jetEta[kMaxTau];   //[Tau_]
   Double_t        Tau_jetPhi[kMaxTau];   //[Tau_]
   Float_t         Tau_emFraction[kMaxTau];   //[Tau_]
   Float_t         Tau_maximumHCALPFClusterEt[kMaxTau];   //[Tau_]
   Float_t         Tau_ecalStripSumEOverPLead[kMaxTau];   //[Tau_]
   Float_t         Tau_bremsRecoveryEOverPLead[kMaxTau];   //[Tau_]
   Float_t         Tau_hcalTotOverPLead[kMaxTau];   //[Tau_]
   Float_t         Tau_hcalMaxOverPLead[kMaxTau];   //[Tau_]
   Float_t         Tau_hcal3x3OverPLead[kMaxTau];   //[Tau_]
   Float_t         Tau_etaetaMoment[kMaxTau];   //[Tau_]
   Float_t         Tau_phiphiMoment[kMaxTau];   //[Tau_]
   Float_t         Tau_etaphiMoment[kMaxTau];   //[Tau_]
   Double_t        Tau_vx[kMaxTau];   //[Tau_]
   Double_t        Tau_vy[kMaxTau];   //[Tau_]
   Double_t        Tau_vz[kMaxTau];   //[Tau_]
   Double_t        Tau_zvertex[kMaxTau];   //[Tau_]
   Double_t        Tau_ltsipt[kMaxTau];   //[Tau_]
   Double_t        Tau_TauVertexPos_x[kMaxTau];   //[Tau_]
   Double_t        Tau_TauVertexPos_y[kMaxTau];   //[Tau_]
   Double_t        Tau_TauVertexPos_z[kMaxTau];   //[Tau_]
   Double_t        Tau_PVPos_x[kMaxTau];   //[Tau_]
   Double_t        Tau_PVPos_y[kMaxTau];   //[Tau_]
   Double_t        Tau_PVPos_z[kMaxTau];   //[Tau_]
   Double_t        Tau_FlightLengh[kMaxTau];   //[Tau_]
   Double_t        Tau_FlightLenghtSig[kMaxTau];   //[Tau_]
   Double_t        Tau_DXY[kMaxTau];   //[Tau_]
   Double_t        Tau_DXY_err[kMaxTau];   //[Tau_]
   Double_t        Tau_DXY_sig[kMaxTau];   //[Tau_]
   Int_t           Tau_NumChHad[kMaxTau];   //[Tau_]
   Int_t           Tau_selbit[kMaxTau];   //[Tau_]
   Int_t           nTau;
   vector<int>     *l1physbits;
   vector<int>     *l1techbits;
   vector<string>  *hltpaths;
   vector<int>     *hltresults;
   vector<int>     *hltprescales;

   // List of branches
   TBranch        *b_Event_;   //!
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_run;   //!
   TBranch        *b_Event_event;   //!
   TBranch        *b_Event_lumis;   //!
   TBranch        *b_Event_bunch;   //!
   TBranch        *b_Event_orbit;   //!
   TBranch        *b_Event_time;   //!
   TBranch        *b_Event_isdata;   //!
   TBranch        *b_Event_isPhysDeclared;   //!
   TBranch        *b_Event_isBPTX0;   //!
   TBranch        *b_Event_isBSCMinBias;   //!
   TBranch        *b_Event_isBSCBeamHalo;   //!
   TBranch        *b_Event_isPrimaryVertex;   //!
   TBranch        *b_Event_isBeamScraping;   //!
   TBranch        *b_Event_passHBHENoiseFilter;   //!
   TBranch        *b_Event_nPU;   //!
   TBranch        *b_Event_bunchCrossing;   //!
   TBranch        *b_Event_trueNInt;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_bunchCrossing;   //!
   TBranch        *b_trueNInt;   //!
   TBranch        *b_Vertex_;   //!
   TBranch        *b_Vertex_fUniqueID;   //!
   TBranch        *b_Vertex_fBits;   //!
   TBranch        *b_Vertex_x;   //!
   TBranch        *b_Vertex_y;   //!
   TBranch        *b_Vertex_z;   //!
   TBranch        *b_Vertex_xErr;   //!
   TBranch        *b_Vertex_yErr;   //!
   TBranch        *b_Vertex_zErr;   //!
   TBranch        *b_Vertex_rho;   //!
   TBranch        *b_Vertex_chi2;   //!
   TBranch        *b_Vertex_ndf;   //!
   TBranch        *b_Vertex_ntracks;   //!
   TBranch        *b_Vertex_ntracksw05;   //!
   TBranch        *b_Vertex_isfake;   //!
   TBranch        *b_Vertex_isvalid;   //!
   TBranch        *b_Vertex_sumPt;   //!
   TBranch        *b_Vertex_selbit;   //!
   TBranch        *b_fnVertex;   //!
   TBranch        *b_Electron_;   //!
   TBranch        *b_Electron_fUniqueID;   //!
   TBranch        *b_Electron_fBits;   //!
   TBranch        *b_Electron_eta;   //!
   TBranch        *b_Electron_phi;   //!
   TBranch        *b_Electron_pt;   //!
   TBranch        *b_Electron_ecalDriven;   //!
   TBranch        *b_Electron_hasGsfTrack;   //!
   TBranch        *b_Electron_trackPt;   //!
   TBranch        *b_Electron_trackPtError;   //!
   TBranch        *b_Electron_energy;   //!
   TBranch        *b_Electron_caloEnergy;   //!
   TBranch        *b_Electron_caloEnergyError;   //!
   TBranch        *b_Electron_charge;   //!
   TBranch        *b_Electron_pixHits;   //!
   TBranch        *b_Electron_trkHits;   //!
   TBranch        *b_Electron_nValidHits;   //!
   TBranch        *b_Electron_trkD0;   //!
   TBranch        *b_Electron_trkD0Error;   //!
   TBranch        *b_Electron_hoe;   //!
   TBranch        *b_Electron_hoeDepth1;   //!
   TBranch        *b_Electron_eop;   //!
   TBranch        *b_Electron_sigmaEtaEta;   //!
   TBranch        *b_Electron_sigmaIEtaIEta;   //!
   TBranch        *b_Electron_deltaPhiTrkSC;   //!
   TBranch        *b_Electron_deltaEtaTrkSC;   //!
   TBranch        *b_Electron_classif;   //!
   TBranch        *b_Electron_e1x5overe5x5;   //!
   TBranch        *b_Electron_e2x5overe5x5;   //!
   TBranch        *b_Electron_isoEcal03;   //!
   TBranch        *b_Electron_isoHcal03;   //!
   TBranch        *b_Electron_isoTrk03;   //!
   TBranch        *b_Electron_isoEcal04;   //!
   TBranch        *b_Electron_isoHcal04;   //!
   TBranch        *b_Electron_isoTrk04;   //!
   TBranch        *b_Electron_isoRel03;   //!
   TBranch        *b_Electron_isoRel04;   //!
   TBranch        *b_Electron_vx;   //!
   TBranch        *b_Electron_vy;   //!
   TBranch        *b_Electron_vz;   //!
   TBranch        *b_Electron_scEn;   //!
   TBranch        *b_Electron_scEta;   //!
   TBranch        *b_Electron_scPhi;   //!
   TBranch        *b_Electron_scET;   //!
   TBranch        *b_Electron_scRawEnergy;   //!
   TBranch        *b_Electron_vtxDist3D;   //!
   TBranch        *b_Electron_vtxIndex;   //!
   TBranch        *b_Electron_vtxDistZ;   //!
   TBranch        *b_Electron_relIso;   //!
   TBranch        *b_Electron_pfRelIso;   //!
   TBranch        *b_Electron_chargedHadronIso;   //!
   TBranch        *b_Electron_neutralHadronIso;   //!
   TBranch        *b_Electron_photonIso;   //!
   TBranch        *b_Electron_missingHits;   //!
   TBranch        *b_Electron_dB;   //!
   TBranch        *b_Electron_edB;   //!
   TBranch        *b_Electron_dB3d;   //!
   TBranch        *b_Electron_edB3d;   //!
   TBranch        *b_Electron_scE1E9;   //!
   TBranch        *b_Electron_scS4S1;   //!
   TBranch        *b_Electron_sckOutOfTime;   //!
   TBranch        *b_Electron_scEcalIso;   //!
   TBranch        *b_Electron_scHEEPEcalIso;   //!
   TBranch        *b_Electron_scHEEPTrkIso;   //!
   TBranch        *b_Electron_nBrems;   //!
   TBranch        *b_Electron_fbrem;   //!
   TBranch        *b_Electron_dist_vec;   //!
   TBranch        *b_Electron_dCotTheta;   //!
   TBranch        *b_Electron_hasMatchedConv;   //!
   TBranch        *b_Electron_mva;   //!
   TBranch        *b_Electron_mvaPOGTrig;   //!
   TBranch        *b_Electron_mvaPOGNonTrig;   //!
   TBranch        *b_Electron_mvaPreselection;   //!
   TBranch        *b_Electron_isTriggerElectron;   //!
   TBranch        *b_Electron_isoMVA;   //!
   TBranch        *b_Electron_pfRelIso03v1;   //!
   TBranch        *b_Electron_pfRelIso03v2;   //!
   TBranch        *b_Electron_pfRelIsoDB03v1;   //!
   TBranch        *b_Electron_pfRelIsoDB03v2;   //!
   TBranch        *b_Electron_pfRelIsoDB03v3;   //!
   TBranch        *b_Electron_pfRelIso04v1;   //!
   TBranch        *b_Electron_pfRelIso04v2;   //!
   TBranch        *b_Electron_pfRelIsoDB04v1;   //!
   TBranch        *b_Electron_pfRelIsoDB04v2;   //!
   TBranch        *b_Electron_pfRelIsoDB04v3;   //!
   TBranch        *b_Electron_pfRelIso03;   //!
   TBranch        *b_Electron_pfRelIso04;   //!
   TBranch        *b_Electron_pfRelIsoDB03;   //!
   TBranch        *b_Electron_pfRelIsoDB04;   //!
   TBranch        *b_Electron_selbit;   //!
   TBranch        *b_Electron_fidFlag;   //!
   TBranch        *b_fnElectron;   //!
   TBranch        *b_GenParticle_;   //!
   TBranch        *b_GenParticle_fUniqueID;   //!
   TBranch        *b_GenParticle_fBits;   //!
   TBranch        *b_GenParticle_index;   //!
   TBranch        *b_GenParticle_eta;   //!
   TBranch        *b_GenParticle_phi;   //!
   TBranch        *b_GenParticle_p;   //!
   TBranch        *b_GenParticle_px;   //!
   TBranch        *b_GenParticle_py;   //!
   TBranch        *b_GenParticle_pz;   //!
   TBranch        *b_GenParticle_pt;   //!
   TBranch        *b_GenParticle_energy;   //!
   TBranch        *b_GenParticle_pdgId;   //!
   TBranch        *b_GenParticle_vx;   //!
   TBranch        *b_GenParticle_vy;   //!
   TBranch        *b_GenParticle_vz;   //!
   TBranch        *b_GenParticle_status;   //!
   TBranch        *b_GenParticle_charge;   //!
   TBranch        *b_GenParticle_numDaught;   //!
   TBranch        *b_GenParticle_numMother;   //!
   TBranch        *b_GenParticle_motherIndex_1;   //!
   TBranch        *b_GenParticle_motherIndex_2;   //!
   TBranch        *b_fnGenParticle;   //!
   TBranch        *b_Jet_;   //!
   TBranch        *b_Jet_fUniqueID;   //!
   TBranch        *b_Jet_fBits;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_pt_raw;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_energy_raw;   //!
   TBranch        *b_Jet_jecUnc;   //!
   TBranch        *b_Jet_resJEC;   //!
   TBranch        *b_Jet_partonFlavour;   //!
   TBranch        *b_Jet_puIdMVA;   //!
   TBranch        *b_Jet_puIdFlag;   //!
   TBranch        *b_Jet_puIdBits;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedMuEnergyFraction;   //!
   TBranch        *b_Jet_electronEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEnergyFraction;   //!
   TBranch        *b_Jet_neutralHadronEnergyFraction;   //!
   TBranch        *b_Jet_photonEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronMultiplicity;   //!
   TBranch        *b_Jet_chargedMultiplicity;   //!
   TBranch        *b_Jet_electronMultiplicity;   //!
   TBranch        *b_Jet_muonMultiplicity;   //!
   TBranch        *b_Jet_neutralHadronMultiplicity;   //!
   TBranch        *b_Jet_neutralMultiplicity;   //!
   TBranch        *b_Jet_photonMultiplicity;   //!
   TBranch        *b_Jet_nConstituents;   //!
   TBranch        *b_Jet_trackCountingHighEffBTag;   //!
   TBranch        *b_Jet_trackCountingHighPurBTag;   //!
   TBranch        *b_Jet_simpleSecondaryVertexHighEffBTag;   //!
   TBranch        *b_Jet_simpleSecondaryVertexHighPurBTag;   //!
   TBranch        *b_Jet_jetProbabilityBTag;   //!
   TBranch        *b_Jet_jetBProbabilityBTag;   //!
   TBranch        *b_Jet_combinedSecondaryVertexBTag;   //!
   TBranch        *b_Jet_combinedSecondaryVertexMVABTag;   //!
   TBranch        *b_Jet_combinedInclusiveSecondaryVertexBTag;   //!
   TBranch        *b_Jet_combinedMVABTag;   //!
   TBranch        *b_Jet_passLooseID;   //!
   TBranch        *b_Jet_passTightID;   //!
   TBranch        *b_Jet_selbit;   //!
   TBranch        *b_fnJet;   //!
   TBranch        *b_MET_;   //!
   TBranch        *b_MET_fUniqueID;   //!
   TBranch        *b_MET_fBits;   //!
   TBranch        *b_MET_met;   //!
   TBranch        *b_MET_metphi;   //!
   TBranch        *b_MET_sumet;   //!
   TBranch        *b_MET_metuncorr;   //!
   TBranch        *b_MET_metphiuncorr;   //!
   TBranch        *b_MET_sumetuncorr;   //!
   TBranch        *b_fnPFMET;   //!
   TBranch        *b_corrMET_;   //!
   TBranch        *b_corrMET_fUniqueID;   //!
   TBranch        *b_corrMET_fBits;   //!
   TBranch        *b_corrMET_met;   //!
   TBranch        *b_corrMET_metphi;   //!
   TBranch        *b_corrMET_sumet;   //!
   TBranch        *b_corrMET_metuncorr;   //!
   TBranch        *b_corrMET_metphiuncorr;   //!
   TBranch        *b_corrMET_sumetuncorr;   //!
   TBranch        *b_fnCorrMET;   //!
   TBranch        *b_mvaMET_;   //!
   TBranch        *b_mvaMET_fUniqueID;   //!
   TBranch        *b_mvaMET_fBits;   //!
   TBranch        *b_mvaMET_met;   //!
   TBranch        *b_mvaMET_metphi;   //!
   TBranch        *b_mvaMET_sumet;   //!
   TBranch        *b_mvaMET_metuncorr;   //!
   TBranch        *b_mvaMET_metphiuncorr;   //!
   TBranch        *b_mvaMET_sumetuncorr;   //!
   TBranch        *b_fnMVAMET;   //!
   TBranch        *b_Muon_;   //!
   TBranch        *b_Muon_fUniqueID;   //!
   TBranch        *b_Muon_fBits;   //!
   TBranch        *b_Muon_isTrackerMuon;   //!
   TBranch        *b_Muon_isPFMuon;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_ptError;   //!
   TBranch        *b_Muon_p;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_trkD0;   //!
   TBranch        *b_Muon_trkD0Error;   //!
   TBranch        *b_Muon_trkDz;   //!
   TBranch        *b_Muon_trkDzError;   //!
   TBranch        *b_Muon_globalChi2;   //!
   TBranch        *b_Muon_trkIso;   //!
   TBranch        *b_Muon_ecalIso;   //!
   TBranch        *b_Muon_hcalIso;   //!
   TBranch        *b_Muon_hoIso;   //!
   TBranch        *b_Muon_relIso;   //!
   TBranch        *b_Muon_passID;   //!
   TBranch        *b_Muon_vtxDist3D;   //!
   TBranch        *b_Muon_vtxIndex;   //!
   TBranch        *b_Muon_vtxDistZ;   //!
   TBranch        *b_Muon_pixHits;   //!
   TBranch        *b_Muon_trkHits;   //!
   TBranch        *b_Muon_muoHits;   //!
   TBranch        *b_Muon_matches;   //!
   TBranch        *b_Muon_trackerLayersWithMeasurement;   //!
   TBranch        *b_Muon_pfRelIso;   //!
   TBranch        *b_Muon_vx;   //!
   TBranch        *b_Muon_vy;   //!
   TBranch        *b_Muon_vz;   //!
   TBranch        *b_Muon_dB;   //!
   TBranch        *b_Muon_edB;   //!
   TBranch        *b_Muon_dB3d;   //!
   TBranch        *b_Muon_edB3d;   //!
   TBranch        *b_Muon_isGlobalMuonPromptTight;   //!
   TBranch        *b_Muon_isAllArbitrated;   //!
   TBranch        *b_Muon_nChambers;   //!
   TBranch        *b_Muon_nMatches;   //!
   TBranch        *b_Muon_nMatchedStations;   //!
   TBranch        *b_Muon_stationMask;   //!
   TBranch        *b_Muon_stationGapMaskDistance;   //!
   TBranch        *b_Muon_stationGapMaskPull;   //!
   TBranch        *b_Muon_muonID;   //!
   TBranch        *b_Muon_idMVA;   //!
   TBranch        *b_Muon_isoRingsMVA;   //!
   TBranch        *b_Muon_isoRingsRadMVA;   //!
   TBranch        *b_Muon_idIsoCombMVA;   //!
   TBranch        *b_Muon_pfRelIso03v1;   //!
   TBranch        *b_Muon_pfRelIso03v2;   //!
   TBranch        *b_Muon_pfRelIsoDB03v1;   //!
   TBranch        *b_Muon_pfRelIsoDB03v2;   //!
   TBranch        *b_Muon_pfRelIso04v1;   //!
   TBranch        *b_Muon_pfRelIso04v2;   //!
   TBranch        *b_Muon_pfRelIsoDB04v1;   //!
   TBranch        *b_Muon_pfRelIsoDB04v2;   //!
   TBranch        *b_Muon_selbit;   //!
   TBranch        *b_fnMuon;   //!
   TBranch        *b_Tau_;   //!
   TBranch        *b_Tau_fUniqueID;   //!
   TBranch        *b_Tau_fBits;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_mass;   //!
   TBranch        *b_Tau_leadTrkPt;   //!
   TBranch        *b_Tau_leadTrkPtError;   //!
   TBranch        *b_Tau_leadTrkEta;   //!
   TBranch        *b_Tau_leadTrkPhi;   //!
   TBranch        *b_Tau_leadTrkCharge;   //!
   TBranch        *b_Tau_leadTrkD0;   //!
   TBranch        *b_Tau_leadTrkD0Error;   //!
   TBranch        *b_Tau_leadTrkDz;   //!
   TBranch        *b_Tau_leadTrkDzError;   //!
   TBranch        *b_Tau_vtxIndex;   //!
   TBranch        *b_Tau_vtxDxy;   //!
   TBranch        *b_Tau_vtxDz;   //!
   TBranch        *b_Tau_leadChargedParticlePt;   //!
   TBranch        *b_Tau_leadNeutralParticlePt;   //!
   TBranch        *b_Tau_leadParticlePt;   //!
   TBranch        *b_Tau_numChargedHadronsSignalCone;   //!
   TBranch        *b_Tau_numNeutralHadronsSignalCone;   //!
   TBranch        *b_Tau_numPhotonsSignalCone;   //!
   TBranch        *b_Tau_numParticlesSignalCone;   //!
   TBranch        *b_Tau_numChargedHadronsIsoCone;   //!
   TBranch        *b_Tau_numNeutralHadronsIsoCone;   //!
   TBranch        *b_Tau_numPhotonsIsoCone;   //!
   TBranch        *b_Tau_numParticlesIsoCone;   //!
   TBranch        *b_Tau_ptSumPFChargedHadronsIsoCone;   //!
   TBranch        *b_Tau_ptSumPFNeutralHadronsIsoCone;   //!
   TBranch        *b_Tau_ptSumPhotonsIsoCone;   //!
   TBranch        *b_Tau_sigChHadCandPt;   //!
   TBranch        *b_Tau_sigChHadCandEta;   //!
   TBranch        *b_Tau_sigChHadCandPhi;   //!
   TBranch        *b_Tau_sigNeHadCandPt;   //!
   TBranch        *b_Tau_sigNeHadCandEta;   //!
   TBranch        *b_Tau_sigNeHadCandPhi;   //!
   TBranch        *b_Tau_sigGammaCandPt;   //!
   TBranch        *b_Tau_sigGammaCandEta;   //!
   TBranch        *b_Tau_sigGammaCandPhi;   //!
   TBranch        *b_Tau_isoChHadCandPt;   //!
   TBranch        *b_Tau_isoChHadCandEta;   //!
   TBranch        *b_Tau_isoChHadCandPhi;   //!
   TBranch        *b_Tau_isoNeHadCandPt;   //!
   TBranch        *b_Tau_isoNeHadCandEta;   //!
   TBranch        *b_Tau_isoNeHadCandPhi;   //!
   TBranch        *b_Tau_isoGammaCandPt;   //!
   TBranch        *b_Tau_isoGammaCandEta;   //!
   TBranch        *b_Tau_isoGammaCandPhi;   //!
   TBranch        *b_Tau_decayModeFinding;   //!
   TBranch        *b_Tau_decayModeFindingNewDMs;   //!
   TBranch        *b_Tau_decayModeFindingOldDMs;   //!
   TBranch        *b_Tau_chargedIsoPtSum;   //!
   TBranch        *b_Tau_neutralIsoPtSum;   //!
   TBranch        *b_Tau_puCorrPtSum;   //!
   TBranch        *b_Tau_ChHadCandPt3Prong_1track;   //!
   TBranch        *b_Tau_ChHadCandPt3Prong_2track;   //!
   TBranch        *b_Tau_ChHadCandPt3Prong_3track;   //!
   TBranch        *b_Tau_ChHadCandEta3Prong_1track;   //!
   TBranch        *b_Tau_ChHadCandEta3Prong_2track;   //!
   TBranch        *b_Tau_ChHadCandEta3Prong_3track;   //!
   TBranch        *b_Tau_ChHadCandPhi3Prong_1track;   //!
   TBranch        *b_Tau_ChHadCandPhi3Prong_2track;   //!
   TBranch        *b_Tau_ChHadCandPhi3Prong_3track;   //!
   TBranch        *b_Tau_ChHadCandPt1Prong;   //!
   TBranch        *b_Tau_ChHadCandEta1Prong;   //!
   TBranch        *b_Tau_ChHadCandPhi1Prong;   //!
   TBranch        *b_Tau_CombinedIsolationDeltaBetaCorrRaw3Hits;   //!
   TBranch        *b_Tau_CombinedIsolationDeltaBetaCorrRaw;   //!
   TBranch        *b_Tau_byIsolationMVA3newDMwLTraw;   //!
   TBranch        *b_Tau_byIsolationMVA3newDMwoLTraw;   //!
   TBranch        *b_Tau_byIsolationMVA3oldDMwLTraw;   //!
   TBranch        *b_Tau_byIsolationMVA3oldDMwoLTraw;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byLooseIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byMediumIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byTightIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_againstMuonLoose;   //!
   TBranch        *b_Tau_againstMuonLoose2;   //!
   TBranch        *b_Tau_againstMuonLoose3;   //!
   TBranch        *b_Tau_againstMuonLooseMVA;   //!
   TBranch        *b_Tau_againstMuonTight;   //!
   TBranch        *b_Tau_againstMuonTight2;   //!
   TBranch        *b_Tau_againstMuonTight3;   //!
   TBranch        *b_Tau_againstMuonTightMVA;   //!
   TBranch        *b_Tau_againstElectronLoose;   //!
   TBranch        *b_Tau_againstElectronMedium;   //!
   TBranch        *b_Tau_againstElectronTight;   //!
   TBranch        *b_Tau_pfElectronMVA;   //!
   TBranch        *b_Tau_againstElectronMVALooseMVA5;   //!
   TBranch        *b_Tau_againstElectronMVAMediumMVA5;   //!
   TBranch        *b_Tau_againstElectronMVATightMVA5;   //!
   TBranch        *b_Tau_byVLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byLooseIsolationMVA;   //!
   TBranch        *b_Tau_byMediumIsolationMVA;   //!
   TBranch        *b_Tau_byTightIsolationMVA;   //!
   TBranch        *b_Tau_jetPt;   //!
   TBranch        *b_Tau_jetEta;   //!
   TBranch        *b_Tau_jetPhi;   //!
   TBranch        *b_Tau_emFraction;   //!
   TBranch        *b_Tau_maximumHCALPFClusterEt;   //!
   TBranch        *b_Tau_ecalStripSumEOverPLead;   //!
   TBranch        *b_Tau_bremsRecoveryEOverPLead;   //!
   TBranch        *b_Tau_hcalTotOverPLead;   //!
   TBranch        *b_Tau_hcalMaxOverPLead;   //!
   TBranch        *b_Tau_hcal3x3OverPLead;   //!
   TBranch        *b_Tau_etaetaMoment;   //!
   TBranch        *b_Tau_phiphiMoment;   //!
   TBranch        *b_Tau_etaphiMoment;   //!
   TBranch        *b_Tau_vx;   //!
   TBranch        *b_Tau_vy;   //!
   TBranch        *b_Tau_vz;   //!
   TBranch        *b_Tau_zvertex;   //!
   TBranch        *b_Tau_ltsipt;   //!
   TBranch        *b_Tau_TauVertexPos_x;   //!
   TBranch        *b_Tau_TauVertexPos_y;   //!
   TBranch        *b_Tau_TauVertexPos_z;   //!
   TBranch        *b_Tau_PVPos_x;   //!
   TBranch        *b_Tau_PVPos_y;   //!
   TBranch        *b_Tau_PVPos_z;   //!
   TBranch        *b_Tau_FlightLengh;   //!
   TBranch        *b_Tau_FlightLenghtSig;   //!
   TBranch        *b_Tau_DXY;   //!
   TBranch        *b_Tau_DXY_err;   //!
   TBranch        *b_Tau_DXY_sig;   //!
   TBranch        *b_Tau_NumChHad;   //!
   TBranch        *b_Tau_selbit;   //!
   TBranch        *b_fnTau;   //!
   TBranch        *b_l1physbits;   //!
   TBranch        *b_l1techbits;   //!
   TBranch        *b_hltpaths;   //!
   TBranch        *b_hltresults;   //!
   TBranch        *b_hltprescales;   //!

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/data/mc/analysis/Tree_RadionM500.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/data/mc/analysis/Tree_RadionM500.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/data/mc/analysis/Tree_RadionM500.root:/treeCreator");
      dir->GetObject("vhtree",tree);

   }
   Init(tree);
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   nPU = 0;
   bunchCrossing = 0;
   trueNInt = 0;
   l1physbits = 0;
   l1techbits = 0;
   hltpaths = 0;
   hltresults = 0;
   hltprescales = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event_, &b_Event_);
   fChain->SetBranchAddress("Event.fUniqueID", Event_fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("Event.fBits", Event_fBits, &b_Event_fBits);
   fChain->SetBranchAddress("Event.run", Event_run, &b_Event_run);
   fChain->SetBranchAddress("Event.event", Event_event, &b_Event_event);
   fChain->SetBranchAddress("Event.lumis", Event_lumis, &b_Event_lumis);
   fChain->SetBranchAddress("Event.bunch", Event_bunch, &b_Event_bunch);
   fChain->SetBranchAddress("Event.orbit", Event_orbit, &b_Event_orbit);
   fChain->SetBranchAddress("Event.time", Event_time, &b_Event_time);
   fChain->SetBranchAddress("Event.isdata", Event_isdata, &b_Event_isdata);
   fChain->SetBranchAddress("Event.isPhysDeclared", Event_isPhysDeclared, &b_Event_isPhysDeclared);
   fChain->SetBranchAddress("Event.isBPTX0", Event_isBPTX0, &b_Event_isBPTX0);
   fChain->SetBranchAddress("Event.isBSCMinBias", Event_isBSCMinBias, &b_Event_isBSCMinBias);
   fChain->SetBranchAddress("Event.isBSCBeamHalo", Event_isBSCBeamHalo, &b_Event_isBSCBeamHalo);
   fChain->SetBranchAddress("Event.isPrimaryVertex", Event_isPrimaryVertex, &b_Event_isPrimaryVertex);
   fChain->SetBranchAddress("Event.isBeamScraping", Event_isBeamScraping, &b_Event_isBeamScraping);
   fChain->SetBranchAddress("Event.passHBHENoiseFilter", Event_passHBHENoiseFilter, &b_Event_passHBHENoiseFilter);
   fChain->SetBranchAddress("Event.nPU", Event_nPU, &b_Event_nPU);
   fChain->SetBranchAddress("Event.bunchCrossing", Event_bunchCrossing, &b_Event_bunchCrossing);
   fChain->SetBranchAddress("Event.trueNInt", Event_trueNInt, &b_Event_trueNInt);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("bunchCrossing", &bunchCrossing, &b_bunchCrossing);
   fChain->SetBranchAddress("trueNInt", &trueNInt, &b_trueNInt);
   fChain->SetBranchAddress("Vertex", &Vertex_, &b_Vertex_);
   fChain->SetBranchAddress("Vertex.fUniqueID", Vertex_fUniqueID, &b_Vertex_fUniqueID);
   fChain->SetBranchAddress("Vertex.fBits", Vertex_fBits, &b_Vertex_fBits);
   fChain->SetBranchAddress("Vertex.x", Vertex_x, &b_Vertex_x);
   fChain->SetBranchAddress("Vertex.y", Vertex_y, &b_Vertex_y);
   fChain->SetBranchAddress("Vertex.z", Vertex_z, &b_Vertex_z);
   fChain->SetBranchAddress("Vertex.xErr", Vertex_xErr, &b_Vertex_xErr);
   fChain->SetBranchAddress("Vertex.yErr", Vertex_yErr, &b_Vertex_yErr);
   fChain->SetBranchAddress("Vertex.zErr", Vertex_zErr, &b_Vertex_zErr);
   fChain->SetBranchAddress("Vertex.rho", Vertex_rho, &b_Vertex_rho);
   fChain->SetBranchAddress("Vertex.chi2", Vertex_chi2, &b_Vertex_chi2);
   fChain->SetBranchAddress("Vertex.ndf", Vertex_ndf, &b_Vertex_ndf);
   fChain->SetBranchAddress("Vertex.ntracks", Vertex_ntracks, &b_Vertex_ntracks);
   fChain->SetBranchAddress("Vertex.ntracksw05", Vertex_ntracksw05, &b_Vertex_ntracksw05);
   fChain->SetBranchAddress("Vertex.isfake", Vertex_isfake, &b_Vertex_isfake);
   fChain->SetBranchAddress("Vertex.isvalid", Vertex_isvalid, &b_Vertex_isvalid);
   fChain->SetBranchAddress("Vertex.sumPt", Vertex_sumPt, &b_Vertex_sumPt);
   fChain->SetBranchAddress("Vertex.selbit", Vertex_selbit, &b_Vertex_selbit);
   fChain->SetBranchAddress("nVertex", &nVertex, &b_fnVertex);
   fChain->SetBranchAddress("Electron", &Electron_, &b_Electron_);
   fChain->SetBranchAddress("Electron.fUniqueID", Electron_fUniqueID, &b_Electron_fUniqueID);
   fChain->SetBranchAddress("Electron.fBits", Electron_fBits, &b_Electron_fBits);
   fChain->SetBranchAddress("Electron.eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron.phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron.pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron.ecalDriven", Electron_ecalDriven, &b_Electron_ecalDriven);
   fChain->SetBranchAddress("Electron.hasGsfTrack", Electron_hasGsfTrack, &b_Electron_hasGsfTrack);
   fChain->SetBranchAddress("Electron.trackPt", Electron_trackPt, &b_Electron_trackPt);
   fChain->SetBranchAddress("Electron.trackPtError", Electron_trackPtError, &b_Electron_trackPtError);
   fChain->SetBranchAddress("Electron.energy", Electron_energy, &b_Electron_energy);
   fChain->SetBranchAddress("Electron.caloEnergy", Electron_caloEnergy, &b_Electron_caloEnergy);
   fChain->SetBranchAddress("Electron.caloEnergyError", Electron_caloEnergyError, &b_Electron_caloEnergyError);
   fChain->SetBranchAddress("Electron.charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron.pixHits", Electron_pixHits, &b_Electron_pixHits);
   fChain->SetBranchAddress("Electron.trkHits", Electron_trkHits, &b_Electron_trkHits);
   fChain->SetBranchAddress("Electron.nValidHits", Electron_nValidHits, &b_Electron_nValidHits);
   fChain->SetBranchAddress("Electron.trkD0", Electron_trkD0, &b_Electron_trkD0);
   fChain->SetBranchAddress("Electron.trkD0Error", Electron_trkD0Error, &b_Electron_trkD0Error);
   fChain->SetBranchAddress("Electron.hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron.hoeDepth1", Electron_hoeDepth1, &b_Electron_hoeDepth1);
   fChain->SetBranchAddress("Electron.eop", Electron_eop, &b_Electron_eop);
   fChain->SetBranchAddress("Electron.sigmaEtaEta", Electron_sigmaEtaEta, &b_Electron_sigmaEtaEta);
   fChain->SetBranchAddress("Electron.sigmaIEtaIEta", Electron_sigmaIEtaIEta, &b_Electron_sigmaIEtaIEta);
   fChain->SetBranchAddress("Electron.deltaPhiTrkSC", Electron_deltaPhiTrkSC, &b_Electron_deltaPhiTrkSC);
   fChain->SetBranchAddress("Electron.deltaEtaTrkSC", Electron_deltaEtaTrkSC, &b_Electron_deltaEtaTrkSC);
   fChain->SetBranchAddress("Electron.classif", Electron_classif, &b_Electron_classif);
   fChain->SetBranchAddress("Electron.e1x5overe5x5", Electron_e1x5overe5x5, &b_Electron_e1x5overe5x5);
   fChain->SetBranchAddress("Electron.e2x5overe5x5", Electron_e2x5overe5x5, &b_Electron_e2x5overe5x5);
   fChain->SetBranchAddress("Electron.isoEcal03", Electron_isoEcal03, &b_Electron_isoEcal03);
   fChain->SetBranchAddress("Electron.isoHcal03", Electron_isoHcal03, &b_Electron_isoHcal03);
   fChain->SetBranchAddress("Electron.isoTrk03", Electron_isoTrk03, &b_Electron_isoTrk03);
   fChain->SetBranchAddress("Electron.isoEcal04", Electron_isoEcal04, &b_Electron_isoEcal04);
   fChain->SetBranchAddress("Electron.isoHcal04", Electron_isoHcal04, &b_Electron_isoHcal04);
   fChain->SetBranchAddress("Electron.isoTrk04", Electron_isoTrk04, &b_Electron_isoTrk04);
   fChain->SetBranchAddress("Electron.isoRel03", Electron_isoRel03, &b_Electron_isoRel03);
   fChain->SetBranchAddress("Electron.isoRel04", Electron_isoRel04, &b_Electron_isoRel04);
   fChain->SetBranchAddress("Electron.vx", Electron_vx, &b_Electron_vx);
   fChain->SetBranchAddress("Electron.vy", Electron_vy, &b_Electron_vy);
   fChain->SetBranchAddress("Electron.vz", Electron_vz, &b_Electron_vz);
   fChain->SetBranchAddress("Electron.scEn", Electron_scEn, &b_Electron_scEn);
   fChain->SetBranchAddress("Electron.scEta", Electron_scEta, &b_Electron_scEta);
   fChain->SetBranchAddress("Electron.scPhi", Electron_scPhi, &b_Electron_scPhi);
   fChain->SetBranchAddress("Electron.scET", Electron_scET, &b_Electron_scET);
   fChain->SetBranchAddress("Electron.scRawEnergy", Electron_scRawEnergy, &b_Electron_scRawEnergy);
   fChain->SetBranchAddress("Electron.vtxDist3D", Electron_vtxDist3D, &b_Electron_vtxDist3D);
   fChain->SetBranchAddress("Electron.vtxIndex", Electron_vtxIndex, &b_Electron_vtxIndex);
   fChain->SetBranchAddress("Electron.vtxDistZ", Electron_vtxDistZ, &b_Electron_vtxDistZ);
   fChain->SetBranchAddress("Electron.relIso", Electron_relIso, &b_Electron_relIso);
   fChain->SetBranchAddress("Electron.pfRelIso", Electron_pfRelIso, &b_Electron_pfRelIso);
   fChain->SetBranchAddress("Electron.chargedHadronIso", Electron_chargedHadronIso, &b_Electron_chargedHadronIso);
   fChain->SetBranchAddress("Electron.neutralHadronIso", Electron_neutralHadronIso, &b_Electron_neutralHadronIso);
   fChain->SetBranchAddress("Electron.photonIso", Electron_photonIso, &b_Electron_photonIso);
   fChain->SetBranchAddress("Electron.missingHits", Electron_missingHits, &b_Electron_missingHits);
   fChain->SetBranchAddress("Electron.dB", Electron_dB, &b_Electron_dB);
   fChain->SetBranchAddress("Electron.edB", Electron_edB, &b_Electron_edB);
   fChain->SetBranchAddress("Electron.dB3d", Electron_dB3d, &b_Electron_dB3d);
   fChain->SetBranchAddress("Electron.edB3d", Electron_edB3d, &b_Electron_edB3d);
   fChain->SetBranchAddress("Electron.scE1E9", Electron_scE1E9, &b_Electron_scE1E9);
   fChain->SetBranchAddress("Electron.scS4S1", Electron_scS4S1, &b_Electron_scS4S1);
   fChain->SetBranchAddress("Electron.sckOutOfTime", Electron_sckOutOfTime, &b_Electron_sckOutOfTime);
   fChain->SetBranchAddress("Electron.scEcalIso", Electron_scEcalIso, &b_Electron_scEcalIso);
   fChain->SetBranchAddress("Electron.scHEEPEcalIso", Electron_scHEEPEcalIso, &b_Electron_scHEEPEcalIso);
   fChain->SetBranchAddress("Electron.scHEEPTrkIso", Electron_scHEEPTrkIso, &b_Electron_scHEEPTrkIso);
   fChain->SetBranchAddress("Electron.nBrems", Electron_nBrems, &b_Electron_nBrems);
   fChain->SetBranchAddress("Electron.fbrem", Electron_fbrem, &b_Electron_fbrem);
   fChain->SetBranchAddress("Electron.dist_vec", Electron_dist_vec, &b_Electron_dist_vec);
   fChain->SetBranchAddress("Electron.dCotTheta", Electron_dCotTheta, &b_Electron_dCotTheta);
   fChain->SetBranchAddress("Electron.hasMatchedConv", Electron_hasMatchedConv, &b_Electron_hasMatchedConv);
   fChain->SetBranchAddress("Electron.mva", Electron_mva, &b_Electron_mva);
   fChain->SetBranchAddress("Electron.mvaPOGTrig", Electron_mvaPOGTrig, &b_Electron_mvaPOGTrig);
   fChain->SetBranchAddress("Electron.mvaPOGNonTrig", Electron_mvaPOGNonTrig, &b_Electron_mvaPOGNonTrig);
   fChain->SetBranchAddress("Electron.mvaPreselection", Electron_mvaPreselection, &b_Electron_mvaPreselection);
   fChain->SetBranchAddress("Electron.isTriggerElectron", Electron_isTriggerElectron, &b_Electron_isTriggerElectron);
   fChain->SetBranchAddress("Electron.isoMVA", Electron_isoMVA, &b_Electron_isoMVA);
   fChain->SetBranchAddress("Electron.pfRelIso03v1", Electron_pfRelIso03v1, &b_Electron_pfRelIso03v1);
   fChain->SetBranchAddress("Electron.pfRelIso03v2", Electron_pfRelIso03v2, &b_Electron_pfRelIso03v2);
   fChain->SetBranchAddress("Electron.pfRelIsoDB03v1", Electron_pfRelIsoDB03v1, &b_Electron_pfRelIsoDB03v1);
   fChain->SetBranchAddress("Electron.pfRelIsoDB03v2", Electron_pfRelIsoDB03v2, &b_Electron_pfRelIsoDB03v2);
   fChain->SetBranchAddress("Electron.pfRelIsoDB03v3", Electron_pfRelIsoDB03v3, &b_Electron_pfRelIsoDB03v3);
   fChain->SetBranchAddress("Electron.pfRelIso04v1", Electron_pfRelIso04v1, &b_Electron_pfRelIso04v1);
   fChain->SetBranchAddress("Electron.pfRelIso04v2", Electron_pfRelIso04v2, &b_Electron_pfRelIso04v2);
   fChain->SetBranchAddress("Electron.pfRelIsoDB04v1", Electron_pfRelIsoDB04v1, &b_Electron_pfRelIsoDB04v1);
   fChain->SetBranchAddress("Electron.pfRelIsoDB04v2", Electron_pfRelIsoDB04v2, &b_Electron_pfRelIsoDB04v2);
   fChain->SetBranchAddress("Electron.pfRelIsoDB04v3", Electron_pfRelIsoDB04v3, &b_Electron_pfRelIsoDB04v3);
   fChain->SetBranchAddress("Electron.pfRelIso03", Electron_pfRelIso03, &b_Electron_pfRelIso03);
   fChain->SetBranchAddress("Electron.pfRelIso04", Electron_pfRelIso04, &b_Electron_pfRelIso04);
   fChain->SetBranchAddress("Electron.pfRelIsoDB03", Electron_pfRelIsoDB03, &b_Electron_pfRelIsoDB03);
   fChain->SetBranchAddress("Electron.pfRelIsoDB04", Electron_pfRelIsoDB04, &b_Electron_pfRelIsoDB04);
   fChain->SetBranchAddress("Electron.selbit", Electron_selbit, &b_Electron_selbit);
   fChain->SetBranchAddress("Electron.fidFlag", Electron_fidFlag, &b_Electron_fidFlag);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_fnElectron);
   fChain->SetBranchAddress("GenParticle", &GenParticle_, &b_GenParticle_);
   fChain->SetBranchAddress("GenParticle.fUniqueID", GenParticle_fUniqueID, &b_GenParticle_fUniqueID);
   fChain->SetBranchAddress("GenParticle.fBits", GenParticle_fBits, &b_GenParticle_fBits);
   fChain->SetBranchAddress("GenParticle.index", GenParticle_index, &b_GenParticle_index);
   fChain->SetBranchAddress("GenParticle.eta", GenParticle_eta, &b_GenParticle_eta);
   fChain->SetBranchAddress("GenParticle.phi", GenParticle_phi, &b_GenParticle_phi);
   fChain->SetBranchAddress("GenParticle.p", GenParticle_p, &b_GenParticle_p);
   fChain->SetBranchAddress("GenParticle.px", GenParticle_px, &b_GenParticle_px);
   fChain->SetBranchAddress("GenParticle.py", GenParticle_py, &b_GenParticle_py);
   fChain->SetBranchAddress("GenParticle.pz", GenParticle_pz, &b_GenParticle_pz);
   fChain->SetBranchAddress("GenParticle.pt", GenParticle_pt, &b_GenParticle_pt);
   fChain->SetBranchAddress("GenParticle.energy", GenParticle_energy, &b_GenParticle_energy);
   fChain->SetBranchAddress("GenParticle.pdgId", GenParticle_pdgId, &b_GenParticle_pdgId);
   fChain->SetBranchAddress("GenParticle.vx", GenParticle_vx, &b_GenParticle_vx);
   fChain->SetBranchAddress("GenParticle.vy", GenParticle_vy, &b_GenParticle_vy);
   fChain->SetBranchAddress("GenParticle.vz", GenParticle_vz, &b_GenParticle_vz);
   fChain->SetBranchAddress("GenParticle.status", GenParticle_status, &b_GenParticle_status);
   fChain->SetBranchAddress("GenParticle.charge", GenParticle_charge, &b_GenParticle_charge);
   fChain->SetBranchAddress("GenParticle.numDaught", GenParticle_numDaught, &b_GenParticle_numDaught);
   fChain->SetBranchAddress("GenParticle.numMother", GenParticle_numMother, &b_GenParticle_numMother);
   fChain->SetBranchAddress("GenParticle.motherIndex_1", GenParticle_motherIndex_1, &b_GenParticle_motherIndex_1);
   fChain->SetBranchAddress("GenParticle.motherIndex_2", GenParticle_motherIndex_2, &b_GenParticle_motherIndex_2);
   fChain->SetBranchAddress("nGenParticle", &nGenParticle, &b_fnGenParticle);
   fChain->SetBranchAddress("Jet", &Jet_, &b_Jet_);
   fChain->SetBranchAddress("Jet.fUniqueID", Jet_fUniqueID, &b_Jet_fUniqueID);
   fChain->SetBranchAddress("Jet.fBits", Jet_fBits, &b_Jet_fBits);
   fChain->SetBranchAddress("Jet.eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet.phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet.pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet.pt_raw", Jet_pt_raw, &b_Jet_pt_raw);
   fChain->SetBranchAddress("Jet.energy", Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet.energy_raw", Jet_energy_raw, &b_Jet_energy_raw);
   fChain->SetBranchAddress("Jet.jecUnc", Jet_jecUnc, &b_Jet_jecUnc);
   fChain->SetBranchAddress("Jet.resJEC", Jet_resJEC, &b_Jet_resJEC);
   fChain->SetBranchAddress("Jet.partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Jet.puIdMVA", Jet_puIdMVA, &b_Jet_puIdMVA);
   fChain->SetBranchAddress("Jet.puIdFlag", Jet_puIdFlag, &b_Jet_puIdFlag);
   fChain->SetBranchAddress("Jet.puIdBits", Jet_puIdBits, &b_Jet_puIdBits);
   fChain->SetBranchAddress("Jet.chargedEmEnergyFraction", Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet.chargedHadronEnergyFraction", Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet.chargedMuEnergyFraction", Jet_chargedMuEnergyFraction, &b_Jet_chargedMuEnergyFraction);
   fChain->SetBranchAddress("Jet.electronEnergyFraction", Jet_electronEnergyFraction, &b_Jet_electronEnergyFraction);
   fChain->SetBranchAddress("Jet.muonEnergyFraction", Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
   fChain->SetBranchAddress("Jet.neutralEmEnergyFraction", Jet_neutralEmEnergyFraction, &b_Jet_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jet.neutralHadronEnergyFraction", Jet_neutralHadronEnergyFraction, &b_Jet_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("Jet.photonEnergyFraction", Jet_photonEnergyFraction, &b_Jet_photonEnergyFraction);
   fChain->SetBranchAddress("Jet.chargedHadronMultiplicity", Jet_chargedHadronMultiplicity, &b_Jet_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jet.chargedMultiplicity", Jet_chargedMultiplicity, &b_Jet_chargedMultiplicity);
   fChain->SetBranchAddress("Jet.electronMultiplicity", Jet_electronMultiplicity, &b_Jet_electronMultiplicity);
   fChain->SetBranchAddress("Jet.muonMultiplicity", Jet_muonMultiplicity, &b_Jet_muonMultiplicity);
   fChain->SetBranchAddress("Jet.neutralHadronMultiplicity", Jet_neutralHadronMultiplicity, &b_Jet_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jet.neutralMultiplicity", Jet_neutralMultiplicity, &b_Jet_neutralMultiplicity);
   fChain->SetBranchAddress("Jet.photonMultiplicity", Jet_photonMultiplicity, &b_Jet_photonMultiplicity);
   fChain->SetBranchAddress("Jet.nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet.trackCountingHighEffBTag", Jet_trackCountingHighEffBTag, &b_Jet_trackCountingHighEffBTag);
   fChain->SetBranchAddress("Jet.trackCountingHighPurBTag", Jet_trackCountingHighPurBTag, &b_Jet_trackCountingHighPurBTag);
   fChain->SetBranchAddress("Jet.simpleSecondaryVertexHighEffBTag", Jet_simpleSecondaryVertexHighEffBTag, &b_Jet_simpleSecondaryVertexHighEffBTag);
   fChain->SetBranchAddress("Jet.simpleSecondaryVertexHighPurBTag", Jet_simpleSecondaryVertexHighPurBTag, &b_Jet_simpleSecondaryVertexHighPurBTag);
   fChain->SetBranchAddress("Jet.jetProbabilityBTag", Jet_jetProbabilityBTag, &b_Jet_jetProbabilityBTag);
   fChain->SetBranchAddress("Jet.jetBProbabilityBTag", Jet_jetBProbabilityBTag, &b_Jet_jetBProbabilityBTag);
   fChain->SetBranchAddress("Jet.combinedSecondaryVertexBTag", Jet_combinedSecondaryVertexBTag, &b_Jet_combinedSecondaryVertexBTag);
   fChain->SetBranchAddress("Jet.combinedSecondaryVertexMVABTag", Jet_combinedSecondaryVertexMVABTag, &b_Jet_combinedSecondaryVertexMVABTag);
   fChain->SetBranchAddress("Jet.combinedInclusiveSecondaryVertexBTag", Jet_combinedInclusiveSecondaryVertexBTag, &b_Jet_combinedInclusiveSecondaryVertexBTag);
   fChain->SetBranchAddress("Jet.combinedMVABTag", Jet_combinedMVABTag, &b_Jet_combinedMVABTag);
   fChain->SetBranchAddress("Jet.passLooseID", Jet_passLooseID, &b_Jet_passLooseID);
   fChain->SetBranchAddress("Jet.passTightID", Jet_passTightID, &b_Jet_passTightID);
   fChain->SetBranchAddress("Jet.selbit", Jet_selbit, &b_Jet_selbit);
   fChain->SetBranchAddress("nJet", &nJet, &b_fnJet);
   fChain->SetBranchAddress("MET", &MET_, &b_MET_);
   fChain->SetBranchAddress("MET.fUniqueID", MET_fUniqueID, &b_MET_fUniqueID);
   fChain->SetBranchAddress("MET.fBits", MET_fBits, &b_MET_fBits);
   fChain->SetBranchAddress("MET.met", MET_met, &b_MET_met);
   fChain->SetBranchAddress("MET.metphi", MET_metphi, &b_MET_metphi);
   fChain->SetBranchAddress("MET.sumet", MET_sumet, &b_MET_sumet);
   fChain->SetBranchAddress("MET.metuncorr", MET_metuncorr, &b_MET_metuncorr);
   fChain->SetBranchAddress("MET.metphiuncorr", MET_metphiuncorr, &b_MET_metphiuncorr);
   fChain->SetBranchAddress("MET.sumetuncorr", MET_sumetuncorr, &b_MET_sumetuncorr);
   fChain->SetBranchAddress("nMET", &nMET, &b_fnPFMET);
   fChain->SetBranchAddress("corrMET", &corrMET_, &b_corrMET_);
   fChain->SetBranchAddress("corrMET.fUniqueID", corrMET_fUniqueID, &b_corrMET_fUniqueID);
   fChain->SetBranchAddress("corrMET.fBits", corrMET_fBits, &b_corrMET_fBits);
   fChain->SetBranchAddress("corrMET.met", corrMET_met, &b_corrMET_met);
   fChain->SetBranchAddress("corrMET.metphi", corrMET_metphi, &b_corrMET_metphi);
   fChain->SetBranchAddress("corrMET.sumet", corrMET_sumet, &b_corrMET_sumet);
   fChain->SetBranchAddress("corrMET.metuncorr", corrMET_metuncorr, &b_corrMET_metuncorr);
   fChain->SetBranchAddress("corrMET.metphiuncorr", corrMET_metphiuncorr, &b_corrMET_metphiuncorr);
   fChain->SetBranchAddress("corrMET.sumetuncorr", corrMET_sumetuncorr, &b_corrMET_sumetuncorr);
   fChain->SetBranchAddress("corrnMET", &corrnMET, &b_fnCorrMET);
   fChain->SetBranchAddress("mvaMET", &mvaMET_, &b_mvaMET_);
   fChain->SetBranchAddress("mvaMET.fUniqueID", mvaMET_fUniqueID, &b_mvaMET_fUniqueID);
   fChain->SetBranchAddress("mvaMET.fBits", mvaMET_fBits, &b_mvaMET_fBits);
   fChain->SetBranchAddress("mvaMET.met", mvaMET_met, &b_mvaMET_met);
   fChain->SetBranchAddress("mvaMET.metphi", mvaMET_metphi, &b_mvaMET_metphi);
   fChain->SetBranchAddress("mvaMET.sumet", mvaMET_sumet, &b_mvaMET_sumet);
   fChain->SetBranchAddress("mvaMET.metuncorr", mvaMET_metuncorr, &b_mvaMET_metuncorr);
   fChain->SetBranchAddress("mvaMET.metphiuncorr", mvaMET_metphiuncorr, &b_mvaMET_metphiuncorr);
   fChain->SetBranchAddress("mvaMET.sumetuncorr", mvaMET_sumetuncorr, &b_mvaMET_sumetuncorr);
   fChain->SetBranchAddress("mvanMET", &mvanMET, &b_fnMVAMET);
   fChain->SetBranchAddress("Muon", &Muon_, &b_Muon_);
   fChain->SetBranchAddress("Muon.fUniqueID", Muon_fUniqueID, &b_Muon_fUniqueID);
   fChain->SetBranchAddress("Muon.fBits", Muon_fBits, &b_Muon_fBits);
   fChain->SetBranchAddress("Muon.isTrackerMuon", Muon_isTrackerMuon, &b_Muon_isTrackerMuon);
   fChain->SetBranchAddress("Muon.isPFMuon", Muon_isPFMuon, &b_Muon_isPFMuon);
   fChain->SetBranchAddress("Muon.eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon.phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon.pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon.ptError", Muon_ptError, &b_Muon_ptError);
   fChain->SetBranchAddress("Muon.p", Muon_p, &b_Muon_p);
   fChain->SetBranchAddress("Muon.energy", Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon.charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon.trkD0", Muon_trkD0, &b_Muon_trkD0);
   fChain->SetBranchAddress("Muon.trkD0Error", Muon_trkD0Error, &b_Muon_trkD0Error);
   fChain->SetBranchAddress("Muon.trkDz", Muon_trkDz, &b_Muon_trkDz);
   fChain->SetBranchAddress("Muon.trkDzError", Muon_trkDzError, &b_Muon_trkDzError);
   fChain->SetBranchAddress("Muon.globalChi2", Muon_globalChi2, &b_Muon_globalChi2);
   fChain->SetBranchAddress("Muon.trkIso", Muon_trkIso, &b_Muon_trkIso);
   fChain->SetBranchAddress("Muon.ecalIso", Muon_ecalIso, &b_Muon_ecalIso);
   fChain->SetBranchAddress("Muon.hcalIso", Muon_hcalIso, &b_Muon_hcalIso);
   fChain->SetBranchAddress("Muon.hoIso", Muon_hoIso, &b_Muon_hoIso);
   fChain->SetBranchAddress("Muon.relIso", Muon_relIso, &b_Muon_relIso);
   fChain->SetBranchAddress("Muon.passID", Muon_passID, &b_Muon_passID);
   fChain->SetBranchAddress("Muon.vtxDist3D", Muon_vtxDist3D, &b_Muon_vtxDist3D);
   fChain->SetBranchAddress("Muon.vtxIndex", Muon_vtxIndex, &b_Muon_vtxIndex);
   fChain->SetBranchAddress("Muon.vtxDistZ", Muon_vtxDistZ, &b_Muon_vtxDistZ);
   fChain->SetBranchAddress("Muon.pixHits", Muon_pixHits, &b_Muon_pixHits);
   fChain->SetBranchAddress("Muon.trkHits", Muon_trkHits, &b_Muon_trkHits);
   fChain->SetBranchAddress("Muon.muoHits", Muon_muoHits, &b_Muon_muoHits);
   fChain->SetBranchAddress("Muon.matches", Muon_matches, &b_Muon_matches);
   fChain->SetBranchAddress("Muon.trackerLayersWithMeasurement", Muon_trackerLayersWithMeasurement, &b_Muon_trackerLayersWithMeasurement);
   fChain->SetBranchAddress("Muon.pfRelIso", Muon_pfRelIso, &b_Muon_pfRelIso);
   fChain->SetBranchAddress("Muon.vx", Muon_vx, &b_Muon_vx);
   fChain->SetBranchAddress("Muon.vy", Muon_vy, &b_Muon_vy);
   fChain->SetBranchAddress("Muon.vz", Muon_vz, &b_Muon_vz);
   fChain->SetBranchAddress("Muon.dB", Muon_dB, &b_Muon_dB);
   fChain->SetBranchAddress("Muon.edB", Muon_edB, &b_Muon_edB);
   fChain->SetBranchAddress("Muon.dB3d", Muon_dB3d, &b_Muon_dB3d);
   fChain->SetBranchAddress("Muon.edB3d", Muon_edB3d, &b_Muon_edB3d);
   fChain->SetBranchAddress("Muon.isGlobalMuonPromptTight", Muon_isGlobalMuonPromptTight, &b_Muon_isGlobalMuonPromptTight);
   fChain->SetBranchAddress("Muon.isAllArbitrated", Muon_isAllArbitrated, &b_Muon_isAllArbitrated);
   fChain->SetBranchAddress("Muon.nChambers", Muon_nChambers, &b_Muon_nChambers);
   fChain->SetBranchAddress("Muon.nMatches", Muon_nMatches, &b_Muon_nMatches);
   fChain->SetBranchAddress("Muon.nMatchedStations", Muon_nMatchedStations, &b_Muon_nMatchedStations);
   fChain->SetBranchAddress("Muon.stationMask", Muon_stationMask, &b_Muon_stationMask);
   fChain->SetBranchAddress("Muon.stationGapMaskDistance", Muon_stationGapMaskDistance, &b_Muon_stationGapMaskDistance);
   fChain->SetBranchAddress("Muon.stationGapMaskPull", Muon_stationGapMaskPull, &b_Muon_stationGapMaskPull);
   fChain->SetBranchAddress("Muon.muonID", Muon_muonID, &b_Muon_muonID);
   fChain->SetBranchAddress("Muon.idMVA", Muon_idMVA, &b_Muon_idMVA);
   fChain->SetBranchAddress("Muon.isoRingsMVA", Muon_isoRingsMVA, &b_Muon_isoRingsMVA);
   fChain->SetBranchAddress("Muon.isoRingsRadMVA", Muon_isoRingsRadMVA, &b_Muon_isoRingsRadMVA);
   fChain->SetBranchAddress("Muon.idIsoCombMVA", Muon_idIsoCombMVA, &b_Muon_idIsoCombMVA);
   fChain->SetBranchAddress("Muon.pfRelIso03v1", Muon_pfRelIso03v1, &b_Muon_pfRelIso03v1);
   fChain->SetBranchAddress("Muon.pfRelIso03v2", Muon_pfRelIso03v2, &b_Muon_pfRelIso03v2);
   fChain->SetBranchAddress("Muon.pfRelIsoDB03v1", Muon_pfRelIsoDB03v1, &b_Muon_pfRelIsoDB03v1);
   fChain->SetBranchAddress("Muon.pfRelIsoDB03v2", Muon_pfRelIsoDB03v2, &b_Muon_pfRelIsoDB03v2);
   fChain->SetBranchAddress("Muon.pfRelIso04v1", Muon_pfRelIso04v1, &b_Muon_pfRelIso04v1);
   fChain->SetBranchAddress("Muon.pfRelIso04v2", Muon_pfRelIso04v2, &b_Muon_pfRelIso04v2);
   fChain->SetBranchAddress("Muon.pfRelIsoDB04v1", Muon_pfRelIsoDB04v1, &b_Muon_pfRelIsoDB04v1);
   fChain->SetBranchAddress("Muon.pfRelIsoDB04v2", Muon_pfRelIsoDB04v2, &b_Muon_pfRelIsoDB04v2);
   fChain->SetBranchAddress("Muon.selbit", Muon_selbit, &b_Muon_selbit);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_fnMuon);
   fChain->SetBranchAddress("Tau", &Tau_, &b_Tau_);
   fChain->SetBranchAddress("Tau.fUniqueID", Tau_fUniqueID, &b_Tau_fUniqueID);
   fChain->SetBranchAddress("Tau.fBits", Tau_fBits, &b_Tau_fBits);
   fChain->SetBranchAddress("Tau.eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau.phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau.pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau.energy", Tau_energy, &b_Tau_energy);
   fChain->SetBranchAddress("Tau.charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau.mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau.leadTrkPt", Tau_leadTrkPt, &b_Tau_leadTrkPt);
   fChain->SetBranchAddress("Tau.leadTrkPtError", Tau_leadTrkPtError, &b_Tau_leadTrkPtError);
   fChain->SetBranchAddress("Tau.leadTrkEta", Tau_leadTrkEta, &b_Tau_leadTrkEta);
   fChain->SetBranchAddress("Tau.leadTrkPhi", Tau_leadTrkPhi, &b_Tau_leadTrkPhi);
   fChain->SetBranchAddress("Tau.leadTrkCharge", Tau_leadTrkCharge, &b_Tau_leadTrkCharge);
   fChain->SetBranchAddress("Tau.leadTrkD0", Tau_leadTrkD0, &b_Tau_leadTrkD0);
   fChain->SetBranchAddress("Tau.leadTrkD0Error", Tau_leadTrkD0Error, &b_Tau_leadTrkD0Error);
   fChain->SetBranchAddress("Tau.leadTrkDz", Tau_leadTrkDz, &b_Tau_leadTrkDz);
   fChain->SetBranchAddress("Tau.leadTrkDzError", Tau_leadTrkDzError, &b_Tau_leadTrkDzError);
   fChain->SetBranchAddress("Tau.vtxIndex", Tau_vtxIndex, &b_Tau_vtxIndex);
   fChain->SetBranchAddress("Tau.vtxDxy", Tau_vtxDxy, &b_Tau_vtxDxy);
   fChain->SetBranchAddress("Tau.vtxDz", Tau_vtxDz, &b_Tau_vtxDz);
   fChain->SetBranchAddress("Tau.leadChargedParticlePt", Tau_leadChargedParticlePt, &b_Tau_leadChargedParticlePt);
   fChain->SetBranchAddress("Tau.leadNeutralParticlePt", Tau_leadNeutralParticlePt, &b_Tau_leadNeutralParticlePt);
   fChain->SetBranchAddress("Tau.leadParticlePt", Tau_leadParticlePt, &b_Tau_leadParticlePt);
   fChain->SetBranchAddress("Tau.numChargedHadronsSignalCone", Tau_numChargedHadronsSignalCone, &b_Tau_numChargedHadronsSignalCone);
   fChain->SetBranchAddress("Tau.numNeutralHadronsSignalCone", Tau_numNeutralHadronsSignalCone, &b_Tau_numNeutralHadronsSignalCone);
   fChain->SetBranchAddress("Tau.numPhotonsSignalCone", Tau_numPhotonsSignalCone, &b_Tau_numPhotonsSignalCone);
   fChain->SetBranchAddress("Tau.numParticlesSignalCone", Tau_numParticlesSignalCone, &b_Tau_numParticlesSignalCone);
   fChain->SetBranchAddress("Tau.numChargedHadronsIsoCone", Tau_numChargedHadronsIsoCone, &b_Tau_numChargedHadronsIsoCone);
   fChain->SetBranchAddress("Tau.numNeutralHadronsIsoCone", Tau_numNeutralHadronsIsoCone, &b_Tau_numNeutralHadronsIsoCone);
   fChain->SetBranchAddress("Tau.numPhotonsIsoCone", Tau_numPhotonsIsoCone, &b_Tau_numPhotonsIsoCone);
   fChain->SetBranchAddress("Tau.numParticlesIsoCone", Tau_numParticlesIsoCone, &b_Tau_numParticlesIsoCone);
   fChain->SetBranchAddress("Tau.ptSumPFChargedHadronsIsoCone", Tau_ptSumPFChargedHadronsIsoCone, &b_Tau_ptSumPFChargedHadronsIsoCone);
   fChain->SetBranchAddress("Tau.ptSumPFNeutralHadronsIsoCone", Tau_ptSumPFNeutralHadronsIsoCone, &b_Tau_ptSumPFNeutralHadronsIsoCone);
   fChain->SetBranchAddress("Tau.ptSumPhotonsIsoCone", Tau_ptSumPhotonsIsoCone, &b_Tau_ptSumPhotonsIsoCone);
   fChain->SetBranchAddress("Tau.sigChHadCandPt[40]", Tau_sigChHadCandPt, &b_Tau_sigChHadCandPt);
   fChain->SetBranchAddress("Tau.sigChHadCandEta[40]", Tau_sigChHadCandEta, &b_Tau_sigChHadCandEta);
   fChain->SetBranchAddress("Tau.sigChHadCandPhi[40]", Tau_sigChHadCandPhi, &b_Tau_sigChHadCandPhi);
   fChain->SetBranchAddress("Tau.sigNeHadCandPt[20]", Tau_sigNeHadCandPt, &b_Tau_sigNeHadCandPt);
   fChain->SetBranchAddress("Tau.sigNeHadCandEta[20]", Tau_sigNeHadCandEta, &b_Tau_sigNeHadCandEta);
   fChain->SetBranchAddress("Tau.sigNeHadCandPhi[20]", Tau_sigNeHadCandPhi, &b_Tau_sigNeHadCandPhi);
   fChain->SetBranchAddress("Tau.sigGammaCandPt[20]", Tau_sigGammaCandPt, &b_Tau_sigGammaCandPt);
   fChain->SetBranchAddress("Tau.sigGammaCandEta[20]", Tau_sigGammaCandEta, &b_Tau_sigGammaCandEta);
   fChain->SetBranchAddress("Tau.sigGammaCandPhi[20]", Tau_sigGammaCandPhi, &b_Tau_sigGammaCandPhi);
   fChain->SetBranchAddress("Tau.isoChHadCandPt[40]", Tau_isoChHadCandPt, &b_Tau_isoChHadCandPt);
   fChain->SetBranchAddress("Tau.isoChHadCandEta[40]", Tau_isoChHadCandEta, &b_Tau_isoChHadCandEta);
   fChain->SetBranchAddress("Tau.isoChHadCandPhi[40]", Tau_isoChHadCandPhi, &b_Tau_isoChHadCandPhi);
   fChain->SetBranchAddress("Tau.isoNeHadCandPt[20]", Tau_isoNeHadCandPt, &b_Tau_isoNeHadCandPt);
   fChain->SetBranchAddress("Tau.isoNeHadCandEta[20]", Tau_isoNeHadCandEta, &b_Tau_isoNeHadCandEta);
   fChain->SetBranchAddress("Tau.isoNeHadCandPhi[20]", Tau_isoNeHadCandPhi, &b_Tau_isoNeHadCandPhi);
   fChain->SetBranchAddress("Tau.isoGammaCandPt[20]", Tau_isoGammaCandPt, &b_Tau_isoGammaCandPt);
   fChain->SetBranchAddress("Tau.isoGammaCandEta[20]", Tau_isoGammaCandEta, &b_Tau_isoGammaCandEta);
   fChain->SetBranchAddress("Tau.isoGammaCandPhi[20]", Tau_isoGammaCandPhi, &b_Tau_isoGammaCandPhi);
   fChain->SetBranchAddress("Tau.decayModeFinding", Tau_decayModeFinding, &b_Tau_decayModeFinding);
   fChain->SetBranchAddress("Tau.decayModeFindingNewDMs", Tau_decayModeFindingNewDMs, &b_Tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("Tau.decayModeFindingOldDMs", Tau_decayModeFindingOldDMs, &b_Tau_decayModeFindingOldDMs);
   fChain->SetBranchAddress("Tau.chargedIsoPtSum", Tau_chargedIsoPtSum, &b_Tau_chargedIsoPtSum);
   fChain->SetBranchAddress("Tau.neutralIsoPtSum", Tau_neutralIsoPtSum, &b_Tau_neutralIsoPtSum);
   fChain->SetBranchAddress("Tau.puCorrPtSum", Tau_puCorrPtSum, &b_Tau_puCorrPtSum);
   fChain->SetBranchAddress("Tau.ChHadCandPt3Prong_1track", Tau_ChHadCandPt3Prong_1track, &b_Tau_ChHadCandPt3Prong_1track);
   fChain->SetBranchAddress("Tau.ChHadCandPt3Prong_2track", Tau_ChHadCandPt3Prong_2track, &b_Tau_ChHadCandPt3Prong_2track);
   fChain->SetBranchAddress("Tau.ChHadCandPt3Prong_3track", Tau_ChHadCandPt3Prong_3track, &b_Tau_ChHadCandPt3Prong_3track);
   fChain->SetBranchAddress("Tau.ChHadCandEta3Prong_1track", Tau_ChHadCandEta3Prong_1track, &b_Tau_ChHadCandEta3Prong_1track);
   fChain->SetBranchAddress("Tau.ChHadCandEta3Prong_2track", Tau_ChHadCandEta3Prong_2track, &b_Tau_ChHadCandEta3Prong_2track);
   fChain->SetBranchAddress("Tau.ChHadCandEta3Prong_3track", Tau_ChHadCandEta3Prong_3track, &b_Tau_ChHadCandEta3Prong_3track);
   fChain->SetBranchAddress("Tau.ChHadCandPhi3Prong_1track", Tau_ChHadCandPhi3Prong_1track, &b_Tau_ChHadCandPhi3Prong_1track);
   fChain->SetBranchAddress("Tau.ChHadCandPhi3Prong_2track", Tau_ChHadCandPhi3Prong_2track, &b_Tau_ChHadCandPhi3Prong_2track);
   fChain->SetBranchAddress("Tau.ChHadCandPhi3Prong_3track", Tau_ChHadCandPhi3Prong_3track, &b_Tau_ChHadCandPhi3Prong_3track);
   fChain->SetBranchAddress("Tau.ChHadCandPt1Prong", Tau_ChHadCandPt1Prong, &b_Tau_ChHadCandPt1Prong);
   fChain->SetBranchAddress("Tau.ChHadCandEta1Prong", Tau_ChHadCandEta1Prong, &b_Tau_ChHadCandEta1Prong);
   fChain->SetBranchAddress("Tau.ChHadCandPhi1Prong", Tau_ChHadCandPhi1Prong, &b_Tau_ChHadCandPhi1Prong);
   fChain->SetBranchAddress("Tau.CombinedIsolationDeltaBetaCorrRaw3Hits", Tau_CombinedIsolationDeltaBetaCorrRaw3Hits, &b_Tau_CombinedIsolationDeltaBetaCorrRaw3Hits);
   fChain->SetBranchAddress("Tau.CombinedIsolationDeltaBetaCorrRaw", Tau_CombinedIsolationDeltaBetaCorrRaw, &b_Tau_CombinedIsolationDeltaBetaCorrRaw);
   fChain->SetBranchAddress("Tau.byIsolationMVA3newDMwLTraw", Tau_byIsolationMVA3newDMwLTraw, &b_Tau_byIsolationMVA3newDMwLTraw);
   fChain->SetBranchAddress("Tau.byIsolationMVA3newDMwoLTraw", Tau_byIsolationMVA3newDMwoLTraw, &b_Tau_byIsolationMVA3newDMwoLTraw);
   fChain->SetBranchAddress("Tau.byIsolationMVA3oldDMwLTraw", Tau_byIsolationMVA3oldDMwLTraw, &b_Tau_byIsolationMVA3oldDMwLTraw);
   fChain->SetBranchAddress("Tau.byIsolationMVA3oldDMwoLTraw", Tau_byIsolationMVA3oldDMwoLTraw, &b_Tau_byIsolationMVA3oldDMwoLTraw);
   fChain->SetBranchAddress("Tau.byLooseCombinedIsolationDeltaBetaCorr", Tau_byLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau.byLooseCombinedIsolationDeltaBetaCorr3Hits", Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau.byMediumCombinedIsolationDeltaBetaCorr", Tau_byMediumCombinedIsolationDeltaBetaCorr, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau.byMediumCombinedIsolationDeltaBetaCorr3Hits", Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau.byTightCombinedIsolationDeltaBetaCorr", Tau_byTightCombinedIsolationDeltaBetaCorr, &b_Tau_byTightCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau.byTightCombinedIsolationDeltaBetaCorr3Hits", Tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau.byLooseIsolationMVA3newDMwLT", Tau_byLooseIsolationMVA3newDMwLT, &b_Tau_byLooseIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("Tau.byLooseIsolationMVA3newDMwoLT", Tau_byLooseIsolationMVA3newDMwoLT, &b_Tau_byLooseIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau.byLooseIsolationMva3oldDMwLT", Tau_byLooseIsolationMva3oldDMwLT, &b_Tau_byLooseIsolationMva3oldDMwLT);
   fChain->SetBranchAddress("Tau.byLooseIsolationMVA3oldDMwoLT", Tau_byLooseIsolationMVA3oldDMwoLT, &b_Tau_byLooseIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau.byMediumIsolationMVA3newDMwLT", Tau_byMediumIsolationMVA3newDMwLT, &b_Tau_byMediumIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("Tau.byMediumIsolationMVA3newDMwoLT", Tau_byMediumIsolationMVA3newDMwoLT, &b_Tau_byMediumIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau.byMediumIsolationMva3oldDMwLT", Tau_byMediumIsolationMva3oldDMwLT, &b_Tau_byMediumIsolationMva3oldDMwLT);
   fChain->SetBranchAddress("Tau.byMediumIsolationMVA3oldDMwoLT", Tau_byMediumIsolationMVA3oldDMwoLT, &b_Tau_byMediumIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau.byTightIsolationMVA3newDMwLT", Tau_byTightIsolationMVA3newDMwLT, &b_Tau_byTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("Tau.byTightIsolationMVA3newDMwoLT", Tau_byTightIsolationMVA3newDMwoLT, &b_Tau_byTightIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau.byTightIsolationMva3oldDMwLT", Tau_byTightIsolationMva3oldDMwLT, &b_Tau_byTightIsolationMva3oldDMwLT);
   fChain->SetBranchAddress("Tau.byTightIsolationMVA3oldDMwoLT", Tau_byTightIsolationMVA3oldDMwoLT, &b_Tau_byTightIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau.againstMuonLoose", Tau_againstMuonLoose, &b_Tau_againstMuonLoose);
   fChain->SetBranchAddress("Tau.againstMuonLoose2", Tau_againstMuonLoose2, &b_Tau_againstMuonLoose2);
   fChain->SetBranchAddress("Tau.againstMuonLoose3", Tau_againstMuonLoose3, &b_Tau_againstMuonLoose3);
   fChain->SetBranchAddress("Tau.againstMuonLooseMVA", Tau_againstMuonLooseMVA, &b_Tau_againstMuonLooseMVA);
   fChain->SetBranchAddress("Tau.againstMuonTight", Tau_againstMuonTight, &b_Tau_againstMuonTight);
   fChain->SetBranchAddress("Tau.againstMuonTight2", Tau_againstMuonTight2, &b_Tau_againstMuonTight2);
   fChain->SetBranchAddress("Tau.againstMuonTight3", Tau_againstMuonTight3, &b_Tau_againstMuonTight3);
   fChain->SetBranchAddress("Tau.againstMuonTightMVA", Tau_againstMuonTightMVA, &b_Tau_againstMuonTightMVA);
   fChain->SetBranchAddress("Tau.againstElectronLoose", Tau_againstElectronLoose, &b_Tau_againstElectronLoose);
   fChain->SetBranchAddress("Tau.againstElectronMedium", Tau_againstElectronMedium, &b_Tau_againstElectronMedium);
   fChain->SetBranchAddress("Tau.againstElectronTight", Tau_againstElectronTight, &b_Tau_againstElectronTight);
   fChain->SetBranchAddress("Tau.pfElectronMVA", Tau_pfElectronMVA, &b_Tau_pfElectronMVA);
   fChain->SetBranchAddress("Tau.againstElectronMVALooseMVA5", Tau_againstElectronMVALooseMVA5, &b_Tau_againstElectronMVALooseMVA5);
   fChain->SetBranchAddress("Tau.againstElectronMVAMediumMVA5", Tau_againstElectronMVAMediumMVA5, &b_Tau_againstElectronMVAMediumMVA5);
   fChain->SetBranchAddress("Tau.againstElectronMVATightMVA5", Tau_againstElectronMVATightMVA5, &b_Tau_againstElectronMVATightMVA5);
   fChain->SetBranchAddress("Tau.byVLooseCombinedIsolationDeltaBetaCorr", Tau_byVLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byVLooseCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau.byLooseIsolationMVA", Tau_byLooseIsolationMVA, &b_Tau_byLooseIsolationMVA);
   fChain->SetBranchAddress("Tau.byMediumIsolationMVA", Tau_byMediumIsolationMVA, &b_Tau_byMediumIsolationMVA);
   fChain->SetBranchAddress("Tau.byTightIsolationMVA", Tau_byTightIsolationMVA, &b_Tau_byTightIsolationMVA);
   fChain->SetBranchAddress("Tau.jetPt", Tau_jetPt, &b_Tau_jetPt);
   fChain->SetBranchAddress("Tau.jetEta", Tau_jetEta, &b_Tau_jetEta);
   fChain->SetBranchAddress("Tau.jetPhi", Tau_jetPhi, &b_Tau_jetPhi);
   fChain->SetBranchAddress("Tau.emFraction", Tau_emFraction, &b_Tau_emFraction);
   fChain->SetBranchAddress("Tau.maximumHCALPFClusterEt", Tau_maximumHCALPFClusterEt, &b_Tau_maximumHCALPFClusterEt);
   fChain->SetBranchAddress("Tau.ecalStripSumEOverPLead", Tau_ecalStripSumEOverPLead, &b_Tau_ecalStripSumEOverPLead);
   fChain->SetBranchAddress("Tau.bremsRecoveryEOverPLead", Tau_bremsRecoveryEOverPLead, &b_Tau_bremsRecoveryEOverPLead);
   fChain->SetBranchAddress("Tau.hcalTotOverPLead", Tau_hcalTotOverPLead, &b_Tau_hcalTotOverPLead);
   fChain->SetBranchAddress("Tau.hcalMaxOverPLead", Tau_hcalMaxOverPLead, &b_Tau_hcalMaxOverPLead);
   fChain->SetBranchAddress("Tau.hcal3x3OverPLead", Tau_hcal3x3OverPLead, &b_Tau_hcal3x3OverPLead);
   fChain->SetBranchAddress("Tau.etaetaMoment", Tau_etaetaMoment, &b_Tau_etaetaMoment);
   fChain->SetBranchAddress("Tau.phiphiMoment", Tau_phiphiMoment, &b_Tau_phiphiMoment);
   fChain->SetBranchAddress("Tau.etaphiMoment", Tau_etaphiMoment, &b_Tau_etaphiMoment);
   fChain->SetBranchAddress("Tau.vx", Tau_vx, &b_Tau_vx);
   fChain->SetBranchAddress("Tau.vy", Tau_vy, &b_Tau_vy);
   fChain->SetBranchAddress("Tau.vz", Tau_vz, &b_Tau_vz);
   fChain->SetBranchAddress("Tau.zvertex", Tau_zvertex, &b_Tau_zvertex);
   fChain->SetBranchAddress("Tau.ltsipt", Tau_ltsipt, &b_Tau_ltsipt);
   fChain->SetBranchAddress("Tau.TauVertexPos_x", Tau_TauVertexPos_x, &b_Tau_TauVertexPos_x);
   fChain->SetBranchAddress("Tau.TauVertexPos_y", Tau_TauVertexPos_y, &b_Tau_TauVertexPos_y);
   fChain->SetBranchAddress("Tau.TauVertexPos_z", Tau_TauVertexPos_z, &b_Tau_TauVertexPos_z);
   fChain->SetBranchAddress("Tau.PVPos_x", Tau_PVPos_x, &b_Tau_PVPos_x);
   fChain->SetBranchAddress("Tau.PVPos_y", Tau_PVPos_y, &b_Tau_PVPos_y);
   fChain->SetBranchAddress("Tau.PVPos_z", Tau_PVPos_z, &b_Tau_PVPos_z);
   fChain->SetBranchAddress("Tau.FlightLengh", Tau_FlightLengh, &b_Tau_FlightLengh);
   fChain->SetBranchAddress("Tau.FlightLenghtSig", Tau_FlightLenghtSig, &b_Tau_FlightLenghtSig);
   fChain->SetBranchAddress("Tau.DXY", Tau_DXY, &b_Tau_DXY);
   fChain->SetBranchAddress("Tau.DXY_err", Tau_DXY_err, &b_Tau_DXY_err);
   fChain->SetBranchAddress("Tau.DXY_sig", Tau_DXY_sig, &b_Tau_DXY_sig);
   fChain->SetBranchAddress("Tau.NumChHad", Tau_NumChHad, &b_Tau_NumChHad);
   fChain->SetBranchAddress("Tau.selbit", Tau_selbit, &b_Tau_selbit);
   fChain->SetBranchAddress("nTau", &nTau, &b_fnTau);
   fChain->SetBranchAddress("l1physbits", &l1physbits, &b_l1physbits);
   fChain->SetBranchAddress("l1techbits", &l1techbits, &b_l1techbits);
   fChain->SetBranchAddress("hltpaths", &hltpaths, &b_hltpaths);
   fChain->SetBranchAddress("hltresults", &hltresults, &b_hltresults);
   fChain->SetBranchAddress("hltprescales", &hltprescales, &b_hltprescales);
   Notify();
}

Bool_t EventTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventTree_cxx
