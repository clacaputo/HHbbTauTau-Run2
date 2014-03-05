#include <iostream>
#include "TTree.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TauPFSpecific.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Utilities/General/interface/FileInPath.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


#include "HHbbTauTau/TreeProduction/plugins/TauBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

namespace tb 
{
  template<typename T>
  bool isValidRef(const edm::Ref<T>& ref)
  //  bool isValidRef(const PFCandidatePtr ref))
  {
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
  }
}


TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc"))
{
}
TauBlock::~TauBlock() { }
void TauBlock::beginJob() 
{

  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneTau = new TClonesArray("vhtm::Tau");
  tree->Branch("Tau", &cloneTau, 32000, 2);
  tree->Branch("nTau", &fnTau, "fnTau/I");
}
void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneTau->Clear();
  fnTau = 0;

  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  edm::Handle<pat::TauCollection> tausHandle;
  iEvent.getByLabel(_inputTag,tausHandle);
  const pat::TauCollection* taus = tausHandle.product();


  //  if (taus.isValid()) {
    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    
    for (std::vector<pat::Tau>::const_iterator it  = taus->begin(); 
                                               it != taus->end(); ++it) {
      tauB = new ((*cloneTau)[fnTau++]) vhtm::Tau();

      // Store Tau variables
      tauB->eta    = it->eta();
      tauB->phi    = it->phi();
      tauB->pt     = it->pt();
      tauB->energy = it->energy();
      tauB->charge = it->charge();
      
      

      // Leading particle pT
      tauB->leadChargedParticlePt = it->leadPFChargedHadrCand().isNonnull() 
	? it->leadPFChargedHadrCand()->pt(): 0.;
      tauB->leadNeutralParticlePt = it->leadPFNeutralCand().isNonnull() 
                                      ? it->leadPFNeutralCand()->et(): 0.;
      tauB->leadParticlePt        = it->leadPFCand().isNonnull() 
	? it->leadPFCand()->et(): 0.;     
      
      // Number of charged/neutral candidates and photons in different cones
      tauB->numChargedHadronsSignalCone = it->signalPFChargedHadrCands().size();
      tauB->numNeutralHadronsSignalCone = it->signalPFNeutrHadrCands().size();
      tauB->numPhotonsSignalCone        = it->signalPFGammaCands().size();
      tauB->numParticlesSignalCone      = it->signalPFCands().size();
      
      tauB->numChargedHadronsIsoCone = it->isolationPFChargedHadrCands().size();
      tauB->numNeutralHadronsIsoCone = it->isolationPFNeutrHadrCands().size();
      tauB->numPhotonsIsoCone        = it->isolationPFGammaCands().size();
      tauB->numParticlesIsoCone      = it->isolationPFCands().size();
      
      tauB->ptSumPFChargedHadronsIsoCone = it->isolationPFChargedHadrCandsPtSum();
      tauB->ptSumPhotonsIsoCone          = it->isolationPFGammaCandsEtSum();
      
      // Signal Constituents
      // Charged hadrons
      //int indx = 0;
     
      // for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iCand  = it->signalPFChargedHadrCands().begin(); 
      //    iCand != it->signalPFChargedHadrCands().end(); ++iCand,indx++) {
      // for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFChargedHadrCands().begin(); 
// 	   iCand != it->signalPFChargedHadrCands().end(); ++iCand,indx++) {
//         const reco::PFCandidate& cand = (**iCand);
//         if (indx < vhtm::Tau::kMaxPFChargedCand) {
//           tauB->sigChHadCandPt[indx]  = cand.pt();
//           tauB->sigChHadCandEta[indx] = cand.eta();
//           tauB->sigChHadCandPhi[indx] = cand.phi();

//         }
//       }  
//       // Neutral hadrons
//       indx = 0; // reset
//       //for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iCand  = it->signalPFNeutrHadrCands().begin(); 
//       //                                              iCand != it->signalPFNeutrHadrCands().end(); ++iCand,indx++) {
// 	for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFNeutrHadrCands().begin(); 
// 	     iCand != it->signalPFNeutrHadrCands().end(); ++iCand,indx++) {
// 	  const reco::PFCandidate& cand = (**iCand);
// 	  if (indx < vhtm::Tau::kMaxPFNeutralCand) {
// 	    tauB->sigNeHadCandPt[indx]  = cand.pt();
// 	    tauB->sigNeHadCandEta[indx] = cand.eta();
// 	    tauB->sigNeHadCandPhi[indx] = cand.phi();
// 	  }
// 	}  
// 	// Photons
// 	indx = 0; // reset
// 	for (reco::PFCandidateRefVector::const_iterator iCand  = it->signalPFGammaCands().begin(); 
// 	     iCand != it->signalPFGammaCands().end(); ++iCand,indx++) 
// 	  //for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iCand  = it->signalPFGammaCands().begin(); 
// 	  //     iCand != it->signalPFGammaCands().end(); ++iCand,indx++){
	  
	  
// 	  if (indx < vhtm::Tau::kMaxPFNeutralCand) {
// 	    const reco::PFCandidate& cand = (**iCand);
// 	    tauB->sigGammaCandPt[indx]  = cand.pt();
// 	    tauB->sigGammaCandEta[indx] = cand.eta();
// 	    tauB->sigGammaCandPhi[indx] = cand.phi();
// 	  }
//     } 
    // Isolation Constituents
      // Charged hadrons
   //  int  indx = 0; // reset
//     for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFChargedHadrCands().begin(); 
// 	 iCand != it->isolationPFChargedHadrCands().end(); ++iCand,indx++) {
//       //      for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iCand = it->isolationPFChargedHadrCands().begin(); 
//       //	   iCand != it->isolationPFChargedHadrCands().end(); ++iCand,indx++) {
//       const reco::PFCandidate& cand = (**iCand);
//       if (indx < vhtm::Tau::kMaxPFChargedCand) {
// 	tauB->isoChHadCandPt[indx]  = cand.pt();
// 	tauB->isoChHadCandEta[indx] = cand.eta();
// 	tauB->isoChHadCandPhi[indx] = cand.phi();
//       }
//     }  
    
//     // Neutral hadrons
//     double ptSum = 0;
//     indx = 0; // reset
//     for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFNeutrHadrCands().begin(); 
// 	 iCand != it->isolationPFNeutrHadrCands().end(); ++iCand,indx++) {
//       //      for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iCand = it->isolationPFNeutrHadrCands().begin(); 
//       // iCand != it->isolationPFNeutrHadrCands().end(); ++iCand,indx++) {
//       const reco::PFCandidate& cand = (**iCand);
//       if (indx < vhtm::Tau::kMaxPFNeutralCand) {
// 	tauB->isoNeHadCandPt[indx]  = cand.pt();
//           tauB->isoNeHadCandEta[indx] = cand.eta();
//           tauB->isoNeHadCandPhi[indx] = cand.phi();
//       }
//       ptSum += std::abs(cand.pt());
//     }  
//     tauB->ptSumPFNeutralHadronsIsoCone = ptSum;
    
//     // Photons
//     indx = 0; // reset
//     for (reco::PFCandidateRefVector::const_iterator iCand  = it->isolationPFGammaCands().begin(); 
// 	 iCand != it->isolationPFGammaCands().end(); ++iCand,indx++) {
//       //      for (std::vector<edm::Ptr<reco::PFCandidate> >::const_iterator iCand = it->isolationPFGammaCands().begin(); 
//       //  iCand != it->isolationPFGammaCands().end(); ++iCand,indx++) {
//       const reco::PFCandidate& cand = (**iCand);
//         if (indx < vhtm::Tau::kMaxPFNeutralCand) {
//           tauB->isoGammaCandPt[indx]  = cand.pt();
//           tauB->isoGammaCandEta[indx] = cand.eta();
//           tauB->isoGammaCandPhi[indx] = cand.phi();
//         }
//       }  

      // tau id. discriminators
      tauB->decayModeFinding = it->tauID("decayModeFinding");
      // tauB->decayModeFindingNewDMs = it->tauID("decayModeFindingNewDMs");
      //      tauB->decayModeFindingOldDMs = it->tauID("decayModeFindingOldDMs");

      // tauB->chargedIsoPtSum  =it->tauID("chargedIsoPtSum");
      //      tauB->neutralIsoPtSum =it->tauID("neutralIsoPtSum");
      //      tauB->puCorrPtSum =it->tauID("puCorrPtSum");

      tauB->CombinedIsolationDeltaBetaCorrRaw3Hits   = it->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
      tauB->CombinedIsolationDeltaBetaCorrRaw   = it->tauID("byCombinedIsolationDeltaBetaCorrRaw");
      // tauB->byIsolationMVA3newDMwLTraw = it->tauID("byIsolationMVA3newDMwLTraw");
      //tauB->byIsolationMVA3newDMwoLTraw=  it->tauID("byIsolationMVA3newDMwoLTraw");
      //      tauB->byIsolationMVA3oldDMwLTraw=  it->tauID("byIsolationMVA3oldDMwLTraw");
      //      tauB-> byIsolationMVA3oldDMwoLTraw = it->tauID("byIsolationMVA3oldDMwoLTraw");

      tauB->byLooseCombinedIsolationDeltaBetaCorr=it->tauID("byLooseCombinedIsolationDeltaBetaCorr");
      tauB->byLooseCombinedIsolationDeltaBetaCorr3Hits=it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
 
      tauB->byMediumCombinedIsolationDeltaBetaCorr=it->tauID("byMediumCombinedIsolationDeltaBetaCorr");
      tauB->byMediumCombinedIsolationDeltaBetaCorr3Hits=it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");


      tauB->byTightCombinedIsolationDeltaBetaCorr=it->tauID("byTightCombinedIsolationDeltaBetaCorr");
      tauB->byTightCombinedIsolationDeltaBetaCorr3Hits=it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");

      //  tauB->byLooseIsolationMVA3newDMwLT=it->tauID("byLooseIsolationMVA3newDMwLT");
      //  tauB->byLooseIsolationMVA3newDMwoLT=it->tauID("byLooseIsolationMVA3newDMwoLT");
      // tauB->byLooseIsolationMva3oldDMwLT=it->tauID("byLooseIsolationMVA3oldDMwLT");
      // tauB->byLooseIsolationMVA3oldDMwoLT=it->tauID("byLooseIsolationMVA3oldDMwoLT");

      //tauB->byMediumIsolationMVA3newDMwLT=it->tauID("byMediumIsolationMVA3newDMwLT");
      //tauB->byMediumIsolationMVA3newDMwoLT=it->tauID("byMediumIsolationMVA3newDMwoLT");
      //tauB->byMediumIsolationMva3oldDMwLT=it->tauID("byMediumIsolationMVA3oldDMwLT");
      //tauB->byMediumIsolationMVA3oldDMwoLT=it->tauID("byMediumIsolationMVA3oldDMwoLT");

      //tauB->byTightIsolationMVA3newDMwLT=it->tauID("byTightIsolationMVA3newDMwLT");
      //tauB->byTightIsolationMVA3newDMwoLT=it->tauID("byTightIsolationMVA3newDMwoLT");
      //tauB->byTightIsolationMva3oldDMwLT=it->tauID("byTightIsolationMVA3oldDMwLT");
      //tauB->byTightIsolationMVA3oldDMwoLT=it->tauID("byTightIsolationMVA3oldDMwoLT");


      // discriminators against electrons/muons
      tauB->againstMuonLoose      = it->tauID("againstMuonLoose");
      tauB->againstMuonLoose2      = it->tauID("againstMuonLoose2");
      tauB->againstMuonLoose3      = it->tauID("againstMuonLoose3");
      //      tauB->againstMuonLooseMVA      = it->tauID("againstMuonLooseMVA");
      tauB->againstMuonTight      = it->tauID("againstMuonTight");
      tauB->againstMuonTight2      = it->tauID("againstMuonTight2");
      tauB->againstMuonTight3      = it->tauID("againstMuonTight3");
      //      tauB->againstMuonTightMVA      = it->tauID("againstMuonTightMVA");
      tauB->againstElectronLoose  = it->tauID("againstElectronLoose");
      tauB->againstElectronMedium = it->tauID("againstElectronMedium");
      tauB->againstElectronTight  = it->tauID("againstElectronTight");
      tauB->pfElectronMVA         = it->leadPFCand().isNonnull() 
                                ? it->leadPFCand()->mva_e_pi() : 1.;
      
      // ElectronIDMVA, electron faking tau
      //tauB->againstElectronMVALooseMVA5    = it->tauID("againstElectronLooseMVA5");
      //tauB->againstElectronMVAMediumMVA5    = it->tauID("againstElectronLooseMVA5");
      //tauB->againstElectronMVATightMVA5    = it->tauID("againstElectronLooseMVA5");

      tauB->byVLooseCombinedIsolationDeltaBetaCorr = it->tauID("byVLooseCombinedIsolationDeltaBetaCorr");
  

      // NEW quantities
      tauB->emFraction              = it->emFraction(); 
      tauB->maximumHCALPFClusterEt  = it->maximumHCALPFClusterEt();
      tauB->ecalStripSumEOverPLead  = it->ecalStripSumEOverPLead();
      tauB->bremsRecoveryEOverPLead = it->bremsRecoveryEOverPLead();
      tauB->hcalTotOverPLead        = it->hcalTotOverPLead();
      tauB->hcalMaxOverPLead        = it->hcalMaxOverPLead();
      tauB->hcal3x3OverPLead        = it->hcal3x3OverPLead();

      tauB->etaetaMoment = it->etaetaMoment();
      tauB->phiphiMoment = it->phiphiMoment();
      tauB->phiphiMoment = it->etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = it->vertex();
      tauB->vx = vertex.x();             
      tauB->vy = vertex.y();             
      tauB->vz = vertex.z();             

      tauB->zvertex = it->vz(); // distance from the primary vertex
      tauB->mass    = it->p4().M();
      //      tauB->ltsipt  = TMath::Abs(it->leadPFChargedHadrCandsignedSipt());




     //  float PV_x_std=0,PV_y_std=0,PV_z_std=0;
//       PV_x_std=   it->primaryVertexPos().x();
//       PV_y_std=   it->primaryVertexPos().y();
//       PV_z_std=   it->primaryVertexPos().z();

//       float SV_x_std=0,SV_y_std=0,SV_z_std=0;
//       SV_x_std=   it->secondaryVertexPos().x();
//       SV_y_std=   it->secondaryVertexPos().y();
//       SV_z_std=   it->secondaryVertexPos().z();
//       float distSV_PVRef_std=TMath::Sqrt((PV_x_std- SV_x_std)*(PV_x_std- SV_x_std) + (PV_y_std- SV_y_std)*(PV_y_std- SV_y_std) + (PV_z_std-SV_z_std)*(PV_z_std-SV_z_std));

//       tauB->TauVertexPos_x= it->secondaryVertexPos().x();
//       tauB->TauVertexPos_y= it->secondaryVertexPos().y();
//       tauB->TauVertexPos_z= it->secondaryVertexPos().z();

//       tauB->PVPos_x=it->primaryVertexPos().x();
//       tauB->PVPos_y=it->primaryVertexPos().y();
//       tauB->PVPos_z=it->primaryVertexPos().z();

//       tauB->FlightLengh=distSV_PVRef_std;
//       tauB->FlightLenghtSig=it->flightLengthSig();
//       tauB->DXY=it->dxy();
//       tauB->DXY_err=it->dxy_error();
//       tauB->DXY_sig=(it->dxy())/it->dxy_error();

      edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
       iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
      
      std::vector<reco::TransientTrack> tracks;
      
      reco::TransientTrack transtrack0;
      reco::TransientTrack transtrack1;
      reco::TransientTrack transtrack2;
   
            std::vector<edm::Ptr<reco::PFCandidate> > signalPFChargedHadrCandidates =it->signalPFChargedHadrCands();
      //edm::RefVector<std::vector<reco::PFCandidate> > signalPFChargedHadrCandidates =it->signalPFChargedHadrCands();
      if(signalPFChargedHadrCandidates.size()==3){//controlla                                                                                                                                                                           
	if(it->signalPFChargedHadrCands()[0]->trackRef().isNonnull()) transtrack0 = trackBuilderHandle->build(it->signalPFChargedHadrCands()[0]->trackRef());
	if(it->signalPFChargedHadrCands()[1]->trackRef().isNonnull()) transtrack1 = trackBuilderHandle->build(it->signalPFChargedHadrCands()[1]->trackRef());
	if(it->signalPFChargedHadrCands()[2]->trackRef().isNonnull()) transtrack2 = trackBuilderHandle->build(it->signalPFChargedHadrCands()[2]->trackRef());
      }

      //std::cout<<" Num Sig Ch Had  "<<signalPFChargedHadrCandidates.size()<<std::endl;

      tauB->NumChHad=it->signalPFChargedHadrCands().size();

      if(signalPFChargedHadrCandidates.size()==1){
	if(it->signalPFChargedHadrCands()[0]->trackRef().isNonnull()) {
	  tauB->ChHadCandPt1Prong=it->signalPFChargedHadrCands()[0]->trackRef()->pt();
	  tauB->ChHadCandEta1Prong=it->signalPFChargedHadrCands()[0]->trackRef()->eta();
	  tauB->ChHadCandPhi1Prong=it->signalPFChargedHadrCands()[0]->trackRef()->phi();
	  // std::cout<<" ********track eta 1p "<<it->signalPFChargedHadrCands()[0]->trackRef()->eta()<<std::endl;
	}
      }

      if( ( it->signalPFChargedHadrCands().size() )==3) {
	  if(it->signalPFChargedHadrCands()[0]->trackRef().isNonnull()) { 
	    tauB->ChHadCandPt3Prong_1track=it->signalPFChargedHadrCands()[0]->trackRef()->pt();
	    tauB->ChHadCandEta3Prong_1track=it->signalPFChargedHadrCands()[0]->trackRef()->eta();
	    // std::cout<<" ******track eta 0 "<<it->signalPFChargedHadrCands()[0]->trackRef()->eta()<<std::endl;
	    tauB->ChHadCandPhi3Prong_1track=it->signalPFChargedHadrCands()[0]->trackRef()->phi();
	  }
	  if(it->signalPFChargedHadrCands()[1]->trackRef().isNonnull()) {
	    tauB->ChHadCandPt3Prong_2track=it->signalPFChargedHadrCands()[1]->trackRef()->pt();
	    tauB->ChHadCandEta3Prong_2track=it->signalPFChargedHadrCands()[1]->trackRef()->eta();
	    tauB->ChHadCandPhi3Prong_2track=it->signalPFChargedHadrCands()[1]->trackRef()->phi();
	    // std::cout<<" ********track eta 1 "<<it->signalPFChargedHadrCands()[1]->trackRef()->eta()<<std::endl;
	    
	  }
	  if(it->signalPFChargedHadrCands()[2]->trackRef().isNonnull())	{ 
	    tauB->ChHadCandPt3Prong_3track=it->signalPFChargedHadrCands()[2]->trackRef()->pt();
	    tauB->ChHadCandEta3Prong_3track=it->signalPFChargedHadrCands()[2]->trackRef()->eta();
	    tauB->ChHadCandPhi3Prong_3track=it->signalPFChargedHadrCands()[2]->trackRef()->phi();
	    // std::cout<<"*********** track eta 2 "<<it->signalPFChargedHadrCands()[2]->trackRef()->eta()<<std::endl;
	  }
	}

  

    }
    //}
    //  else {
    //    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: " 
    //                        << _inputTag;
    //}
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
