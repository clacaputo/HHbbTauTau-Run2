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

#define TAU_DISCRIMINATOR(name) \
    tauB->name = patTau.tauID(#name) \
    /**/

namespace {
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
  if(!taus)
      throw std::runtime_error("Tau collection not found.");

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    
    for (const pat::Tau& patTau : *taus) {
      vhtm::Tau* tauB = new ((*cloneTau)[fnTau++]) vhtm::Tau();

      // Store Tau variables
      tauB->eta    = patTau.eta();
      tauB->phi    = patTau.phi();
      tauB->pt     = patTau.pt();
      tauB->energy = patTau.energy();
      tauB->charge = patTau.charge();
      
      

      // Leading particle pT
      tauB->leadChargedParticlePt = patTau.leadPFChargedHadrCand().isNonnull()
                                  ? patTau.leadPFChargedHadrCand()->pt() : 0.;
      tauB->leadNeutralParticlePt = patTau.leadPFNeutralCand().isNonnull()
                                  ? patTau.leadPFNeutralCand()->et(): 0.;
      tauB->leadParticlePt        = patTau.leadPFCand().isNonnull()
                                  ? patTau.leadPFCand()->et(): 0.;
      
      // Number of charged/neutral candidates and photons in different cones
      tauB->numChargedHadronsSignalCone = patTau.signalPFChargedHadrCands().size();
      tauB->numNeutralHadronsSignalCone = patTau.signalPFNeutrHadrCands().size();
      tauB->numPhotonsSignalCone        = patTau.signalPFGammaCands().size();
      tauB->numParticlesSignalCone      = patTau.signalPFCands().size();
      
      tauB->numChargedHadronsIsoCone = patTau.isolationPFChargedHadrCands().size();
      tauB->numNeutralHadronsIsoCone = patTau.isolationPFNeutrHadrCands().size();
      tauB->numPhotonsIsoCone        = patTau.isolationPFGammaCands().size();
      tauB->numParticlesIsoCone      = patTau.isolationPFCands().size();
      
      tauB->ptSumPFChargedHadronsIsoCone = patTau.isolationPFChargedHadrCandsPtSum();
      tauB->ptSumPhotonsIsoCone          = patTau.isolationPFGammaCandsEtSum();


      // tau id. discriminators
      // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
      // for discriminator names see PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py
      TAU_DISCRIMINATOR(decayModeFinding);

      TAU_DISCRIMINATOR(againstElectronLoose);
      TAU_DISCRIMINATOR(againstElectronMedium);
      TAU_DISCRIMINATOR(againstElectronTight);
      TAU_DISCRIMINATOR(againstElectronLooseMVA5);
      TAU_DISCRIMINATOR(againstElectronMediumMVA5);
      TAU_DISCRIMINATOR(againstElectronTightMVA5);
      TAU_DISCRIMINATOR(againstElectronVTightMVA5);

      TAU_DISCRIMINATOR(againstMuonLoose);
      TAU_DISCRIMINATOR(againstMuonMedium);
      TAU_DISCRIMINATOR(againstMuonTight);
      TAU_DISCRIMINATOR(againstMuonLoose3);
      TAU_DISCRIMINATOR(againstMuonTight3);
      TAU_DISCRIMINATOR(againstMuonLooseMVA);
      TAU_DISCRIMINATOR(againstMuonMediumMVA);
      TAU_DISCRIMINATOR(againstMuonTightMVA);
      TAU_DISCRIMINATOR(againstMuonMVAraw);

      TAU_DISCRIMINATOR(byVLooseCombinedIsolationDeltaBetaCorr);
      TAU_DISCRIMINATOR(byLooseCombinedIsolationDeltaBetaCorr);
      TAU_DISCRIMINATOR(byMediumCombinedIsolationDeltaBetaCorr);
      TAU_DISCRIMINATOR(byTightCombinedIsolationDeltaBetaCorr);
      TAU_DISCRIMINATOR(byLooseCombinedIsolationDeltaBetaCorr3Hits);
      TAU_DISCRIMINATOR(byMediumCombinedIsolationDeltaBetaCorr3Hits);
      TAU_DISCRIMINATOR(byTightCombinedIsolationDeltaBetaCorr3Hits);
      TAU_DISCRIMINATOR(byIsolationMVA3oldDMwoLTraw);
      TAU_DISCRIMINATOR(byVLooseIsolationMVA3oldDMwoLT);
      TAU_DISCRIMINATOR(byLooseIsolationMVA3oldDMwoLT);
      TAU_DISCRIMINATOR(byMediumIsolationMVA3oldDMwoLT);
      TAU_DISCRIMINATOR(byTightIsolationMVA3oldDMwoLT);
      TAU_DISCRIMINATOR(byVTightIsolationMVA3oldDMwoLT);
      TAU_DISCRIMINATOR(byVVTightIsolationMVA3oldDMwoLT);
      TAU_DISCRIMINATOR(byIsolationMVA3oldDMwLTraw);
      TAU_DISCRIMINATOR(byVLooseIsolationMVA3oldDMwLT);
      TAU_DISCRIMINATOR(byLooseIsolationMVA3oldDMwLT);
      TAU_DISCRIMINATOR(byMediumIsolationMVA3oldDMwLT);
      TAU_DISCRIMINATOR(byTightIsolationMVA3oldDMwLT);
      TAU_DISCRIMINATOR(byVTightIsolationMVA3oldDMwLT);
      TAU_DISCRIMINATOR(byVVTightIsolationMVA3oldDMwLT);
      TAU_DISCRIMINATOR(byIsolationMVA3newDMwoLTraw);
      TAU_DISCRIMINATOR(byVLooseIsolationMVA3newDMwoLT);
      TAU_DISCRIMINATOR(byLooseIsolationMVA3newDMwoLT);
      TAU_DISCRIMINATOR(byMediumIsolationMVA3newDMwoLT);
      TAU_DISCRIMINATOR(byTightIsolationMVA3newDMwoLT);
      TAU_DISCRIMINATOR(byVTightIsolationMVA3newDMwoLT);
      TAU_DISCRIMINATOR(byVVTightIsolationMVA3newDMwoLT);
      TAU_DISCRIMINATOR(byIsolationMVA3newDMwLTraw);
      TAU_DISCRIMINATOR(byVLooseIsolationMVA3newDMwLT);
      TAU_DISCRIMINATOR(byLooseIsolationMVA3newDMwLT);
      TAU_DISCRIMINATOR(byMediumIsolationMVA3newDMwLT);
      TAU_DISCRIMINATOR(byTightIsolationMVA3newDMwLT);
      TAU_DISCRIMINATOR(byVTightIsolationMVA3newDMwLT);
      TAU_DISCRIMINATOR(byVVTightIsolationMVA3newDMwLT);


      tauB->pfElectronMVA = patTau.leadPFCand().isNonnull() ? patTau.leadPFCand()->mva_e_pi() : 1.;

      // NEW quantities
      tauB->emFraction              = patTau.emFraction();
      tauB->maximumHCALPFClusterEt  = patTau.maximumHCALPFClusterEt();
      tauB->ecalStripSumEOverPLead  = patTau.ecalStripSumEOverPLead();
      tauB->bremsRecoveryEOverPLead = patTau.bremsRecoveryEOverPLead();
      tauB->hcalTotOverPLead        = patTau.hcalTotOverPLead();
      tauB->hcalMaxOverPLead        = patTau.hcalMaxOverPLead();
      tauB->hcal3x3OverPLead        = patTau.hcal3x3OverPLead();

      tauB->etaetaMoment = patTau.etaetaMoment();
      tauB->phiphiMoment = patTau.phiphiMoment();
      tauB->phiphiMoment = patTau.etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = patTau.vertex();
      tauB->vx = vertex.x();             
      tauB->vy = vertex.y();             
      tauB->vz = vertex.z();             

      tauB->zvertex = patTau.vz(); // distance from the primary vertex
      tauB->mass    = patTau.p4().M();

      edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
       iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
      
      std::vector<reco::TransientTrack> tracks;
      
      reco::TransientTrack transtrack0;
      reco::TransientTrack transtrack1;
      reco::TransientTrack transtrack2;
   
      std::vector<edm::Ptr<reco::PFCandidate> > signalPFChargedHadrCandidates =patTau.signalPFChargedHadrCands();
      if(signalPFChargedHadrCandidates.size()==3){//controlla                                                                                                                                                                           
        if(patTau.signalPFChargedHadrCands()[0]->trackRef().isNonnull())
            transtrack0 = trackBuilderHandle->build(patTau.signalPFChargedHadrCands()[0]->trackRef());
        if(patTau.signalPFChargedHadrCands()[1]->trackRef().isNonnull())
            transtrack1 = trackBuilderHandle->build(patTau.signalPFChargedHadrCands()[1]->trackRef());
        if(patTau.signalPFChargedHadrCands()[2]->trackRef().isNonnull())
            transtrack2 = trackBuilderHandle->build(patTau.signalPFChargedHadrCands()[2]->trackRef());
      }

      tauB->NumChHad=patTau.signalPFChargedHadrCands().size();

      if(signalPFChargedHadrCandidates.size()==1){
        if(patTau.signalPFChargedHadrCands()[0]->trackRef().isNonnull()) {
          tauB->ChHadCandPt1Prong=patTau.signalPFChargedHadrCands()[0]->trackRef()->pt();
          tauB->ChHadCandEta1Prong=patTau.signalPFChargedHadrCands()[0]->trackRef()->eta();
          tauB->ChHadCandPhi1Prong=patTau.signalPFChargedHadrCands()[0]->trackRef()->phi();
        }
      }

      if( ( patTau.signalPFChargedHadrCands().size() )==3) {
          if(patTau.signalPFChargedHadrCands()[0]->trackRef().isNonnull()) {
            tauB->ChHadCandPt3Prong_1track=patTau.signalPFChargedHadrCands()[0]->trackRef()->pt();
            tauB->ChHadCandEta3Prong_1track=patTau.signalPFChargedHadrCands()[0]->trackRef()->eta();
            tauB->ChHadCandPhi3Prong_1track=patTau.signalPFChargedHadrCands()[0]->trackRef()->phi();
          }
          if(patTau.signalPFChargedHadrCands()[1]->trackRef().isNonnull()) {
            tauB->ChHadCandPt3Prong_2track=patTau.signalPFChargedHadrCands()[1]->trackRef()->pt();
            tauB->ChHadCandEta3Prong_2track=patTau.signalPFChargedHadrCands()[1]->trackRef()->eta();
            tauB->ChHadCandPhi3Prong_2track=patTau.signalPFChargedHadrCands()[1]->trackRef()->phi();

          }
          if(patTau.signalPFChargedHadrCands()[2]->trackRef().isNonnull())	{
            tauB->ChHadCandPt3Prong_3track=patTau.signalPFChargedHadrCands()[2]->trackRef()->pt();
            tauB->ChHadCandEta3Prong_3track=patTau.signalPFChargedHadrCands()[2]->trackRef()->eta();
            tauB->ChHadCandPhi3Prong_3track=patTau.signalPFChargedHadrCands()[2]->trackRef()->phi();
          }
      }
    }
}

#undef TAU_DISCRIMINATOR

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
