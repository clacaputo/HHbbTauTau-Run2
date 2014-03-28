#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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

#define SMART_TREE_FOR_CMSSW
#include "HHbbTauTau/TreeProduction/interface/Tau.h"

#define SIMPLE_VAR(type, name, default_value) tauTree.name() = patTau.tauID(#name);

class TauBlock : public edm::EDAnalyzer {
public:
    explicit TauBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("patTauSrc")),
        _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexSrc")) {}

private:
    virtual void endJob() { tauTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;
    edm::InputTag _vtxInputTag;

    const TransientTrackBuilder* trackBuilder_;
    ntuple::TauTree tauTree;
};

void TauBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    tauTree.EventId() = iEvent.id().event();

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel(_vtxInputTag, primaryVertices);

    edm::Handle<pat::TauCollection> tausHandle;
    iEvent.getByLabel(_inputTag,tausHandle);
    const pat::TauCollection* taus = tausHandle.product();
    if(!taus)
        throw std::runtime_error("Tau collection not found.");

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    
    for (const pat::Tau& patTau : *taus) {
      // Store Tau variables
      tauTree.eta()    = patTau.eta();
      tauTree.phi()    = patTau.phi();
      tauTree.pt()     = patTau.pt();
      tauTree.energy() = patTau.energy();
      tauTree.charge() = patTau.charge();
      
      

      // Leading particle pT
      tauTree.leadChargedParticlePt() = patTau.leadPFChargedHadrCand().isNonnull()
                                  ? patTau.leadPFChargedHadrCand()->pt() : 0.;
      tauTree.leadNeutralParticleEt() = patTau.leadPFNeutralCand().isNonnull()
                                  ? patTau.leadPFNeutralCand()->et(): 0.;
      tauTree.leadParticleEt()        = patTau.leadPFCand().isNonnull()
                                  ? patTau.leadPFCand()->et(): 0.;
      
      // Number of charged/neutral candidates and photons in different cones
      tauTree.numChargedHadronsSignalCone() = patTau.signalPFChargedHadrCands().size();
      tauTree.numNeutralHadronsSignalCone() = patTau.signalPFNeutrHadrCands().size();
      tauTree.numPhotonsSignalCone()        = patTau.signalPFGammaCands().size();
      tauTree.numParticlesSignalCone()      = patTau.signalPFCands().size();
      
      tauTree.numChargedHadronsIsoCone() = patTau.isolationPFChargedHadrCands().size();
      tauTree.numNeutralHadronsIsoCone() = patTau.isolationPFNeutrHadrCands().size();
      tauTree.numPhotonsIsoCone()        = patTau.isolationPFGammaCands().size();
      tauTree.numParticlesIsoCone()      = patTau.isolationPFCands().size();
      
      tauTree.ptSumPFChargedHadronsIsoCone() = patTau.isolationPFChargedHadrCandsPtSum();
      tauTree.etSumPhotonsIsoCone()          = patTau.isolationPFGammaCandsEtSum();

      // tau id. discriminators
      // see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
      // for discriminator names see PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py
      TAU_DISCRIMINATOR_DATA()

      tauTree.leadPFCand_mva_e_pi() = patTau.leadPFCand().isNonnull() ? patTau.leadPFCand()->mva_e_pi() : 1.;

      // NEW quantities
      tauTree.emFraction()              = patTau.emFraction();
      tauTree.maximumHCALPFClusterEt()  = patTau.maximumHCALPFClusterEt();
      tauTree.ecalStripSumEOverPLead()  = patTau.ecalStripSumEOverPLead();
      tauTree.bremsRecoveryEOverPLead() = patTau.bremsRecoveryEOverPLead();
      tauTree.hcalTotOverPLead()        = patTau.hcalTotOverPLead();
      tauTree.hcalMaxOverPLead()        = patTau.hcalMaxOverPLead();
      tauTree.hcal3x3OverPLead()        = patTau.hcal3x3OverPLead();

      tauTree.etaetaMoment() = patTau.etaetaMoment();
      tauTree.phiphiMoment() = patTau.phiphiMoment();
      tauTree.etaphiMoment() = patTau.etaphiMoment();
      
      // Vertex information
      const reco::Candidate::Point& vertex = patTau.vertex();
      tauTree.vx() = vertex.x();
      tauTree.vy() = vertex.y();
      tauTree.vz() = vertex.z();

      edm::ESHandle<TransientTrackBuilder> trackBuilderHandle;
       iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilderHandle);
        
      for(size_t n = 0; n < patTau.signalPFChargedHadrCands().size(); ++n) {
          if(patTau.signalPFChargedHadrCands()[n]->trackRef().isNonnull()) {
            tauTree.ChHadCand_Pt().push_back(patTau.signalPFChargedHadrCands()[n]->trackRef()->pt());
            tauTree.ChHadCand_Eta().push_back(patTau.signalPFChargedHadrCands()[n]->trackRef()->eta());
            tauTree.ChHadCand_Phi().push_back(patTau.signalPFChargedHadrCands()[n]->trackRef()->phi());
          }
      }

      tauTree.Fill();
    }
}

#undef SIMPLE_VAR
#undef TAU_DISCRIMINATOR_DATA

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
