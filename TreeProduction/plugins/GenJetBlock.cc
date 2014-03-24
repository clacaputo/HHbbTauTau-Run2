#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#define SMART_TREE_FOR_CMSSW
#include "HHbbTauTau/TreeProduction/interface/GenJet.h"

class GenJetBlock : public edm::EDAnalyzer {
public:
    explicit GenJetBlock(const edm::ParameterSet& iConfig) :
        _verbosity(iConfig.getParameter<int>("verbosity")),
        _inputTag(iConfig.getParameter<edm::InputTag>("genJetSrc")) {}
private:
    virtual void endJob() { genEventTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;

    ntuple::GenJetTree genJetTree;
};

void GenJetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    genJetTree.EventId() = iEvent.id().event();

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel(_inputTag, genJets);

    if (!genJets.isValid()) {
        edm::LogError("GenJetBlock") << "Error >> Failed to get GenJetCollection for label: "
                                     << _inputTag;
        throw std::runtime_error("Failed to get GenJetCollection");
    }

    edm::LogInfo("GenJetBlock") << "Total # GenJets: " << genJets->size();

    for (const reco::GenJet& genJet : *genJets) {
        genJetTree.eta()    = genJet.eta();
        genJetTree.phi()    = genJet.phi();
        genJetTree.p()      = genJet.p();
        genJetTree.pt()     = genJet.pt();
        double energy   = genJet.energy();
        genJetTree.energy() = energy;
        genJetTree.emf()    = (energy>0) ? genJet.emEnergy()/energy : 0;
        genJetTree.hadf()   = (energy>0) ? genJet.hadEnergy()/energy : 0;

        genJetTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenJetBlock);
