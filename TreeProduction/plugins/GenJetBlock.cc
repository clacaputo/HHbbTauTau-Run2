/*!
 * \file GenJetBlock.cc
 * \author Original author: Subir Sarkar
 * \author Contributing author: Konstantin Androsov (INFN Pisa, Siena University)
 * \author Contributing author: Maria Teresa Grippo (INFN Pisa, Siena University)
 */


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
    virtual void endJob() { genJetTree.Write(); }
    virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
    int _verbosity;
    edm::InputTag _inputTag;

    ntuple::GenJetTree genJetTree;
};

void GenJetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    if (iEvent.isRealData()) return;

    genJetTree.RunId() = iEvent.id().run();
    genJetTree.LumiBlock() = iEvent.id().luminosityBlock();
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
        genJetTree.pt()     = genJet.pt();
        genJetTree.eta()    = genJet.eta();
        genJetTree.phi()    = genJet.phi();
        genJetTree.energy() = genJet.energy();
        genJetTree.emEnergy()    = genJet.emEnergy();
        genJetTree.hadEnergy()   = genJet.hadEnergy();

        genJetTree.Fill();
    }
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenJetBlock);
