#include <iostream>

#include "TTree.h"
#include "TClonesArray.h"
#include "HHbbTauTau/TreeProduction/plugins/GenJetBlock.h"
#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

GenJetBlock::GenJetBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genJetSrc"))
{
}
void GenJetBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneGenJet = new TClonesArray("vhtm::GenJet");
  tree->Branch("GenJet", &cloneGenJet, 32000, 2);
  tree->Branch("nGenJet", &fnGenJet,  "fnGenJet/I");
}
void GenJetBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneGenJet->Clear();
  fnGenJet = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel(_inputTag, genJets);

    if (genJets.isValid()) {
      edm::LogInfo("GenJetBlock") << "Total # GenJets: " << genJets->size();

      for (reco::GenJetCollection::const_iterator it = genJets->begin(); 
                                                 it != genJets->end(); ++it) {

        if (fnGenJet == kMaxGenJet) {
	  edm::LogInfo("GenJetBlock") << "Too many Gen Jets, fnGenJet = " << fnGenJet; 
	  break;
        }
        genJetB = new ((*cloneGenJet)[fnGenJet++]) vhtm::GenJet();

        // fill in all the vectors
        genJetB->eta    = it->eta();
        genJetB->phi    = it->phi();
        genJetB->p      = it->p();
        genJetB->pt     = it->pt();
        double energy   = it->energy();
        genJetB->energy = energy;
        genJetB->emf    = (energy>0) ? it->emEnergy()/energy : 0;
        genJetB->hadf   = (energy>0) ? it->hadEnergy()/energy : 0;
      }
    } 
    else {
      edm::LogError("GenJetBlock") << "Error >> Failed to get GenJetCollection for label: " 
                                   << _inputTag;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenJetBlock);
