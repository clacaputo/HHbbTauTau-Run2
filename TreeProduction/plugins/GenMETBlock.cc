#include "TTree.h"
#include "TClonesArray.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"

#include "HHbbTauTau/TreeProduction/plugins/GenMETBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

GenMETBlock::GenMETBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _inputTag(iConfig.getParameter<edm::InputTag>("genMETSrc"))
{
}
void GenMETBlock::beginJob() 
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneGenMET = new TClonesArray("vhtm::GenMET");
  tree->Branch("GenMET", &cloneGenMET, 32000, 2);
  tree->Branch("nGenMET", &fnGenMET,  "fnGenMET/I");
}
void GenMETBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  cloneGenMET->Clear();
  fnGenMET = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenMETCollection> mets;
    iEvent.getByLabel(_inputTag, mets);
    if (mets.isValid()) {
      edm::LogInfo("GenMETBlock") << "Total # GenMETs: " << mets->size();
      for (reco::GenMETCollection::const_iterator it = mets->begin(); 
                                                 it != mets->end(); ++it ) {
        if (fnGenMET == kMaxGenMET) {
	  edm::LogInfo("GenMETBlock") << "Too many GenMET, fnGenMET = " << fnGenMET; 
	  break;
        }
        genMetB = new ((*cloneGenMET)[fnGenMET++]) vhtm::GenMET();

        // fill in all the vectors
        genMetB->met    = it->pt();
        genMetB->metphi = it->phi();
        genMetB->sumet  = it->sumEt();
      }
    } 
    else {
      edm::LogError("GenMETBlock") << "Error >>  Failed to get GenMETCollection for label: " 
                                   << _inputTag;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenMETBlock);
