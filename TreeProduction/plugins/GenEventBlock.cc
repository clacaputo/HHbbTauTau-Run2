#include <algorithm>
#include <iostream>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TTree.h"
#include "TClonesArray.h"

#include "HHbbTauTau/TreeProduction/plugins/GenEventBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

GenEventBlock::GenEventBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  _genEvtInfoInputTag(iConfig.getParameter<edm::InputTag>("GenEventInfoInputTag")),
  _storePDFWeights(iConfig.getParameter<bool>("StorePDFWeights")),
  _pdfWeightsInputTag(iConfig.getParameter<edm::InputTag>("PDFWeightsInputTag"))
{}
GenEventBlock::~GenEventBlock() {
  delete _pdfWeights;
}
void GenEventBlock::beginJob() 
{
  _pdfWeights = new std::vector<double>();

  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneGenEvent = new TClonesArray("vhtm::GenEvent");
  tree->Branch("GenEvent", &cloneGenEvent, 32000, 2);
  tree->Branch("pdfWeights", "vector<double>", &_pdfWeights);
}
void GenEventBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray 
  cloneGenEvent->Clear();

  // Clear the two independent vectors
  _pdfWeights->clear();

  if (!iEvent.isRealData()) {
    // Create Event Object
    genEventB = new ( (*cloneGenEvent)[0] ) vhtm::GenEvent();

    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtInfoProduct;
    iEvent.getByLabel(_genEvtInfoInputTag, genEvtInfoProduct);
    if (genEvtInfoProduct.isValid()) {
      edm::LogInfo("GenEventBlock") << "Success >> obtained GenEventInfoProduct for label:" 
                                    << _genEvtInfoInputTag;
      genEventB->processID = genEvtInfoProduct->signalProcessID();
      genEventB->ptHat     = genEvtInfoProduct->hasBinningValues() 
                         ? genEvtInfoProduct->binningValues()[0] : 0.;
    } 
    else {
      edm::LogError("GenEventBlock") << "Error >> Failed to get GenEventInfoProduct for label: " 
                                     << _genEvtInfoInputTag;
    }
    // PDF Weights Part
    if (_storePDFWeights) {
      edm::Handle<std::vector<double> > pdfWeightsHandle;
      iEvent.getByLabel(_pdfWeightsInputTag, pdfWeightsHandle);
      if (pdfWeightsHandle.isValid()) {
	edm::LogInfo("GenEventBlock") << "Success >> obtained PDF handle for label: " 
                                      << _pdfWeightsInputTag;
	copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), genEventB->pdfWeights.begin());
	copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), _pdfWeights->begin());
      } 
      else {
	edm::LogError("GenEventBlock") << "Error >>  Failed to get PDF handle for label: " 
                                       << _pdfWeightsInputTag;
      }
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventBlock);
