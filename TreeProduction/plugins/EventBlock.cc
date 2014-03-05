#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "TTree.h"
#include "TClonesArray.h"

#include "HHbbTauTau/TreeProduction/plugins/EventBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

EventBlock::EventBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _l1InputTag(iConfig.getParameter<edm::InputTag>("l1InputTag")),
  _vtxInputTag(iConfig.getParameter<edm::InputTag>("vertexInputTag")),
  _trkInputTag(iConfig.getParameter<edm::InputTag>("trkInputTag")),
  _vtxMinNDOF(iConfig.getParameter<unsigned int>("vertexMinimumNDOF")),
  _vtxMaxAbsZ(iConfig.getParameter<double>("vertexMaxAbsZ")),
  _vtxMaxd0(iConfig.getParameter<double>("vertexMaxd0")),
  _numTracks(iConfig.getParameter<unsigned int>("numTracks")),
  _hpTrackThreshold(iConfig.getParameter<double>("hpTrackThreshold"))
{
}
EventBlock::~EventBlock() {
  delete _nPU;
  delete _bunchCrossing;
  delete _trueNInt;
}
void EventBlock::beginJob() {
  _nPU = new std::vector<int>();
  _bunchCrossing = new std::vector<int>();
  _trueNInt = new std::vector<int>();

  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  cloneEvent = new TClonesArray("vhtm::Event");
  tree->Branch("Event", &cloneEvent, 32000, 2);

  tree->Branch("nPU", "std::vector<int>", &_nPU);
  tree->Branch("bunchCrossing", "std::vector<int>", &_bunchCrossing);
  tree->Branch("trueNInt", "std::vector<int>", &_trueNInt);
}
void EventBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the TClonesArray
  cloneEvent->Clear();

  // Clear the two independent vectors
  _nPU->clear();
  _bunchCrossing->clear();
  _trueNInt->clear();

  // Create Event Object
  eventB = new ( (*cloneEvent)[0] ) vhtm::Event();
  eventB->run   = iEvent.id().run();
  eventB->event = iEvent.id().event();
  eventB->lumis = iEvent.id().luminosityBlock();
  eventB->bunch = iEvent.bunchCrossing();
  eventB->orbit = iEvent.orbitNumber();

  double sec  = iEvent.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & iEvent.time().value();
  double conv = 1e6;
  eventB->time = sec+usec/conv;
  eventB->isdata = iEvent.isRealData();
 
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(_l1InputTag, l1GtReadoutRecord);

  // Technical Trigger Part
  if (l1GtReadoutRecord.isValid()) {
    edm::LogInfo("EventBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: " 
                               << _l1InputTag;

    L1GtFdlWord fdlWord = l1GtReadoutRecord->gtFdlWord();
    if (fdlWord.physicsDeclared() == 1)
      eventB->isPhysDeclared = true;

    // BPTX0
    if (l1GtReadoutRecord->technicalTriggerWord()[0])
      eventB->isBPTX0 = true;

    // MinBias
    if (l1GtReadoutRecord->technicalTriggerWord()[40] || l1GtReadoutRecord->technicalTriggerWord()[41])
      eventB->isBSCMinBias = true;

    // BeamHalo
    if ( (l1GtReadoutRecord->technicalTriggerWord()[36] || l1GtReadoutRecord->technicalTriggerWord()[37] ||
          l1GtReadoutRecord->technicalTriggerWord()[38] || l1GtReadoutRecord->technicalTriggerWord()[39]) ||
         ((l1GtReadoutRecord->technicalTriggerWord()[42] && !l1GtReadoutRecord->technicalTriggerWord()[43]) ||
          (l1GtReadoutRecord->technicalTriggerWord()[43] && !l1GtReadoutRecord->technicalTriggerWord()[42])) )
      eventB->isBSCBeamHalo = true;
  } 
  else {
    edm::LogError("EventBlock") << "Error >> Failed to get L1GlobalTriggerReadoutRecord for label:" 
                                << _l1InputTag;
    throw std::runtime_error("Failed to get L1GlobalTriggerReadoutRecord for label.");
  }

  // Good Primary Vertex Part
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(_vtxInputTag, primaryVertices);

  if (primaryVertices.isValid()) {
    edm::LogInfo("EventBlock") << "Total # Primary Vertices: " << primaryVertices->size();
    for (reco::VertexCollection::const_iterator it = primaryVertices->begin(); 
                                               it != primaryVertices->end(); ++it) {
      if (!(it->isFake()) && it->ndof() > _vtxMinNDOF &&
            fabs(it->z()) <= _vtxMaxAbsZ && fabs(it->position().rho()) <= _vtxMaxd0
	  ) 
      {
        eventB->isPrimaryVertex = true;
        break;
      }
    }
  } 
  else {
    edm::LogError("EventBlock") << "Error >> Failed to get VertexCollection for label:" 
                                << _vtxInputTag;
    throw std::runtime_error("Failed to get VertexCollection for label.");
  }

  // Scraping Events Part
  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByLabel(_trkInputTag, tracks);

  if (tracks.isValid()) {
    edm::LogInfo("EventBlock") << "Total # Tracks: " << tracks->size();

    int numhighpurity = 0;
    double fraction   = 1.;
    reco::TrackBase::TrackQuality trackQuality = reco::TrackBase::qualityByName("highPurity");
    if (tracks->size() > _numTracks) {
      for (reco::TrackCollection::const_iterator it = tracks->begin(); 
                                                it != tracks->end(); ++it) {
        if (it->quality(trackQuality)) numhighpurity++;
      }
      fraction = (double)numhighpurity/(double)tracks->size();
      if (fraction < _hpTrackThreshold)
        eventB->isBeamScraping = true;
    }
  } 
  else {
    edm::LogError("EventBlock") << "Error >> Failed to get TrackCollection for label:" 
                                << _trkInputTag;
    throw std::runtime_error("Failed to get TrackCollection for label");
  }
  // Access PU information
  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByLabel("addPileupInfo", PupInfo);

    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for (PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      eventB->bunchCrossing.push_back(PVI->getBunchCrossing());
      _bunchCrossing->push_back(PVI->getBunchCrossing());
      eventB->nPU.push_back(PVI->getPU_NumInteractions());      
      _nPU->push_back(PVI->getPU_NumInteractions());      
      eventB->trueNInt.push_back(PVI->getTrueNumInteractions());      
      _trueNInt->push_back(PVI->getTrueNumInteractions());    
    }

    // More info about PU is here:
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
  } 
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventBlock);
