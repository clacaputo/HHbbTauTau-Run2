#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TPRegexp.h"

#include "HHbbTauTau/TreeProduction/plugins/TriggerBlock.h"
#include "HHbbTauTau/TreeProduction/interface/Utility.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"

static const unsigned int NmaxL1AlgoBit = 128;
static const unsigned int NmaxL1TechBit = 64;

// Constructor
TriggerBlock::TriggerBlock(const edm::ParameterSet& iConfig) :
  _verbosity(iConfig.getParameter<int>("verbosity")),
  _l1InputTag(iConfig.getParameter<edm::InputTag>("l1InputTag")),
  _hltInputTag(iConfig.getParameter<edm::InputTag>("hltInputTag")),
  _hltPathsOfInterest(iConfig.getParameter<std::vector<std::string> > ("hltPathsOfInterest"))
{
}
TriggerBlock::~TriggerBlock() {
  delete _l1physbits;
  delete _l1techbits;
  delete _hltpaths;
  delete _hltresults;
  delete _hltprescales;
}
void TriggerBlock::beginJob() 
{
  _l1physbits = new std::vector<int>();
  _l1techbits = new std::vector<int>();
  _hltpaths = new std::vector<std::string>();
  _hltresults = new std::vector<int>();
  _hltprescales = new std::vector<int>();

  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  tree->Branch("l1physbits", "vector<int>", &_l1physbits);
  tree->Branch("l1techbits", "vector<int>", &_l1techbits);
  tree->Branch("hltpaths", "vector<string>", &_hltpaths);
  tree->Branch("hltresults", "vector<int>", &_hltresults);
  tree->Branch("hltprescales", "vector<int>", &_hltprescales);
}
void TriggerBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig.init(iRun, iSetup, _hltInputTag.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerBlock") << "HLT config with process name " 
                                 << _hltInputTag.process() << " successfully extracted";
  } 
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerBlock") << "Error! HLT config extraction with process name " 
                                  << _hltInputTag.process() << " failed";
    // In this case, all access methods will return empty values!
  }
  
}
void TriggerBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the TClonesArray and the nObj variables
  //cloneTrigger->Clear();
  _l1physbits->clear();
  _l1techbits->clear();
  _hltpaths->clear();
  _hltresults->clear();
  _hltprescales->clear();

  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  iEvent.getByLabel(_l1InputTag, l1GtReadoutRecord);

  if (l1GtReadoutRecord.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: " 
                                 << _l1InputTag;
    for (unsigned int i = 0; i < NmaxL1AlgoBit; ++i) {
      _l1physbits->push_back(l1GtReadoutRecord->decisionWord()[i] ? 1 : 0);
    }
    for (unsigned int i = 0; i < NmaxL1TechBit; ++i) {
      _l1techbits->push_back(l1GtReadoutRecord->technicalTriggerWord()[i] ? 1 : 0 );
    }
  } 
  else {
    edm::LogError("TriggerBlock") << "Error >> Failed to get L1GlobalTriggerReadoutRecord for label: " 
                                  << _l1InputTag;
  }
  edm::Handle<edm::TriggerResults> triggerResults;
  iEvent.getByLabel(_hltInputTag, triggerResults);
  if (triggerResults.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained " << _hltInputTag;
    const std::vector<std::string>& pathList = hltConfig.triggerNames();
    for (std::vector<std::string>::const_iterator it = pathList.begin();
                                                 it != pathList.end(); ++it) {
      if (_hltPathsOfInterest.size()) {
        int nmatch = 0;
        for (std::vector<std::string>::const_iterator kt = _hltPathsOfInterest.begin();
                                                     kt != _hltPathsOfInterest.end(); ++kt) {
	  nmatch += TPRegexp(*kt).Match(*it);
        }
        if (!nmatch) continue; 
      }
      _hltpaths->push_back(*it);

      int fired = 0;
      unsigned int index = hltConfig.triggerIndex(*it);
      if (index < triggerResults->size()) {
        if (triggerResults->accept(index)) fired = 1;
      }
      else {
	edm::LogInfo("TriggerBlock") << "Requested HLT path \"" << (*it) << "\" does not exist";
      }
      _hltresults->push_back(fired);

      int prescale = -1;
      if (hltConfig.prescaleSet(iEvent, iSetup) < 0) {
	edm::LogError("TriggerBlock") << "Error >> The prescale set index number could not be obtained";
      } 
      else {
        prescale = hltConfig.prescaleValue(iEvent, iSetup, *it);
      }
      _hltprescales->push_back(prescale);

      if (_verbosity)
	std::cout << ">>> Path: " << (*it) 
                  << ", prescale: " << prescale 
                  << ", fired: " << fired
                  << std::endl; 
    }
  } else {
    edm::LogError("TriggerBlock") << "Error >> Failed to get TriggerResults for label: " 
                                  << _hltInputTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerBlock);
