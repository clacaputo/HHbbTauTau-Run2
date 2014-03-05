#ifndef __VHTauTau_TreeMaker_TriggerBlock_h
#define __VHTauTau_TreeMaker_TriggerBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class TClonesArray;
//class Trigger;

class TriggerBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup);
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit TriggerBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerBlock();

private:
  //TClonesArray* cloneTrigger; 

  int _verbosity;
  edm::InputTag _inputTag;

  const edm::InputTag _l1InputTag;
  const edm::InputTag _hltInputTag;
  const std::vector<std::string> _hltPathsOfInterest;
  HLTConfigProvider hltConfig;

  std::vector<int> *_l1physbits;
  std::vector<int> *_l1techbits;
  std::vector<std::string> *_hltpaths;
  std::vector<int> *_hltresults;
  std::vector<int> *_hltprescales;
};
#endif
