#ifndef __VHTauTau_TreeMaker_GenMETBlock_h
#define __VHTauTau_TreeMaker_GenMETBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class GenMET;

class GenMETBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit GenMETBlock(const edm::ParameterSet& iConfig);
  virtual ~GenMETBlock() {}

  enum {
    kMaxGenMET = 100
  };

private:
  TClonesArray* cloneGenMET; 
  int  fnGenMET;
  int _verbosity;
  edm::InputTag _inputTag;

  vhtm::GenMET* genMetB;
};
#endif
