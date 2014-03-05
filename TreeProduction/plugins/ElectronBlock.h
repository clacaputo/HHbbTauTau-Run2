
#ifndef __VHTauTau_TreeMaker_ElectronBlock_h
#define __VHTauTau_TreeMaker_ElectronBlock_h

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

#include "HHbbTauTau/TreeProduction/interface/PhysicsObjects.h"

class TClonesArray;
class Electron;

class ElectronBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit ElectronBlock(const edm::ParameterSet& iConfig);
  virtual ~ElectronBlock();

  enum {
    kMaxElectron = 100
  };

private:
  TClonesArray* cloneElectron; 
  int fnElectron;

  int _verbosity;
  edm::InputTag _bsInputTag;
  edm::InputTag _trkInputTag;
  edm::InputTag _vtxInputTag;
  edm::InputTag _electronInputTag;

  vhtm::Electron* electronB;
};
#endif
