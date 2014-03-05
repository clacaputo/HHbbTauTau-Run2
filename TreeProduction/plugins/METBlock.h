#ifndef __VHTauTau_TreeMaker_METBlock_h
#define __VHTauTau_TreeMaker_METBlock_h

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
class MET;

class METBlock : public edm::EDAnalyzer 
{
private:
  virtual void beginJob();
  virtual void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {}
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  virtual void endJob() {}

public:
  explicit METBlock(const edm::ParameterSet& iConfig);
  virtual ~METBlock() {}

  enum {
    kMaxMET = 5
  };

  void fillMET(const edm::Event& iEvent, 
               const edm::EventSetup& iSetup,
               TClonesArray* cloneMET,
               int fnMET,
               edm::InputTag iTag,
               vhtm::MET* metB);

private:
  TClonesArray* clonePFMET; 
  int  fnPFMET;

  TClonesArray* cloneCorrMET; 
  int  fnCorrMET;

  TClonesArray* cloneMVAMET; 
  int  fnMVAMET;

  int _verbosity;
  edm::InputTag _pfinputTag;
  edm::InputTag _corrinputTag;
  edm::InputTag _mvainputTag;

  vhtm::MET* pfmetB;
  vhtm::MET* corrmetB;
  vhtm::MET* mvametB;
};
#endif
