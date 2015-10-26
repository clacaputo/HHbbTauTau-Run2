/*!
 * \file SyncTreeProducer.cc
 * \author author: Claudio Caputo
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * X->HH->bbTauTau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.
 */

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/PatCandidates/interface/Tau.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//HHbbTauTau Framework
#include "HHbbTauTau/AnalysisBase/include/SyncTree.h"
#include "HHbbTauTau/TreeProduction/interface/Tau.h"
#include "HHbbTauTau/AnalysisBase/include/AnalyzerData.h"
#include "HHbbTauTau/AnalysisBase/include/CutTools.h"
#include "HHbbTauTau/Analysis/include/SelectionResults.h"

#include "TTree.h"
#include "Math/VectorUtil.h"


//Analyzer Data Class
//namespace analysis {
//class SyncAnalyzerData : public root_ext::AnalyzerData {
//public:
//    SyncAnalyzerData(TFile outputFile) {
//	std::shared_ptr<TFile> sp_outputFile;
//	sp_outputFile = std::shared_ptr<TFile>(new TFile(outputFile));
//	AnalyzerData(sp_outputFile);
//    }
//
//    SELECTION_ENTRY(Selection)
//
//};
//}

//
// class declaration
//

class SyncTreeProducer : public edm::EDAnalyzer {
   public:
      explicit SyncTreeProducer(const edm::ParameterSet&);
      ~SyncTreeProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum ElectronMatchType {UNMATCHED = 0,
              TRUE_PROMPT_ELECTRON,
              TRUE_ELECTRON_FROM_TAU,
              TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

//      analysis::SyncAnalyzerData& GetAnaData() { return anaData; }

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el,
               const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonElectronMother(const reco::Candidate *particle,
                    int &ancestorPID, int &ancestorStatus);

      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      // ... none ...

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
      edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

      //Tau Tag
      edm::EDGetToken tausMiniAODToken_;

      ntuple::SyncTree syncTree;

//      analysis::SyncAnalyzerData anaData;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SyncTreeProducer::SyncTreeProducer(const edm::ParameterSet& iConfig):
  tausMiniAODToken_(mayConsume<edm::View<pat::Tau> >(iConfig.getParameter<edm::InputTag>("tauSrc"))),
  syncTree(&edm::Service<TFileService>()->file(),false) // anaData(&edm::Service<TFileService>()->file())
{
  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticlesMiniAOD"));

}


SyncTreeProducer::~SyncTreeProducer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SyncTreeProducer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

 // cuts::Cutter cut(&GetAnaData().Selection("events"));

  // Save global info right away
  syncTree.run()  = iEvent.id().run();
  syncTree.lumi() = iEvent.id().luminosityBlock();
  syncTree.evt()  = iEvent.id().event();

  // Get the MC collection
//  iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  //Get Tau collection
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(tausMiniAODToken_, taus);
//Usare ntuple::Muon e Tau definiti in TreeProduction in modo da poter utilizzare i metodi del BaseAnalyzer
  ntuple::TauVector tausV;

 // try{

 //     cut(true,"events");
  for (const pat::Tau &tau : *taus){
      ntuple::Tau tmp_tau;

      //pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(src.leadChargedHadrCand().get());
      //fabs(packedLeadTauCand->dz()) < 0.2;  // The PackedCandidate::dz() method is wrt. the first PV by default

      if(!(tau.pt() > 20 && fabs(tau.eta()) < 2.3 && tau.tauID("decayModeFindingNewDMs") > 0.5)) continue;
      tmp_tau.eta = tau.eta();
      tmp_tau.pt  = tau.pt();
      tmp_tau.phi = tau.phi();
 //     tmp_tau.againstElectronLooseMVA5   = tau.tauID('againstElectronLooseMVA5');
 //     tmp_tau.againstElectronMediumMVA5  = tau.tauID('againstElectronMediumMVA5');
 //     tmp_tau.againstElectronTightMVA5   = tau.tauID('againstElectronTightMVA5');
 //     tmp_tau.againstElectronVTightMVA5  = tau.tauID('againstElectronVTightMVA5');

      tausV.push_back(tmp_tau);

  }

    //if(!tausV.size()) return;
 //   cut(tausV.size(),"taus");
    std::cout<<"Taus size:  "<<tausV.size()<<std::endl;
 // }catch(cuts::cut_failed&){}

  syncTree.pt_1() = tausV.at(0).pt;
  syncTree.Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void
SyncTreeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SyncTreeProducer::endJob()
{
    syncTree.Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void
SyncTreeProducer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
SyncTreeProducer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SyncTreeProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SyncTreeProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SyncTreeProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

int SyncTreeProducer::matchToTruth(const edm::Ptr<reco::GsfElectron> el,
                  const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  //
  // Explicit loop and geometric matching method (advised by Josh Bendavid)
  //

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  //
  int ancestorPID = -999;
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }

  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void SyncTreeProducer::findFirstNonElectronMother(const reco::Candidate *particle,
                         int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SyncTreeProducer);
