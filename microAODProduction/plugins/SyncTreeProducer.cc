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
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

//Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"


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
#include "../../AnalysisBase/include/SyncTree.h"
#include "../../TreeProduction/interface/Tau.h"
#include "../../TreeProduction/interface/Muon.h"
#include "../../AnalysisBase/include/AnalyzerData.h"
#include "../../AnalysisBase/include/CutTools.h"
#include "../interface/Candidate.h"
#include "../../Analysis/include/SelectionResults.h"

#include "../interface/Htautau_2015.h"

#include "TTree.h"
#include "Math/VectorUtil.h"

using namespace analysis;

//Analyzer Data Class
namespace analysis {
class SyncAnalyzerData : public root_ext::AnalyzerData {
public:
    SyncAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}

    SELECTION_ENTRY(Selection)

    TH1D_ENTRY_FIX(N_objects, 1, 500, -0.5)
    TH1D_ENTRY(Mass, 3000, 0.0, 3000.0)

};
}

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

      analysis::SyncAnalyzerData& GetAnaData() { return anaData; }

      analysis::CandidateV2PtrVector FindCompatibleObjects(const CandidateV2PtrVector& objects1,
                                                         const CandidateV2PtrVector& objects2,
                                                         double minDeltaR, CandidateV2::Type type, const std::string& hist_name);

      analysis::CandidateV2PtrVector ApplyTriggerMatch(
             const pat::TriggerObjectStandAloneCollection& triggerObjects,
                                                       const CandidateV2PtrVector& higgses,const edm::TriggerNames &names,
                                                       const std::set<std::string>& hltPaths,
                                                     bool useStandardTriggerMatch);

      int matchToTruth(const edm::Ptr<reco::GsfElectron> el,
               const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);

      void findFirstNonElectronMother(const reco::Candidate *particle,
                    int &ancestorPID, int &ancestorStatus);



      inline bool HaveTriggerMatched(const pat::TriggerObjectStandAloneCollection& triggerObjects,
                                     const edm::TriggerNames &names,
                                     const std::set<std::string>& interestingPaths, const CandidateV2& candidate,
                                     double deltaR_Limit)
      {
          for (const std::string& interestinPath : interestingPaths){
              if (HaveTriggerMatched(triggerObjects,names,interestinPath,candidate, deltaR_Limit)) return true;
          }
          return false;
      }

      inline bool HaveTriggerMatched(const pat::TriggerObjectStandAloneCollection& triggerObjects,
                                     const edm::TriggerNames &names,
                                     const std::string& interestingPath,
                                     const CandidateV2& candidate, double deltaR_Limit)
      {
          if(candidate.GetFinalStateDaughters().size()) {
              for(const auto& daughter : candidate.GetFinalStateDaughters()) {
                  if(!SyncTreeProducer::HaveTriggerMatched(triggerObjects, names, interestingPath, *daughter, deltaR_Limit))
                      return false;
              }
              return true;
          }

          for (auto &triggerObject : triggerObjects){
              triggerObject.unpackPathNames(names);
              TLorentzVector triggerObjectMomentum;
              triggerObjectMomentum.SetPtEtaPhiM(triggerObject.pt(), triggerObject.eta(), triggerObject.phi(), triggerObject.mass());
              for (unsigned n = 0; n < triggerObject.pathNames(false).size(); ++n){
                  const std::string& objectMatchedPath = triggerObject.pathNames(false).at(n);
                  const size_t found = objectMatchedPath.find(interestingPath);
                  bool isBoth = triggerObject.hasPathName( triggerObject.pathNames(false).at(n), true, true );
                  if (found != std::string::npos && isBoth &&
                          triggerObjectMomentum.DeltaR(candidate.GetMomentum()) < deltaR_Limit)
                      return true;
              }
          }
          return false;
      }

    //  bool HaveTriggerFired(const std::set<std::string>& interestinghltPaths) const ;

      // ----------member data ---------------------------

      // Data members that are the same for AOD and miniAOD
      // ... none ...

      // MiniAOD case data members
      edm::EDGetToken electronsMiniAODToken_;
  //    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesMiniAODToken_;

      // ID decisions objects
      edm::EDGetTokenT<edm::ValueMap<bool> > eleMediumIdMapToken_;
      edm::EDGetTokenT<edm::ValueMap<bool> > eleTightIdMapToken_;

      // MVA values and categories (optional)
      edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
      edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;

      //InputTag
      edm::EDGetToken tausMiniAODToken_;
      edm::EDGetToken muonsMiniAODToken_;
      edm::EDGetToken vtxMiniAODToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
      edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;

      ntuple::SyncTree syncTree;
      analysis::SyncAnalyzerData anaData;
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
  muonsMiniAODToken_(mayConsume<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muonSrc"))),
  vtxMiniAODToken_(mayConsume<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vtxSrc"))),
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  syncTree(&edm::Service<TFileService>()->file(),false),
  anaData("cuts.root")
{
//  genParticlesMiniAODToken_ = mayConsume<edm::View<reco::GenParticle> >
//    (iConfig.getParameter<edm::InputTag>
//     ("genParticlesMiniAOD"));

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
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::MuTau;

  cuts::Cutter cut(&GetAnaData().Selection("events"));

  // Save global info right away
  syncTree.run()  = iEvent.id().run();
  syncTree.lumi() = iEvent.id().luminosityBlock();
  syncTree.evt()  = iEvent.id().event();

  // Get the MC collection
//  iEvent.getByToken(genParticlesMiniAODToken_,genParticles);

  //Get collection
  edm::Handle<edm::View<pat::Tau> > taus;
  iEvent.getByToken(tausMiniAODToken_, taus);
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(muonsMiniAODToken_, muons);
  edm::Handle<edm::View<reco::Vertex> > vertexes;
  iEvent.getByToken(vtxMiniAODToken_, vertexes);
  edm::Handle<edm::TriggerResults> triggerBits;
  iEvent.getByToken(triggerBits_, triggerBits);
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjects_, triggerObjects);




   /*   std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
          for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
              obj.unpackPathNames(names);
              std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
              // Print trigger object collection and type
              std::cout << "\t   Collection: " << obj.collection() << std::endl;
              std::cout << "\t   Type IDs:   ";
              for (unsigned h = 0; h < obj.filterIds().size(); ++h) std::cout << " " << obj.filterIds()[h] ;
              std::cout << std::endl;
              // Print associated trigger filters
              std::cout << "\t   Filters:    ";
              for (unsigned h = 0; h < obj.filterLabels().size(); ++h) std::cout << " " << obj.filterLabels()[h];
              std::cout << std::endl;
              std::vector<std::string> pathNamesAll  = obj.pathNames(false);
              std::vector<std::string> pathNamesLast = obj.pathNames(true);
              // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
              // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
              // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
              std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
              for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
                  bool isBoth = obj.hasPathName( pathNamesAll[h], true, true );
                  bool isL3   = obj.hasPathName( pathNamesAll[h], false, true );
                  bool isLF   = obj.hasPathName( pathNamesAll[h], true, false );
                  bool isNone = obj.hasPathName( pathNamesAll[h], false, false );
                  std::cout << "   " << pathNamesAll[h];
                  if (isBoth) std::cout << "(L,3)";
                  if (isL3 && !isBoth) std::cout << "(*,3)";
                  if (isLF && !isBoth) std::cout << "(L,*)";
                  if (isNone && !isBoth && !isL3 && !isLF) std::cout << "(*,*)";
              }
              std::cout << std::endl;
          }
          std::cout << std::endl;

          */

  //Usare ntuple::Muon e Tau definiti in TreeProduction in modo da poter utilizzare i metodi del BaseAnalyzer
  ntuple::TauVector tausV;
  ntuple::MuonVector muonsV;

//  std::map <size_t, ntuple::Muon> muonsMaps;

  try{

      cut(true,"events");

      bool triggerFired = false;
      const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    //  const auto Key = analysis::DataSourceType::Spring15MC;
    //  const auto& hltPaths = MuTau::trigger::hltPathMaps[Key];
      const auto& hltPaths = MuTau::trigger::hltPathMC;
          std::cout << "\n === TRIGGER PATHS === " << std::endl;
          for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
              for (const std::string& triggerPath : hltPaths ){

                  const std::string& objectMatchedPath = names.triggerName(i);
                  size_t found = objectMatchedPath.find(triggerPath);


                  if(found != std::string::npos){
//                      std::cout << "TriggerPath --->   " << triggerPath << std::endl;
//                      std::cout << "Trigger " << names.triggerName(i) <<
//                      ", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
//                      ": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)")
//                      << std::endl;
                      if(triggerPrescales->getPrescaleForIndex(i) == 1 && triggerBits->accept(i)) triggerFired = true;
                }
              }
          }

      cut(triggerFired,"trigger");

      const auto PV = (*vertexes).ptrAt(0); //Deferenzio per passare da edm::Handle a edm::View. Quest'ultimo permette
                                            //di gestire una qualsiasi collezione del tipo passatogli tramite Tamplate.
                                            //Es. edm::View<int> gestisce int semplici, vector<int>, set<int> etc.
      cut(PV.isNonnull(),"vertex");

      //std::cout<<" Vertici    --->   "<<PV->ndof()<<std::endl;

      std::vector<pat::Muon> patMuonsVector;
      CandidateV2PtrVector muonCollection, tauCollection;
      for(const pat::Muon &muon : *muons){
        ntuple::Muon tmp_muon;


        if(!(muon.pt() > muonID::pt  &&
            fabs(muon.eta()) < muonID::eta &&
            muon.isMediumMuon() == muonID::isMediumMuon)) continue;

        bool muonIP = fabs(muon.muonBestTrack()->dxy(PV->position())) < muonID::dB &&
                      fabs(muon.muonBestTrack()->dz(PV->position())) < muonID::dz;
        if(!muonIP) continue;

        tmp_muon.eta = muon.eta();

        const CandidateV2Ptr tmp_candidate(new CandidateV2(muon));
        muonCollection.push_back(tmp_candidate);
        patMuonsVector.push_back(muon);
        muonsV.push_back(tmp_muon);

      }

      cut(muonCollection.size(),"muons");

      const CandidateV2Ptr candidate(new CandidateV2(patMuonsVector.at(0)));
      std::cout<< "Muons Checks ------------------------------------------------ \n"
                  << "Candidate Muon :  "<< candidate->GetMomentum()<<"\n"
                     << "MuonVector  :  "<< (muonCollection.at(0))->GetMomentum()<<"\n";
//                  <<"Eta vettore  :  "<< patMuonsVector.at(0).eta() <<  "  type of muon:  "<<
//                    typeid(patMuonsVector.at(0)).name() << "\n"
//                    << typeid(pat::Muon).name() << "\n"
//                    << "Hash COde : "<<typeid(patMuonsVector.at(0)).hash_code() <<"\n"
//                    << typeid(pat::Muon).hash_code() << "\n"
//                    <<"   type ntuple  :  "<< typeid(muonsV.at(0)).name() << std::endl;

      for (const pat::Tau &tau : *taus){

          ntuple::Tau tmp_tau;

      pat::PackedCandidate const* packedLeadTauCand =
              dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      //fabs(packedLeadTauCand->dz()) < 0.2;  // The PackedCandidate::dz() method is wrt. the first PV by default


          if(!(tau.pt() > tauID::pt && fabs(tau.eta()) < tauID::eta &&
               tau.tauID("decayModeFindingNewDMs") > tauID::decayModeFinding &&
               fabs(packedLeadTauCand->dz()) < tauID::dz )) continue;


        const CandidateV2Ptr tmp_candidate(new CandidateV2(tau));
        tmp_tau.eta = tau.eta();
//      tmp_tau.pt  = tau.pt();
//      tmp_tau.phi = tau.phi();
 //     tmp_tau.againstElectronLooseMVA5   = tau.tauID('againstElectronLooseMVA5');
 //     tmp_tau.againstElectronMediumMVA5  = tau.tauID('againstElectronMediumMVA5');
 //     tmp_tau.againstElectronTightMVA5   = tau.tauID('againstElectronTightMVA5');
 //     tmp_tau.againstElectronVTightMVA5  = tau.tauID('againstElectronVTightMVA5');

        tauCollection.push_back(tmp_candidate);
          tausV.push_back(tmp_tau);
      }

    cut(tauCollection.size(),"taus");

    auto higgses = FindCompatibleObjects(muonCollection,tauCollection,DeltaR_betweenSignalObjects,CandidateV2::Type::Higgs,
                                         "h_mu_tau");

    cut(higgses.size(),"mu_tau");

    auto triggeredHiggses = ApplyTriggerMatch(triggerObjects.product(), higgses,names,MuTau::trigger::hltPathMC,false);

    cut(triggeredHiggses.size(),"triggerMatch");

//  syncTree.pt_1() = tausV.at(0).pt;
//  syncTree.Fill();

  }catch(cuts::cut_failed&){}

  GetAnaData().Selection("events").fill_selection();

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

///////////


analysis::CandidateV2PtrVector SyncTreeProducer::FindCompatibleObjects(const CandidateV2PtrVector& objects1,
                                                   const CandidateV2PtrVector& objects2,
                                                   double minDeltaR, CandidateV2::Type type, const std::string& hist_name)
    {
        CandidateV2PtrVector result;
        for(const CandidateV2Ptr& object1 : objects1) {
            for(const CandidateV2Ptr& object2 : objects2) {
                if(object2->GetMomentum().DeltaR(object1->GetMomentum()) > minDeltaR) {
                    const CandidateV2Ptr candidate(new CandidateV2(type, object1, object2));
//                    if (expectedCharge != CandidateV2::UnknownCharge() && candidate->GetCharge() != expectedCharge)
//                        continue;
                    result.push_back(candidate);
                    GetAnaData().Mass(hist_name).Fill(candidate->GetMomentum().M(),1);
                }
            }
        }

        GetAnaData().N_objects(hist_name).Fill(result.size(),1);
        return result;
    }

CandidateV2PtrVector SyncTreeProducer::ApplyTriggerMatch(const pat::TriggerObjectStandAloneCollection& triggerObjects,
                                                         const CandidateV2PtrVector& higgses, const edm::TriggerNames &names,
                                                         const std::set<std::string>& hltPaths,
                                                         bool useStandardTriggerMatch)
{
    CandidateV2PtrVector triggeredHiggses;
    for (const auto& higgs : higgses){
        if(!useStandardTriggerMatch && SyncTreeProducer::HaveTriggerMatched(triggerObjects,names, hltPaths, *higgs,
                                                                    cuts::Htautau_2015::DeltaR_triggerMatch))
            triggeredHiggses.push_back(higgs);
//        if (useStandardTriggerMatch && SyncTreeProducer::HaveTriggerMatched(hltPaths, *higgs))
//            triggeredHiggses.push_back(higgs);
    }
    return triggeredHiggses;
}



////Forse bisogna passargli anche iEvent,
//bool HaveTriggerFired(const std::set<std::string>& interestinghltPaths) const
//{
//    for (const ntuple::Trigger& trigger : event->triggers()){
//        for (size_t n = 0; HaveTriggerMatched(trigger.hltpaths, interestinghltPaths, n); ++n){
//            if (trigger.hltresults.at(n) == 1 && trigger.hltprescales.at(n) == 1)
//                return true;
//        }
//    }
//    return false;
//}

////Far diventare questa una Lambda
//inline bool HaveTriggerMatched(const std::vector<std::string>& objectMatchedPaths,
//                               const std::set<std::string>& interestinghltPaths, size_t& n)
//{
//    for (; n < objectMatchedPaths.size(); ++n){
//        for (const std::string& interestingPath : interestinghltPaths){
//            const std::string& objectMatchedPath = objectMatchedPaths.at(n);
//            const size_t found = objectMatchedPath.find(interestingPath);
//            if (found != std::string::npos) return true;
//        }
//    }
//    return false;
//}

//define this as a plug-in
DEFINE_FWK_MODULE(SyncTreeProducer);
