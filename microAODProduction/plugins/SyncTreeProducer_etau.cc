/*!
 * \file SyncTreeProducer_mutau.cc
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
#include "HHbbTauTau/microAODProduction/interface/SyncTreeProducer_etau.h"

//
// constructors and destructor
//
SyncTreeProducer_etau::SyncTreeProducer_etau(const edm::ParameterSet& iConfig):
  BaseEDAnalyzer(iConfig),
  syncTree(&edm::Service<TFileService>()->file(),false),
  anaData("etau_cuts.root") {}


SyncTreeProducer_etau::~SyncTreeProducer_etau()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
SyncTreeProducer_etau::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace cuts::Htautau_2015;
  using namespace cuts::Htautau_2015::ETau;

  cuts::Cutter cut(&GetAnaData().Selection("events"));

  //Get collection
  Initialize(iEvent);

  try{

      cut(true,"events");

      //cut(HaveTriggerFired(iEvent),"trigger");

      const auto PV = (*vertices).ptrAt(0); //Deferenzio per passare da edm::Handle a edm::View. Quest'ultimo permette
                                            //di gestire una qualsiasi collezione del tipo passatogli tramite Tamplate.
                                            //Es. edm::View<int> gestisce int semplici, vector<int>, set<int> etc.
      cut(PV.isNonnull(),"vertex");

      const auto selectedMuons = CollectElectrons();

      cut(selectedMuons.size(),"electrons");
      //std::cout<<" Vertici    --->   "<<PV->ndof()<<std::endl;


  }catch(cuts::cut_failed&){}

  GetAnaData().Selection("events").fill_selection();

}


// ------------ method called once each job just before starting event loop  ------------
void
SyncTreeProducer_etau::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
SyncTreeProducer_etau::endJob()
{
    syncTree.Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void
SyncTreeProducer_mutau::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
SyncTreeProducer_mutau::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
SyncTreeProducer_mutau::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
SyncTreeProducer_mutau::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SyncTreeProducer_etau::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}


analysis::CandidateV2Ptr SyncTreeProducer_etau::SelectHiggs(analysis::CandidateV2PtrVector& higgses)
{
    if(!higgses.size())
        throw std::runtime_error("no available higgs candidate to select");
    if(higgses.size()==1) return higgses.front();

    const auto higgsSelector = [&] (const analysis::CandidateV2Ptr& first, const analysis::CandidateV2Ptr& second) -> bool
    {
        const pat::Muon& first_muon = first->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
        const pat::Tau&  first_tau  = first->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();
        const pat::Muon& second_muon = second->GetDaughter(analysis::CandidateV2::Type::Muon)->GetNtupleObject<pat::Muon>();
        const pat::Tau&  second_tau  = second->GetDaughter(analysis::CandidateV2::Type::Tau)->GetNtupleObject<pat::Tau>();

        double iso_mu1 = (first_muon.pfIsolationR03().sumChargedHadronPt + std::max(
                          first_muon.pfIsolationR03().sumNeutralHadronEt +
                          first_muon.pfIsolationR03().sumPhotonEt -
                          0.5 * first_muon.pfIsolationR03().sumPUPt, 0.0)) / first_muon.pt();
        double iso_mu2 = (second_muon.pfIsolationR03().sumChargedHadronPt + std::max(
                          second_muon.pfIsolationR03().sumNeutralHadronEt +
                          second_muon.pfIsolationR03().sumPhotonEt -
                          0.5 * second_muon.pfIsolationR03().sumPUPt, 0.0)) / second_muon.pt();

        double iso_tau1 = first_tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");
        double iso_tau2 = second_tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

        std::cout << "LAMBDA ------------------------  \n IsoMu1 = "<<iso_mu1<<"  IsoMu2 = "<<iso_mu2
                  << "  IsoTau1 = "<<iso_tau1<<"  IsoTau2 = "<<iso_tau2<<std::endl;
        std::cout << " PtMu1 = "<<first_muon.pt()<<"  PtMu2 = "<<second_muon.pt()
                  << "  PtTau1 = "<<first_tau.pt()<<"  PtTau2 = "<<second_tau.pt()<<std::endl;

        if (!iso_mu1 || !iso_mu2 ){
            std::cout << "\t 1st Muon \t ChargedHadronPt = "<<first_muon.pfIsolationR03().sumChargedHadronPt
                      << "\t NeutralHadronEt = "<<first_muon.pfIsolationR03().sumNeutralHadronEt
                      << "\t PhotonEt = "<<first_muon.pfIsolationR03().sumPhotonEt
                      << "\t PUPt = "<<first_muon.pfIsolationR03().sumPUPt<<std::endl;
            std::cout << "\t 2nd Muon \t ChargedHadronPt = "<<second_muon.pfIsolationR03().sumChargedHadronPt
                      << "\t NeutralHadronEt = "<<second_muon.pfIsolationR03().sumNeutralHadronEt
                      << "\t PhotonEt = "<<second_muon.pfIsolationR03().sumPhotonEt
                      << "\t PUPt = "<<second_muon.pfIsolationR03().sumPUPt<<std::endl;
        }
        bool muon = (iso_mu1 < iso_mu2) ||
                    ((iso_mu1 == iso_mu2) ? first_muon.pt() > second_muon.pt() : false);
        bool tau = (iso_tau1 < iso_tau2) ||
                    ((iso_tau1 == iso_tau2) ? first_tau.pt() > second_tau.pt() : false);

//        if ( &first_muon==&second_muon ){
//            if ( &first_tau==&second_tau ) return true;
//            return tau;
//        }

//        if (!muon) return tau;
//        return muon;
        if ( &first_muon!=&second_muon ) return muon;
        if ( &first_tau!=&second_tau ) return tau;

        throw analysis::exception("not found a good criteria for best tau pair");
    };

    std::sort(higgses.begin(), higgses.end(), higgsSelector) ;
    return higgses.front();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SyncTreeProducer_etau);