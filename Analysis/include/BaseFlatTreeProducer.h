/*!
 * \file BaseFlatTreeProducer.h
 * \brief Definition of BaseFlatTreeProducer class, the base class for flat tree producers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-07-11 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>
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

#pragma once

#include "BaseAnalyzer.h"
#include "FlatTree.h"

namespace analysis {

class BaseFlatTreeProducer : public BaseAnalyzer {
public:
    BaseFlatTreeProducer(const std::string& inputFileName, const std::string& outputFileName,
                         const std::string& configFileName, const std::string& _prefix = "none",
                         size_t _maxNumberOfEvents = 0,
                         std::shared_ptr<ntuple::FlatTree> _flatTree = std::shared_ptr<ntuple::FlatTree>())
        : BaseAnalyzer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          flatTree(_flatTree), writeFlatTree(!flatTree)
    {
        if(!flatTree)
            flatTree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree("flatTree"));
    }

    virtual ~BaseFlatTreeProducer() override
    {
        if(writeFlatTree)
            flatTree->Write();
    }

protected:

    virtual void CalculateTriggerWeights(const Candidate& candidate) override
    {
        triggerWeights.clear();
        triggerWeights.push_back(1);
        triggerWeights.push_back(1);
    }

    virtual void CalculateIsoWeights(const Candidate& candidate) override
    {
        IsoWeights.clear();
        IsoWeights.push_back(1);
        IsoWeights.push_back(1);
    }

    virtual void CalculateIdWeights(const Candidate& candidate) override
    {
        IDweights.clear();
        IDweights.push_back(1);
        IDweights.push_back(1);
    }

    virtual void CalculateDMWeights(const Candidate& candidate) override
    {
        DMweights.clear();
        DMweights.push_back(1);
        DMweights.push_back(1);
    }

    virtual void CalculateFakeWeights(const Candidate& candidate) override
    {
        fakeWeights.clear();
        fakeWeights.push_back(1);
        fakeWeights.push_back(1);
    }

    CandidateVector CollectBJets(const CandidateVector& looseJets, bool doReTag)
    {
        using namespace cuts::Htautau_Summer13::btag;
        analysis::CandidateVector bjets;
        for(const Candidate& looseJetCandidate : looseJets) {
            const ntuple::Jet& looseJet = event->jets().at(looseJetCandidate.index);
            if(looseJet.pt <= pt || std::abs(looseJet.eta) >= eta) continue;
            if(doReTag && !analysis::btag::ReTag(looseJet, btag::payload::EPS13, btag::tagger::CSVM, 0, 0, CSV))
                continue;
            else if(!doReTag && looseJet.combinedSecondaryVertexBJetTags <= CSV)
                continue;

            bjets.push_back(looseJetCandidate);
        }

        const auto bjetsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const ntuple::Jet& first_bjet = event->jets().at(first.index);
            const ntuple::Jet& second_bjet = event->jets().at(second.index);

            return first_bjet.combinedSecondaryVertexBJetTags > second_bjet.combinedSecondaryVertexBJetTags;
        };

        std::sort(bjets.begin(), bjets.end(), bjetsSelector);
        return bjets;
    }

    CandidateVector CollectLooseJets()
    {
        const auto base_selector = [&](unsigned id, cuts::ObjectSelector* _objectSelector,
                root_ext::AnalyzerData& _anaData, const std::string& _selection_label) -> Candidate
            { return SelectLooseJet(id, _objectSelector, _anaData, _selection_label); };
        return CollectObjects<Candidate>("jets", base_selector, event->jets().size());
    }

    virtual Candidate SelectLooseJet(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData, const std::string& selection_label)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        cuts::Cutter cut(objectSelector);
        const ntuple::Jet& object = event->jets().at(id);

        cut(true, ">0 jet cand");
        cut(X(pt) > pt_loose, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const bool passLooseID = analysis::passPFLooseId(object);
        cut(Y(passLooseID) == pfLooseID, "pfLooseID");
        const bool passPUlooseID = ntuple::JetID_MVA::PassLooseId(object.puIdBits);
        cut(Y(passPUlooseID) == puLooseID, "puLooseID");

        return Candidate(analysis::Candidate::Jet, id, object);
    }

    CandidateVector CollectJets(const CandidateVector& looseJets)
    {
        using namespace cuts::Htautau_Summer13::jetID;
        analysis::CandidateVector jets;
        for(const Candidate& looseJet : looseJets) {
            if(looseJet.momentum.Pt() > pt)
                jets.push_back(looseJet);
        }
        return jets;
    }

    bool FindAnalysisFinalState(analysis::finalState::TauTau& final_state)
    {
        static const particles::ParticleCodes resonanceCodes = { particles::Higgs, particles::Z };
        static const particles::ParticleCodes resonanceDecay = { particles::tau, particles::tau };

        genEvent.Initialize(event->genParticles());
//        genEvent.Print();

        const analysis::GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);

        analysis::GenParticlePtrVector resonancesToTauTau;

        for (const GenParticle* resonance : resonances){
            analysis::GenParticlePtrVector resonanceDecayProducts;
            if(analysis::FindDecayProducts(*resonance, resonanceDecay,resonanceDecayProducts)){
                final_state.taus = resonanceDecayProducts;
                final_state.resonance = resonance;
                resonancesToTauTau.push_back(resonance);
            }
        }

        if (resonancesToTauTau.size() > 1)
            throw exception("more than one resonance to tautau per event");

        if(resonancesToTauTau.size() == 0){
            if(config.ExpectedOneResonanceToTauTau())
                throw exception("resonance to tautau not found");
            else
                return false;
        }

        for (const analysis::GenParticle* tau_MC : final_state.taus) {
            analysis::GenParticlePtrVector TauProducts;
            if (!analysis::FindDecayProducts(*tau_MC, analysis::TauMuonicDecay, TauProducts)
                    && !analysis::FindDecayProducts(*tau_MC, analysis::TauElectronDecay, TauProducts))
                final_state.hadronic_taus.push_back(tau_MC);
        }
        return true;
    }

    void FillFlatTree(const Candidate& higgs      , double m_sv,
                      const CandidateVector& jets , const CandidateVector& jetsPt20,
                      const CandidateVector& bjets, const analysis::CandidateVector& retagged_bjets,
                      const VertexVector& vertices, const Candidate& leg1, const Candidate& leg2,
                      const ntuple::MET& pfMET)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();
        flatTree->run() = event->eventInfo().run;
        flatTree->lumi() = event->eventInfo().lumis;
        flatTree->evt() = event->eventInfo().EventId;

        flatTree->npv() = vertices.size();
        if (config.ApplyPUreweight()){
            const size_t bxIndex = tools::find_index(event->eventInfo().bunchCrossing, 0);
            if(bxIndex >= event->eventInfo().bunchCrossing.size())
                throw std::runtime_error("in-time BX not found");
            flatTree->npu() = event->eventInfo().trueNInt.at(bxIndex);
        }
        //flatTree->rho();

        //flatTree->mcweight();
        flatTree->puweight() = PUweight;
        flatTree->trigweight_1() = triggerWeights.at(0);
        flatTree->trigweight_2() = triggerWeights.at(1);
        flatTree->idweight_1() = IDweights.at(0);
        flatTree->idweight_2() = IDweights.at(1);
        flatTree->isoweight_1() = IsoWeights.at(0);
        flatTree->isoweight_2() = IsoWeights.at(1);
        //flatTree->effweight();
        flatTree->weight() = eventWeight;
        //flatTree->embeddedWeight();
        //flatTree->signalWeight();

        flatTree->mvis() = higgs.momentum.M();
        flatTree->m_sv() = m_sv;

        flatTree->pt_1()     = leg1.momentum.Pt();
        flatTree->phi_1()    = leg1.momentum.Phi();
        flatTree->eta_1()    = leg1.momentum.Eta();
        flatTree->m_1()      = leg1.momentum.M();
        flatTree->energy_1() = leg1.momentum.E();
        flatTree->q_1()      = leg1.charge;
        flatTree->mt_1()     = analysis::Calculate_MT(leg1.momentum, postRecoilMET.pt, postRecoilMET.phi);
        flatTree->d0_1()     = analysis::Calculate_dxy(leg1.vertexPosition, primaryVertex.position,leg1.momentum);
        flatTree->dZ_1()     = leg1.vertexPosition.Z() - primaryVertex.position.Z();

        // leg1 lepton specific variable should be filled outside. Here all them set to the default value.
        flatTree->pfRelIso_1()                                 = default_value;
        flatTree->decayMode_1()                                = default_value;
        flatTree->againstElectronLooseMVA_1()                  = default_value;
        flatTree->againstElectronMediumMVA_1()                 = default_value;
        flatTree->againstElectronTightMVA_1()                  = default_value;
        flatTree->againstElectronLoose_1()                     = default_value;
        flatTree->againstElectronMedium_1()                    = default_value;
        flatTree->againstElectronTight_1()                     = default_value;
        flatTree->againstMuonLoose_1()                         = default_value;
        flatTree->againstMuonMedium_1()                        = default_value;
        flatTree->againstMuonTight_1()                         = default_value;
        flatTree->passid_1()                                   = default_value;
        flatTree->passiso_1()                                  = default_value;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = default_value;

        flatTree->pt_2()     = leg2.momentum.Pt();
        flatTree->phi_2()    = leg2.momentum.Phi();
        flatTree->eta_2()    = leg2.momentum.Eta();
        flatTree->m_2()      = leg2.momentum.M();
        flatTree->energy_2() = leg2.momentum.E();
        flatTree->q_2()      = leg2.charge;
        flatTree->mt_2()     = analysis::Calculate_MT(leg2.momentum, postRecoilMET.pt, postRecoilMET.phi);
        flatTree->d0_2()     = analysis::Calculate_dxy(leg2.vertexPosition, primaryVertex.position,leg2.momentum);
        flatTree->dZ_2()     = leg2.vertexPosition.Z() - primaryVertex.position.Z();

        // RM: for the three channels, mt, et, tt this leg is always a tau
        const ntuple::Tau& ntuple_tau_leg2 = correctedTaus.at(leg2.index);
        flatTree->decayMode_2()                                = ntuple_tau_leg2.decayMode                               ;
//  needs custom definition
//         flatTree->againstElectronLooseMVA_2()                  = ntuple_tau_leg2.default_value;
//         flatTree->againstElectronMediumMVA_2()                 = ntuple_tau_leg2.default_value;
//         flatTree->againstElectronTightMVA_2()                  = ntuple_tau_leg2.default_value;
        flatTree->againstElectronLoose_2()                     = ntuple_tau_leg2.againstElectronLoose                    ;
        flatTree->againstElectronMedium_2()                    = ntuple_tau_leg2.againstElectronMedium                   ;
        flatTree->againstElectronTight_2()                     = ntuple_tau_leg2.againstElectronTight                    ;
        flatTree->againstMuonLoose_2()                         = ntuple_tau_leg2.againstMuonLoose                        ;
        flatTree->againstMuonMedium_2()                        = ntuple_tau_leg2.againstMuonMedium                       ;
        flatTree->againstMuonTight_2()                         = ntuple_tau_leg2.againstMuonTight                        ;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits;

        TLorentzVector postRecoilMetMomentum;
        postRecoilMetMomentum.SetPtEtaPhiM(postRecoilMET.pt, 0, postRecoilMET.phi, 0.);
        flatTree->pt_tt() = (leg1.momentum + leg2.momentum + postRecoilMetMomentum).Pt();

        flatTree->met() = pfMET.pt;
        flatTree->metphi() = pfMET.phi;
        flatTree->mvamet() = postRecoilMET.pt;
        flatTree->mvametphi() = postRecoilMET.phi;
        //flatTree->pzetavis();
        //flatTree->pzetamiss();
        if(pfMET.significanceMatrix.size()) {
            const TMatrixD metPFcov = ntuple::VectorToSignificanceMatrix(pfMET.significanceMatrix);
            flatTree->metcov00() = metPFcov[0][0];
            flatTree->metcov01() = metPFcov[0][1];
            flatTree->metcov10() = metPFcov[1][0];
            flatTree->metcov11() = metPFcov[1][1];
        }
        const TMatrixD metMVAcov = ntuple::VectorToSignificanceMatrix(postRecoilMET.significanceMatrix);
        flatTree->mvacov00() = metMVAcov[0][0];
        flatTree->mvacov01() = metMVAcov[0][1];
        flatTree->mvacov10() = metMVAcov[1][0];
        flatTree->mvacov11() = metMVAcov[1][1];

        flatTree->njets() = jets.size();

        flatTree->nBjets() = retagged_bjets.size();

        // RM: EXTRA MUONS
        // need to clean the muons from the signal muon!
        const auto muons_bkg = CollectBackgroundMuons() ;
        CandidateVector filtered_muons_bkg ;
        for (unsigned int imuon = 0 ; imuon < muons_bkg.size() ; ++imuon) {
          // for all analyses, leg1 is always the lepton
          if ( leg1.momentum.DeltaR(muons_bkg.at(imuon).momentum) < 0.01 ) continue ;
          filtered_muons_bkg.push_back(muons_bkg.at(imuon)) ;
        }
        // sort the muons by pt
        std::sort( filtered_muons_bkg.begin(), filtered_muons_bkg.end(), []( analysis::Candidate a, analysis::Candidate b ){ return a.momentum.Pt() > b.momentum.Pt(); } ) ;

        // fill the extra muons collection
        for (unsigned int imuon = 0 ; imuon <= 3 ; ++imuon) {
          flatTree->pt_muons()      .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).momentum.Pt()  : default_value) );
          flatTree->eta_muons()     .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).momentum.Eta() : default_value) );
          flatTree->phi_muons()     .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).momentum.Phi() : default_value) );
          flatTree->mass_muons()    .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).momentum.M()   : default_value) );
          flatTree->energy_muons()  .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).momentum.E()   : default_value) );
          flatTree->charge_muons()  .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).charge         : default_value) );
          // flatTree->pfRelIso_muons().push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).pfRelIso()     : default_value) ); // FIXME
          // flatTree->passId_muons()  .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).passID()       : default_value) ); // FIXME
          // flatTree->passIso_muons() .push_back( (imuon < filtered_muons_bkg.size() ? filtered_muons_bkg.at(imuon).passIso()      : default_value) ); // FIXME
        }
    }

protected:
    std::shared_ptr<ntuple::FlatTree> flatTree;
    bool writeFlatTree;
    ntuple::TauVector correctedTaus;
    ntuple::MET correctedMET;
    ntuple::MET postRecoilMET;
};
} // analysis
