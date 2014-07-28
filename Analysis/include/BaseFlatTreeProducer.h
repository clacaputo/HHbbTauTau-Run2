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
#include "GenMatch.h" // is it needed?

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
    virtual analysis::Candidate SelectBackgroundElectron(size_t id, cuts::ObjectSelector* objectSelector,
                                                         root_ext::AnalyzerData& _anaData,
                                                         const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::electronVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Electron& object = event->electrons().at(id);
        const analysis::Candidate electron(analysis::Candidate::Electron, id, object);

        cut(true, ">0 ele cand");
        cut(X(pt) > pt, "pt");
        const double eta = std::abs( X(eta) );
        cut(eta < eta_high, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 ele_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(ele_vertex,primaryVertex.position,electron.momentum); // same as dB
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(pfRelIso) < pFRelIso, "pFRelIso");
        const size_t pt_index = object.pt < ref_pt ? 0 : 1;
        const size_t eta_index = eta < scEta_min[0] ? 0 : (eta < scEta_min[1] ? 1 : 2);
        cut(X(mvaPOGNonTrig) > MVApogNonTrig[pt_index][eta_index], "mva");
        cut(X(missingHits) < missingHits, "missingHits");
        cut(X(hasMatchedConversion) == hasMatchedConversion, "conversion");

        return electron;
    }

    virtual analysis::Candidate SelectBackgroundMuon(size_t id, cuts::ObjectSelector* objectSelector,
                                                     root_ext::AnalyzerData& _anaData,
                                                     const std::string& selection_label) override
    {
        using namespace cuts::Htautau_Summer13::muonVeto;
        cuts::Cutter cut(objectSelector);
        const ntuple::Muon& object = event->muons().at(id);
        const analysis::Candidate muon(analysis::Candidate::Mu, id, object);

        cut(true, ">0 mu cand");
        cut(X(pt) > pt, "pt");
        cut(std::abs( X(eta) ) < eta, "eta");
        const double DeltaZ = std::abs(object.vz - primaryVertex.position.Z());
        cut(Y(DeltaZ)  < dz, "dz");
        const TVector3 mu_vertex(object.vx, object.vy, object.vz);
        const double d0_PV = analysis::Calculate_dxy(mu_vertex,primaryVertex.position,muon.momentum);
        cut(std::abs( Y(d0_PV) ) < d0, "d0");
        cut(X(isGlobalMuonPromptTight) == isGlobalMuonPromptTight, "tight");
        cut(X(isPFMuon) == isPFMuon, "PF");
        cut(X(nMatchedStations) > nMatched_Stations, "stations");
        cut(X(pixHits) > pixHits, "pix_hits");
        cut(X(trackerLayersWithMeasurement) > trackerLayersWithMeasurement, "layers");
        cut(X(pfRelIso) < pfRelIso, "pFRelIso");

        return muon;
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

    //from HH_BaseAnalyzer Master branch
    void FindAnalysisBBTauTauFinalState(finalState::bbTauTau& final_state)
    {
        static const particles::ParticleCodes resonanceCodes = { particles::Radion};
        static const particles::ParticleCodes resonanceDecay = { particles::Higgs, particles::Higgs };
        static const particles::ParticleCodes2D HiggsDecays = { { particles::b, particles::b },
                                                                { particles::tau, particles::tau } };
        static const particles::ParticleCodes SM_ResonanceCodes = { particles::Higgs, particles::Z };
        static const particles::ParticleCodes SM_ResonanceDecay_1 = { particles::tau, particles::tau };
        static const particles::ParticleCodes SM_ResonanceDecay_2 = { particles::b, particles::b };


        genEvent.Initialize(event.genParticles());

        const GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);

        if (resonances.size() > 1)
            throw std::runtime_error("more than 1 resonance per event");

        if (resonances.size() == 1){
            final_state.resonance = *resonances.begin();

            GenParticlePtrVector HiggsBosons;
            if(!FindDecayProducts(*final_state.resonance, resonanceDecay,HiggsBosons))
                throw std::runtime_error("Resonance does not decay into 2 Higgs");

            GenParticleVector2D HiggsDecayProducts;
            GenParticleIndexVector HiggsIndexes;
            if(!FindDecayProducts2D(HiggsBosons,HiggsDecays,HiggsDecayProducts,HiggsIndexes))
                throw std::runtime_error("NOT HH -> bb tautau");

            final_state.b_jets = HiggsDecayProducts.at(0);
            final_state.taus = HiggsDecayProducts.at(1);

            final_state.Higgs_TauTau = HiggsBosons.at(HiggsIndexes.at(1));
            final_state.Higgs_BB = HiggsBosons.at(HiggsIndexes.at(0));

        }

        if (resonances.size() == 0){
            //search H->bb H->tautau
            const analysis::GenParticleSet SM_particles = genEvent.GetParticles(SM_ResonanceCodes);

            analysis::GenParticlePtrVector SM_ResonanceToTauTau;
            analysis::GenParticlePtrVector SM_ResonanceToBB;

            for (const GenParticle* SM_particle : SM_particles){
                analysis::GenParticlePtrVector resonanceDecayProducts;
                if(analysis::FindDecayProducts(*SM_particle, SM_ResonanceDecay_1,resonanceDecayProducts)){
                    final_state.taus = resonanceDecayProducts;
                    final_state.resonance = SM_particle;
                    SM_ResonanceToTauTau.push_back(SM_particle);
                }
                else if (analysis::FindDecayProducts(*SM_particle, SM_ResonanceDecay_2,resonanceDecayProducts)){
                    final_state.b_jets = resonanceDecayProducts;
                    final_state.resonance = SM_particle;
                    SM_ResonanceToBB.push_back(SM_particle);
                }
            }
        }


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

    // you might want to move this somewhere else
    bool ComputeAntiElectronMVA3New(int category, float raw, int WP)
    {
      if (category < 0 ) return false ;
      if (category > 15) return true  ;

      float cutsLoose    [16] = {0.835,0.831,0.849,0.859,0.873,0.823,0.85 ,0.855,0.816,0.861,0.862,0.847,0.893,0.82 ,0.845,0.851} ;
      float cutsMedium   [16] = {0.933,0.921,0.944,0.945,0.918,0.941,0.981,0.943,0.956,0.947,0.951,0.95 ,0.897,0.958,0.955,0.942} ;
      float cutsTight    [16] = { 0.96,0.968,0.971,0.972,0.969,0.959,0.981,0.965,0.975,0.972,0.974,0.971,0.897,0.971,0.961,0.97 } ;
      float cutsVeryTight[16] = {0.978,0.98 ,0.982,0.985,0.977,0.974,0.989,0.977,0.986,0.983,0.984,0.983,0.971,0.987,0.977,0.981} ;

      switch (WP)
      {
        case 0 : return (raw > cutsLoose    [category]) ;
        case 1 : return (raw > cutsMedium   [category]) ;
        case 2 : return (raw > cutsTight    [category]) ;
        case 3 : return (raw > cutsVeryTight[category]) ;
      }

      return false ; // to avoid warnings, smarter solutions exist
    }

    void FillFlatTree(const Candidate& higgs      , double m_sv,
                      double m_sv_Up, double m_sv_Down,
                      const CandidateVector& jets , const CandidateVector& jetsPt20,
                      const CandidateVector& bjets, const CandidateVector& retagged_bjets,
                      const VertexVector& vertices, const Candidate& leg1, const Candidate& leg2,
                      const ntuple::MET& pfMET)
    {
        static const float default_value = ntuple::DefaultFloatFillValueForFlatTree();

        // Event
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


        // Weights
        flatTree->puweight()        = PUweight;
        flatTree->trigweight_1()    = triggerWeights.at(0);
        flatTree->trigweight_2()    = triggerWeights.at(1);
        flatTree->idweight_1()      = IDweights.at(0);
        flatTree->idweight_2()      = IDweights.at(1);
        flatTree->isoweight_1()     = IsoWeights.at(0);
        flatTree->isoweight_2()     = IsoWeights.at(1);
        // flatTree->fakeweight()      = fakeWeights.at(1); // what's this?
        flatTree->weight()          = eventWeight;
        flatTree->embeddedWeight()  = 1.; // FIXME! once we have the embedded samples

        // HTT candidate
        flatTree->mvis()           = higgs.momentum.M();
        flatTree->m_sv()           = m_sv;
        flatTree->m_sv_Up()        = m_sv_Up;
        flatTree->m_sv_Down()      = m_sv_Down;
        flatTree->pt_sv()          = default_value; // SVfit not ready yet
        flatTree->eta_sv()         = default_value; // SVfit not ready yet
        flatTree->phi_sv()         = default_value; // SVfit not ready yet
        flatTree->DeltaR_leptons() = leg1.momentum.DeltaR(leg2.momentum) ;
        flatTree->pt_tt()          = (leg1.momentum + leg2.momentum).Pt();

        // Hhh generator info candidate
        mcmatching::MChiggses mch = mcmatching::GetMChiggses(event->genParticles()) ;
        //flatTree->pdgId_resonance_MC() = mch.heavyH.M() ;
        flatTree->pt_resonance_MC()    = mch.heavyH.Pt()  ;
        flatTree->eta_resonance_MC()   = mch.heavyH.Eta() ;
        flatTree->phi_resonance_MC()   = mch.heavyH.Phi() ;
        flatTree->mass_resonance_MC()  = mch.heavyH.M()   ;
        //flatTree->pdgId_Htt_MC()       = mch.heavyH.M() ;
        flatTree->pt_Htt_MC()          = mch.lightH1.Pt()  ;
        flatTree->eta_Htt_MC()         = mch.lightH1.Eta() ;
        flatTree->phi_Htt_MC()         = mch.lightH1.Phi() ;
        flatTree->mass_Htt_MC()        = mch.lightH1.M()   ;
        //flatTree->pdgId_Hbb_MC()       = mch.heavyH.M() ;
        flatTree->pt_Hbb_MC()          = mch.lightH2.Pt()  ;
        flatTree->eta_Hbb_MC()         = mch.lightH2.Eta() ;
        flatTree->phi_Hbb_MC()         = mch.lightH2.Phi() ;
        flatTree->mass_Hbb_MC()        = mch.lightH2.M()   ;
        flatTree->n_extraJets_MC()     = default_value ; // needs to be filles with NUP! https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py#L51

        // Leg 1, lepton
        flatTree->pt_1()     = leg1.momentum.Pt()  ;
        flatTree->phi_1()    = leg1.momentum.Phi() ;
        flatTree->eta_1()    = leg1.momentum.Eta() ;
        flatTree->m_1()      = leg1.momentum.M()   ;
        flatTree->energy_1() = leg1.momentum.E()   ;
        flatTree->q_1()      = leg1.charge         ;
        flatTree->mt_1()     = analysis::Calculate_MT(leg1.momentum, postRecoilMET.pt, postRecoilMET.phi);
        flatTree->d0_1()     = analysis::Calculate_dxy(leg1.vertexPosition, primaryVertex.position,leg1.momentum);
        flatTree->dZ_1()     = leg1.vertexPosition.Z() - primaryVertex.position.Z();

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatchingExercise
        mcmatching::matchedGenParticle genp1 = mcmatching::MatchToGen(leg1,event->genParticles(),1) ; // stable muons
        flatTree->pdgId_1_MC() = genp1.pdgId    ;
        flatTree->pt_1_MC   () = genp1.matched ? genp1.p4.Pt()  : default_value ;
        flatTree->phi_1_MC  () = genp1.matched ? genp1.p4.Phi() : default_value ;
        flatTree->eta_1_MC  () = genp1.matched ? genp1.p4.Eta() : default_value ;
        flatTree->m_1_MC    () = genp1.matched ? genp1.p4.M()   : default_value ;

        // Leg 2, tau
        flatTree->pt_2()     = leg2.momentum.Pt();
        flatTree->phi_2()    = leg2.momentum.Phi();
        flatTree->eta_2()    = leg2.momentum.Eta();
        flatTree->m_2()      = leg2.momentum.M();
        flatTree->energy_2() = leg2.momentum.E();
        flatTree->q_2()      = leg2.charge;
        flatTree->mt_2()     = analysis::Calculate_MT(leg2.momentum, postRecoilMET.pt, postRecoilMET.phi);
        flatTree->d0_2()     = analysis::Calculate_dxy(leg2.vertexPosition, primaryVertex.position,leg2.momentum);
        flatTree->dZ_2()     = leg2.vertexPosition.Z() - primaryVertex.position.Z();

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATMCMatchingExercise
        mcmatching::matchedGenParticle genp2 = mcmatching::MatchToGen(leg2,event->genParticles(),3) ; // taus before fragmentation (not perfect, but acceptable)
        flatTree->pdgId_2_MC() = genp2.pdgId    ;
        flatTree->pt_2_MC   () = genp2.matched ? genp2.p4.Pt()  : default_value ;
        flatTree->phi_2_MC  () = genp2.matched ? genp2.p4.Phi() : default_value ;
        flatTree->eta_2_MC  () = genp2.matched ? genp2.p4.Eta() : default_value ;
        flatTree->m_2_MC    () = genp2.matched ? genp2.p4.M()   : default_value ;

        // RM: for the three channels, mt, et, tt this leg is always a tau
        const ntuple::Tau& ntuple_tau_leg2 = correctedTaus.at(leg2.index);
        flatTree->decayMode_2()                                = ntuple_tau_leg2.decayMode;
        bool antiEloose  = analysis::BaseFlatTreeProducer::ComputeAntiElectronMVA3New(ntuple_tau_leg2.againstElectronMVA3category,
                                                                                      ntuple_tau_leg2.againstElectronMVA3raw,
                                                                                      0) ;
        bool antiEMedium = analysis::BaseFlatTreeProducer::ComputeAntiElectronMVA3New(ntuple_tau_leg2.againstElectronMVA3category,
                                                                                      ntuple_tau_leg2.againstElectronMVA3raw,
                                                                                      1) ;
        bool antiETight  = analysis::BaseFlatTreeProducer::ComputeAntiElectronMVA3New(ntuple_tau_leg2.againstElectronMVA3category,
                                                                                      ntuple_tau_leg2.againstElectronMVA3raw,
                                                                                      2) ;
        bool antiEVTight = analysis::BaseFlatTreeProducer::ComputeAntiElectronMVA3New(ntuple_tau_leg2.againstElectronMVA3category,
                                                                                      ntuple_tau_leg2.againstElectronMVA3raw,
                                                                                      3) ;
        flatTree->againstElectronLooseMVA_2()                  = antiEloose  ;
        flatTree->againstElectronMediumMVA_2()                 = antiEMedium ;
        flatTree->againstElectronTightMVA_2()                  = antiETight  ;
        flatTree->againstElectronVTightMVA_2()                 = antiEVTight ;
        flatTree->againstElectronLoose_2()                     = ntuple_tau_leg2.againstElectronLoose  ;
        flatTree->againstElectronMedium_2()                    = ntuple_tau_leg2.againstElectronMedium ;
        flatTree->againstElectronTight_2()                     = ntuple_tau_leg2.againstElectronTight  ;
        flatTree->againstMuonLoose_2()                         = ntuple_tau_leg2.againstMuonLoose      ;
        flatTree->againstMuonMedium_2()                        = ntuple_tau_leg2.againstMuonMedium     ;
        flatTree->againstMuonTight_2()                         = ntuple_tau_leg2.againstMuonTight      ;
        flatTree->byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = ntuple_tau_leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits ;

        // MET
        TLorentzVector postRecoilMetMomentum;
        postRecoilMetMomentum.SetPtEtaPhiM(postRecoilMET.pt, 0, postRecoilMET.phi, 0.);
        flatTree->pt_tt_MET() = (leg1.momentum + leg2.momentum + postRecoilMetMomentum).Pt();

        flatTree->met()       = pfMET.pt          ;
        flatTree->metphi()    = pfMET.phi         ;
        flatTree->mvamet()    = postRecoilMET.pt  ;
        flatTree->mvametphi() = postRecoilMET.phi ;
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

        // Jets
        flatTree->njets()     = jets.size()           ;
        flatTree->njetspt20() = jetsPt20.size()       ;
        flatTree->nBjets()    = retagged_bjets.size() ;

        // RM: SAVE ALL THE JETS WITH PT > 20 (AND SAVE THE CSV DISCRIMINATOR)
        std::vector<ntuple::Jet> csv_sorted_jetsPt20 ;
        for (unsigned int ijet = 0 ; ijet < jetsPt20.size() ; ++ijet) {
          ntuple::Jet myJet = event->jets().at(jetsPt20.at(ijet).index) ;
          csv_sorted_jetsPt20.push_back(myJet) ;
        }

        // sort the jets by pt
        std::sort( csv_sorted_jetsPt20.begin(), csv_sorted_jetsPt20.end(), []( ntuple::Jet a, ntuple::Jet b ){ return a.pt > b.pt; } ) ;

        // fill the jet collections
        for (unsigned int ijet = 0 ; ijet < csv_sorted_jetsPt20.size() ; ++ijet) {
          flatTree->pt_jets()      .push_back( csv_sorted_jetsPt20.at(ijet).pt  );
          flatTree->eta_jets()     .push_back( csv_sorted_jetsPt20.at(ijet).eta );
          flatTree->phi_jets()     .push_back( csv_sorted_jetsPt20.at(ijet).phi );
          flatTree->mass_jets()    .push_back( csv_sorted_jetsPt20.at(ijet).mass);
          // flatTree->energy_jets()  .push_back( csv_sorted_jetsPt20.at(ijet).E()); // not in ntuple::Jet
          flatTree->csv_jets()     .push_back( csv_sorted_jetsPt20.at(ijet).combinedSecondaryVertexBJetTags );
          // inspect the flavour of the gen jet
          mcmatching::genjetflavour genjetinfo = mcmatching::IsGenBJet(csv_sorted_jetsPt20.at(ijet),event->genParticles()) ;
          flatTree->isJet_MC_Bjet()                  .push_back( genjetinfo.isBjet                                 );
          flatTree->isJet_MC_Bjet_withLeptonicDecay().push_back( genjetinfo.leptonsInJet.size() > 0 ? true : false );
        }

        // sort the jets by csv discriminator
        std::sort( csv_sorted_jetsPt20.begin(), csv_sorted_jetsPt20.end(), []( ntuple::Jet a, ntuple::Jet b ){ return a.combinedSecondaryVertexBJetTags > b.combinedSecondaryVertexBJetTags; } ) ;

        // fill the bjet collections
        for (unsigned int ibjet = 0 ; ibjet < csv_sorted_jetsPt20.size() ; ++ibjet) {
          flatTree->pt_Bjets()      .push_back( csv_sorted_jetsPt20.at(ibjet).pt  );
          flatTree->eta_Bjets()     .push_back( csv_sorted_jetsPt20.at(ibjet).eta );
          flatTree->phi_Bjets()     .push_back( csv_sorted_jetsPt20.at(ibjet).phi );
          flatTree->mass_Bjets()    .push_back( csv_sorted_jetsPt20.at(ibjet).mass);
          // flatTree->energy_Bjets()  .push_back( csv_sorted_jetsPt20.at(ibjet).E()); // not in ntuple::Jet
          flatTree->csv_Bjets()     .push_back( csv_sorted_jetsPt20.at(ibjet).combinedSecondaryVertexBJetTags );
          // inspect the flavour of the gen jet
          mcmatching::genjetflavour genjetinfo = mcmatching::IsGenBJet(csv_sorted_jetsPt20.at(ibjet),event->genParticles()) ;
          flatTree->isBjet_MC_Bjet()                  .push_back( genjetinfo.isBjet                                 );
          flatTree->isBjet_MC_Bjet_withLeptonicDecay().push_back( genjetinfo.leptonsInJet.size() > 0 ? true : false );
        }

        // PV
        ntuple::Vertex PV = event->vertices().front() ;
        flatTree->x_PV() = PV.x ;
        flatTree->y_PV() = PV.y ;
        flatTree->z_PV() = PV.z ;

    }

protected:
    std::shared_ptr<ntuple::FlatTree> flatTree;
    bool writeFlatTree;
    ntuple::TauVector correctedTaus;
    ntuple::MET correctedMET;
    ntuple::MET postRecoilMET;
};
} // analysis
