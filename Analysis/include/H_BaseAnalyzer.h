/*!
 * \file H_BaseAnalyzer.h
 * \brief Definition of H_BaseAnalyzer class which is the base class for all H->tautau analyzers.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-07 created
 */

#pragma once

#include "BaseAnalyzer.h"
#include "SyncTree.h"

namespace analysis {

class H_BaseAnalyzer : public BaseAnalyzer {
public:
    H_BaseAnalyzer(const std::string& inputFileName, const std::string& outputFileName,
                 const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0, bool _useMCtruth = false,
                 const std::string& reweightFileName = "none")
        : BaseAnalyzer(inputFileName, outputFileName, _prefix, _maxNumberOfEvents, _useMCtruth, reweightFileName),
          syncTree("syncTree")
    {}

    virtual ~H_BaseAnalyzer()
    {
        GetAnaData().getOutputFile().cd();
        syncTree.Write();
    }

protected:
    void ApplyTauCorrections(const finalState::TauTau& mcFinalState)
    {
        static const double deltaR = 0.3;
        if(!useMCtruth) {
            correctedTaus = event.taus();
            correctedMET = event.metMVA();
            return;
        }
        correctedTaus.clear();

        TLorentzVector sumCorrectedTaus, sumTaus;
        for(const ntuple::Tau& tau : event.taus()) {
            TLorentzVector momentum;
            momentum.SetPtEtaPhiE(tau.pt, tau.eta, tau.phi, tau.energy);
            sumTaus += momentum;
            double scaleFactor = 1.0;

            const bool hasMCmatch = analysis::FindMatchedParticles(momentum, mcFinalState.taus, deltaR).size() != 0;
            if(hasMCmatch) {
                if(tau.decayMode == ntuple::tau_id::kOneProng1PiZero)
                    scaleFactor = 1.025 + 0.001 * std::min(std::max(momentum.Pt() - 45.0, 0.0), 10.0);
                else if(tau.decayMode == ntuple::tau_id::kThreeProng0PiZero)
                    scaleFactor = 1.012 + 0.001 * std::min(std::max(momentum.Pt() - 32.0, 0.0), 18.0);
            }
            const TLorentzVector correctedMomentum = momentum * scaleFactor;
            ntuple::Tau correctedTau(tau);
            correctedTau.pt = correctedMomentum.Pt();
            correctedTau.eta = correctedMomentum.Eta();
            correctedTau.phi = correctedMomentum.Phi();
            correctedTau.energy = correctedMomentum.E();
            correctedTaus.push_back(correctedTau);
            sumCorrectedTaus += correctedMomentum;
        }
        TLorentzVector met, metCorrected;
        met.SetPtEtaPhiM(event.metMVA().pt, event.metMVA().eta, event.metMVA().phi, 0.);
        metCorrected = met + sumTaus - sumCorrectedTaus;
        correctedMET = event.metMVA();
        correctedMET.pt = metCorrected.Pt();
        correctedMET.eta = metCorrected.Eta();
        correctedMET.phi = metCorrected.Phi();
    }

    CandidateVector ApplyCorrections(const CandidateVector& higgses, const GenParticle* resonance)
    {
        CandidateVector candidates;
        for (const Candidate& higgs : higgses){
            ntuple::MET postRecoilMET(correctedMET);
            if (useMCtruth){
                postRecoilMET = ApplyPostRecoilCorrection(correctedMET, higgs.momentum, resonance->momentum);
            }
            const Candidate corrected_h_tautau = CorrectMassBySVfit(higgs, postRecoilMET);
            candidates.push_back(corrected_h_tautau);
        }
        return candidates;
    }

    const analysis::Candidate& SelectSemiLeptonicHiggs(const analysis::CandidateVector& higgses)
    {
        if(!higgses.size())
            throw std::runtime_error("no available higgs candidate to select");
        const auto higgsSelector = [&] (const analysis::Candidate& first, const analysis::Candidate& second) -> bool
        {
            const double first_Pt1 = first.daughters.at(0).momentum.Pt();
            const double first_Pt2 = first.daughters.at(1).momentum.Pt();
            const double first_sumPt = first_Pt1 + first_Pt2;
            const double second_Pt1 = second.daughters.at(0).momentum.Pt();
            const double second_Pt2 = second.daughters.at(1).momentum.Pt();
            const double second_sumPt = second_Pt1 + second_Pt2;

            return first_sumPt > second_sumPt;
        };
        return *std::max_element(higgses.begin(), higgses.end(), higgsSelector);
    }


    bool FindAnalysisFinalState(analysis::finalState::TauTau& final_state)
    {
        static const analysis::ParticleCodes resonanceCodes = { particles::Higgs, particles::Z };
        static const analysis::ParticleCodes resonanceDecay = { particles::tau, particles::tau };

        genEvent.Initialize(event.genParticles());

        const analysis::GenParticleSet resonances = genEvent.GetParticles(resonanceCodes);
        if (resonances.size() != 1)
            throw std::runtime_error("not one resonance per event");

        final_state.resonance = *resonances.begin();

        analysis::GenParticlePtrVector resonanceDecayProducts;
        if(!analysis::FindDecayProducts(*final_state.resonance, resonanceDecay,resonanceDecayProducts)) {
            std::cout << "event id = " << event.eventId().eventId << std::endl;
//            genEvent.PrintChain(final_state.resonance);
//            throw std::runtime_error("Resonance does not decay into 2 taus");
            std::cerr << "Resonance does not decay into 2 taus\n";
            return false;
        }

        final_state.taus = resonanceDecayProducts;

        return true;
    }


    virtual void FillSyncTree(const analysis::Candidate& higgs)
    {
        syncTree.run() = event.eventInfo().run;
        syncTree.lumi() = event.eventInfo().lumis;
        syncTree.evt() = event.eventInfo().EventId;
        for (unsigned n = 0; n < event.eventInfo().bunchCrossing.size(); ++n){
            if (event.eventInfo().bunchCrossing.at(n) == 0){
                syncTree.npu() = event.eventInfo().nPU.at(n); //only in-time PU
            }
        }
        syncTree.puweight() = weight;
        Double_t DMweight = 1;
        analysis::Candidate cand1 = higgs.finalStateDaughters.at(0);
        analysis::Candidate cand2 = higgs.finalStateDaughters.at(1);
        if (cand1.momentum.Pt() < cand2.momentum.Pt()){
            cand1 = higgs.finalStateDaughters.at(1);
            cand2 = higgs.finalStateDaughters.at(0);
        }
        const ntuple::Tau& leg1 = correctedTaus.at(cand1.index);
        const ntuple::Tau& leg2 = correctedTaus.at(cand2.index);
        if (leg1.decayMode == 0)
            DMweight *= 0.88;
        if (leg2.decayMode == 0)
            DMweight *= 0.88;
        syncTree.decaymodeweight() = DMweight;
        syncTree.mvis() = higgs.momentum.M();

        syncTree.pt_1() = cand1.momentum.Pt();
        syncTree.eta_1() = cand1.momentum.Eta();
        syncTree.phi_1() = cand1.momentum.Phi();
        syncTree.m_1() = cand1.momentum.M();
        syncTree.q_1() = leg1.charge;
        syncTree.iso_1() = leg1.byIsolationMVAraw;
        syncTree.passid_1() = 1;
        syncTree.passiso_1() = 1;
        syncTree.mt_1() = cand1.momentum.Mt();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_1() = leg1.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_1() = leg1.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_1() = leg1.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_1() = leg1.byIsolationMVA2raw;
        syncTree.againstMuonLoose_1() = leg1.againstMuonLoose;
        syncTree.againstMuonLoose2_1() = leg1.againstMuonLoose2;
        syncTree.againstMuonMedium2_1() = leg1.againstMuonMedium2;
        syncTree.againstMuonTight2_1() = leg1.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_1() = leg1.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_1() = leg1.againstElectronLoose;

        syncTree.pt_2() = cand2.momentum.Pt();
        syncTree.eta_2() = cand2.momentum.Eta();
        syncTree.phi_2() = cand2.momentum.Phi();
        syncTree.m_2() = cand2.momentum.M();
        syncTree.q_2() = leg2.charge;
        syncTree.iso_2() = leg2.byIsolationMVAraw;
        syncTree.passid_2() = 1;
        syncTree.passiso_2() = 1;
        syncTree.mt_2() = cand2.momentum.Mt();
        syncTree.byCombinedIsolationDeltaBetaCorrRaw3Hits_2() = leg2.byCombinedIsolationDeltaBetaCorrRaw3Hits;
        syncTree.againstElectronMVA3raw_2() = leg2.againstElectronMVA3raw;
        syncTree.againstElectronMVA3category_2() = leg2.againstElectronMVA3category;
        syncTree.byIsolationMVA2raw_2() = leg2.byIsolationMVA2raw;
        syncTree.againstMuonLoose_2() = leg2.againstMuonLoose;
        syncTree.againstMuonLoose2_2() = leg2.againstMuonLoose2;
        syncTree.againstMuonMedium2_2() = leg2.againstMuonMedium2;
        syncTree.againstMuonTight2_2() = leg2.againstMuonTight2;
        syncTree.againstElectronLooseMVA3_2() = leg2.againstElectronLooseMVA3;
        syncTree.againstElectronLoose_2() = leg2.againstElectronLoose;

        syncTree.met() = event.metPF().pt_uncorrected; //raw
        syncTree.metphi() = event.metPF().phi_uncorrected; //raw
        syncTree.mvamet() = event.metMVA().pt;
        syncTree.mvametphi() = event.metMVA().phi;
        syncTree.metcov00() = event.metPF().significanceMatrix.at(0);
        syncTree.metcov01() = event.metPF().significanceMatrix.at(1);
        syncTree.metcov10() = event.metPF().significanceMatrix.at(2);
        syncTree.metcov11() = event.metPF().significanceMatrix.at(3);
        syncTree.mvacov00() = event.metMVA().significanceMatrix.at(0);
        syncTree.mvacov01() = event.metMVA().significanceMatrix.at(1);
        syncTree.mvacov10() = event.metMVA().significanceMatrix.at(2);
        syncTree.mvacov11() = event.metMVA().significanceMatrix.at(3);

        syncTree.Fill();
    }

protected:
    ntuple::SyncTree syncTree;
    ntuple::TauVector correctedTaus;
    ntuple::MET correctedMET;
};
} // analysis
