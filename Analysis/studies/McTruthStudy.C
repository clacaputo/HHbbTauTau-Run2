/*!
 * \file McTruthStudy.C
 * \brief Study of base analysis object at level of MC truth.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-24 created
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

#include "Analysis/include/LightBaseBigTreeAnalyzer.h"

class McTruthStudyData : public analysis::BaseAnalyzerData {
public:
    McTruthStudyData(TFile& outputFile) : BaseAnalyzerData(outputFile) {}

    TH1D_ENTRY(DeltaRmin1_visible, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin2_visible, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin1_original, 600, 0, 3)
    TH1D_ENTRY(DeltaRmin2_original, 600, 0, 3)
    TH1D_ENTRY(DeltaRbjets_MC, 600, 0, 6)
    TH1D_ENTRY(deltaRmin_MC, 600, 0, 6)
    TH1D_ENTRY(deltaRmax_visible, 600, 0, 6)
    TH1D_ENTRY(deltaRmax_original, 600, 0, 6)
    TH1D_ENTRY(deltaPtMax, 250, 0, 5)
    TH1D_ENTRY(deltaPtMax_vis, 250, 0, 5)
    TH1D_ENTRY(MinPtBjetsMC, 20, 0, 200)
    TH1D_ENTRY(MassBB_MC, 30, 0, 300)
    TH1D_ENTRY(MassBB_MCvis, 30, 0, 300)
    TH1D_ENTRY(goodElectronsFromZee, 5, -0.5, 4.5)
    TH1D_ENTRY(hardElectronsZee, 5, -0.5, 4.5)
    TH1D_ENTRY(hardMuonsZee, 5, -0.5, 4.5)
    TH1D_ENTRY(goodMuonsFromZmm, 5, -0.5, 4.5)
    TH1D_ENTRY(hardElectronsZmm, 5, -0.5, 4.5)
    TH1D_ENTRY(hardMuonsZmm, 5, -0.5, 4.5)

    TH1D_ENTRY(MC_total_momentum_pt, 20, 0, 200)
    TH1D_ENTRY(MC_total_momentum_eta, 20, 0, 10)
    TH1D_ENTRY(MC_total_momentum_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_total_momentum_energy, 100, 0, 10000)

    TH1D_ENTRY(MC_visible_momentum_pt, 20, 0, 200)
    TH1D_ENTRY(MC_visible_momentum_eta, 20, 0, 10)
    TH1D_ENTRY(MC_visible_momentum_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_visible_momentum_energy, 100, 0, 10000)

    TH1D_ENTRY(MC_invisible_momentum_pt, 20, 0, 200)
    TH1D_ENTRY(MC_invisible_momentum_eta, 20, 0, 10)
    TH1D_ENTRY(MC_invisible_momentum_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_invisible_momentum_energy, 100, 0, 10000)

    TH1D_ENTRY(MC_delta_invisible_momentum_pt, 20, 0, 200)
    TH1D_ENTRY(MC_delta_invisible_momentum_eta, 20, 0, 10)
    TH1D_ENTRY(MC_delta_invisible_momentum_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_delta_invisible_momentum_energy, 100, 0, 10000)

    TH1D_ENTRY(MC_H_invisible_momentum_pt, 20, 0, 200)
    TH1D_ENTRY(MC_H_invisible_momentum_eta, 20, 0, 10)
    TH1D_ENTRY(MC_H_invisible_momentum_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_H_invisible_momentum_energy, 100, 0, 10000)

    TH1D_ENTRY(MC_underlying_invisible_momentum_pt, 20, 0, 200)
    TH1D_ENTRY(MC_underlying_invisible_momentum_eta, 20, 0, 10)
    TH1D_ENTRY(MC_underlying_invisible_momentum_phi, 32, 0, 3.2)
    TH1D_ENTRY(MC_underlying_invisible_momentum_energy, 100, 0, 10000)

    TH1D_ENTRY(MET_predicition_delta_pt, 20, 0, 200)
    TH1D_ENTRY(MET_predicition_delta_phi, 32, 0, 3.2)

    TH1D_ENTRY(MC_invisible_momentum_tau_ratio, 11, 0, 1.1)
};


class McTruthStudy : public analysis::LightBaseBigTreeAnalyzer {
public:
    McTruthStudy(const std::string& channelName, const std::string& inputFileName, const std::string& outputFileName,
                 const std::string& configFileName, const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          LightBaseBigTreeAnalyzer(channelName, inputFileName, outputFileName, configFileName, _prefix,
                                   _maxNumberOfEvents),
          anaData(*outputFile)
    {
    }

protected:
    virtual void AnalyzeSelection(const analysis::SelectionResults& selection) override
    {
        using namespace analysis;

        if(selection.GetFinalStateMC().b_jets.size() < 2 || selection.GetFinalStateMC().taus.size() < 2) return;

        const VisibleGenObject& bjet1_visible = selection.GetFinalStateMC().b_jets.at(0);
        const GenParticle* bjet1_MC = bjet1_visible.origin;

        const VisibleGenObject& bjet2_visible = selection.GetFinalStateMC().b_jets.at(1);
        const GenParticle* bjet2_MC = bjet2_visible.origin;

        const VisibleGenObject tau_1 = selection.GetFinalStateMC().taus.at(0);
        const VisibleGenObject tau_2 = selection.GetFinalStateMC().taus.at(1);
        const GenParticle* tau_MC_1 = tau_1.origin;
        const GenParticle* tau_MC_2 = tau_2.origin;

        TLorentzVector total_momentum, visible_momentum, invisible_momentum;
        for(const GenParticle& genParticle : genEvent.genParticles) {
            if(genParticle.status != particles::FinalStateParticle) continue;

            total_momentum += genParticle.momentum;
            if(particles::neutrinos.count(genParticle.pdg.Code))
                invisible_momentum += genParticle.momentum;
            else
                visible_momentum += genParticle.momentum;
        }
        if(invisible_momentum.Pt() < 40)
            return;
        anaData.MC_total_momentum_pt().Fill(total_momentum.Pt());
        anaData.MC_total_momentum_eta().Fill(std::abs(total_momentum.Eta()));
        anaData.MC_total_momentum_phi().Fill(std::abs(total_momentum.Phi()));
        anaData.MC_total_momentum_energy().Fill(total_momentum.E());

        anaData.MC_visible_momentum_pt().Fill(visible_momentum.Pt());
        anaData.MC_visible_momentum_eta().Fill(std::abs(visible_momentum.Eta()));
        anaData.MC_visible_momentum_phi().Fill(std::abs(visible_momentum.Phi()));
        anaData.MC_visible_momentum_energy().Fill(visible_momentum.E());

        anaData.MC_invisible_momentum_pt().Fill(invisible_momentum.Pt());
        anaData.MC_invisible_momentum_eta().Fill(std::abs(invisible_momentum.Eta()));
        anaData.MC_invisible_momentum_phi().Fill(std::abs(invisible_momentum.Phi()));
        anaData.MC_invisible_momentum_energy().Fill(invisible_momentum.E());

        unsigned n_lept_bjet = 0;
        for(const VisibleGenObject& bjet : selection.GetFinalStateMC().b_jets) {
            if(bjet.finalStateChargedLeptons.size() != 0)
                ++n_lept_bjet;
        }



//            const TLorentzVector MET = MakeLorentzVectorPtEtaPhiM(correctedMET.pt, 0, correctedMET.phi, 0);
        const TLorentzVector MET = MakeLorentzVectorPtEtaPhiM(selection.MET_with_recoil_corrections.pt, 0,
                                                              selection.MET_with_recoil_corrections.phi, 0);
        const TLorentzVector MET_predicition_delta = invisible_momentum - MET;
        anaData.MET_predicition_delta_pt(n_lept_bjet).Fill(MET_predicition_delta.Pt());
        anaData.MET_predicition_delta_phi(n_lept_bjet).Fill(std::abs(MET_predicition_delta.Phi()));

        const VisibleGenObject H(selection.GetFinalStateMC().resonance->mothers.front());
        const GenParticle* genProton = *genEvent.primaryParticles.begin();
        const VisibleGenObject proton(genProton, H.particlesProcessed);
        const TLorentzVector delta_invisible = invisible_momentum - proton.invisibleMomentum - H.invisibleMomentum;

        anaData.MC_delta_invisible_momentum_pt().Fill(delta_invisible.Pt());
        anaData.MC_delta_invisible_momentum_eta().Fill(std::abs(delta_invisible.Eta()));
        anaData.MC_delta_invisible_momentum_phi().Fill(std::abs(delta_invisible.Phi()));
        anaData.MC_delta_invisible_momentum_energy().Fill(delta_invisible.E());

        anaData.MC_H_invisible_momentum_pt().Fill(H.invisibleMomentum.Pt());
        anaData.MC_H_invisible_momentum_eta().Fill(std::abs(H.invisibleMomentum.Eta()));
        anaData.MC_H_invisible_momentum_phi().Fill(std::abs(H.invisibleMomentum.Phi()));
        anaData.MC_H_invisible_momentum_energy().Fill(H.invisibleMomentum.E());

        anaData.MC_underlying_invisible_momentum_pt().Fill(proton.invisibleMomentum.Pt());
        anaData.MC_underlying_invisible_momentum_eta().Fill(std::abs(proton.invisibleMomentum.Eta()));
        anaData.MC_underlying_invisible_momentum_phi().Fill(std::abs(proton.invisibleMomentum.Phi()));
        anaData.MC_underlying_invisible_momentum_energy().Fill(proton.invisibleMomentum.E());

        const VisibleGenObject htt(selection.GetFinalStateMC().Higgs_TauTau->mothers.front());
        const double tau_pt_ratio = htt.invisibleMomentum.Pt() / H.invisibleMomentum.Pt();
        anaData.MC_invisible_momentum_tau_ratio(n_lept_bjet).Fill(tau_pt_ratio);


        if (selection.bjets_all.size() > 1 && bjet1_MC->momentum.Pt() > 20 && bjet2_MC->momentum.Pt() > 20 &&
                std::abs(bjet1_MC->momentum.Eta()) < 2.4 && std::abs(bjet2_MC->momentum.Eta()) < 2.4 ){

            double deltaRmin_firstCouple =
                    std::min(bjet1_MC->momentum.DeltaR(tau_MC_1->momentum),bjet1_MC->momentum.DeltaR(tau_MC_2->momentum));
            double deltaRmin_secondCouple =
                    std::min(bjet2_MC->momentum.DeltaR(tau_MC_1->momentum),bjet2_MC->momentum.DeltaR(tau_MC_2->momentum));
            double deltaRmin_MC = std::min(deltaRmin_firstCouple,deltaRmin_secondCouple);

            //std::cout << "deltaRmin_MC=" << deltaRmin_MC << std::endl;
            anaData.deltaRmin_MC().Fill(deltaRmin_MC);
            anaData.DeltaRbjets_MC().Fill(bjet1_MC->momentum.DeltaR(bjet2_MC->momentum));
            anaData.MinPtBjetsMC().Fill(std::min(bjet1_MC->momentum.Pt(),bjet2_MC->momentum.Pt()));
            TLorentzVector bb_MC = bjet1_MC->momentum + bjet2_MC->momentum;
            TLorentzVector bb_MC_visible = bjet1_visible.visibleMomentum + bjet2_visible.visibleMomentum;

            anaData.MassBB_MC().Fill(bb_MC.M());
            anaData.MassBB_MCvis().Fill(bb_MC_visible.M());

            double deltaRmin1_original = std::numeric_limits<double>::max();
            double deltaRmin2_original = std::numeric_limits<double>::max();
            double deltaRmin1_visible = std::numeric_limits<double>::max();
            double deltaRmin2_visible = std::numeric_limits<double>::max();
            unsigned index_bjet1 = 0;
            unsigned index_bjet2 = 0;
            unsigned index_bjet1_vis = 0;
            unsigned index_bjet2_vis = 0;
            for (unsigned i = 0; i < selection.bjets_all.size(); ++i){
                const Candidate& bjet = selection.bjets_all.at(i);
//                double deltaPt1 = std::abs(bjet.momentum.Pt() - bjet1_MC->momentum.Pt())/bjet1_MC->momentum.Pt();
//                double deltaPt2 = std::abs(bjet.momentum.Pt() - bjet2_MC->momentum.Pt())/bjet2_MC->momentum.Pt();
                if (bjet.momentum.DeltaR(bjet1_MC->momentum) < deltaRmin1_original /*&& deltaPt1 < 0.4*/){
                    deltaRmin1_original = bjet.momentum.DeltaR(bjet1_MC->momentum);
                    index_bjet1 = i;
                }

                if (bjet.momentum.DeltaR(bjet2_MC->momentum) < deltaRmin2_original /*&& deltaPt2 < 0.4*/){
                    deltaRmin2_original = bjet.momentum.DeltaR(bjet2_MC->momentum);
                    index_bjet2 = i;
                }

                if (bjet.momentum.DeltaR(bjet1_visible.visibleMomentum) < deltaRmin1_visible){
                    deltaRmin1_visible = bjet.momentum.DeltaR(bjet1_visible.visibleMomentum);
                    index_bjet1_vis = i;
                }

                if (bjet.momentum.DeltaR(bjet2_visible.visibleMomentum) < deltaRmin2_visible){
                    deltaRmin2_visible = bjet.momentum.DeltaR(bjet2_visible.visibleMomentum);
                    index_bjet2_vis = i;
                }

            }

            anaData.DeltaRmin1_original().Fill(deltaRmin1_original);
            anaData.DeltaRmin2_original().Fill(deltaRmin2_original);

            double deltaRmax_original = std::max(deltaRmin1_original,deltaRmin2_original);
            anaData.deltaRmax_original().Fill(deltaRmax_original);

            anaData.DeltaRmin1_visible().Fill(deltaRmin1_visible);
            anaData.DeltaRmin2_visible().Fill(deltaRmin2_visible);

            double deltaRmax_visible = std::max(deltaRmin1_visible,deltaRmin2_visible);
            anaData.deltaRmax_visible().Fill(deltaRmax_visible);

            double deltaPtMax;
            if (deltaRmax_original == std::numeric_limits<double>::max()){
                deltaPtMax = std::numeric_limits<double>::max();
            }
            else {
                const Candidate& selectedBjets1 = selection.bjets_all.at(index_bjet1);
                const Candidate& selectedBjets2 = selection.bjets_all.at(index_bjet2);
                double deltaPt1 = std::abs(selectedBjets1.momentum.Pt() - bjet1_MC->momentum.Pt());
                double deltaPt2 = std::abs(selectedBjets2.momentum.Pt() - bjet2_MC->momentum.Pt());
                deltaPtMax = std::max(deltaPt1/bjet1_MC->momentum.Pt(),deltaPt2/bjet2_MC->momentum.Pt());
            }
            anaData.deltaPtMax().Fill(deltaPtMax);

            double deltaPtMax_vis;
            if (deltaRmax_visible == std::numeric_limits<double>::max()){
                deltaPtMax_vis = std::numeric_limits<double>::max();
            }
            else {
                const Candidate& selectedBjets1 = selection.bjets_all.at(index_bjet1_vis);
                const Candidate& selectedBjets2 = selection.bjets_all.at(index_bjet2_vis);
                double deltaPt1 = std::abs(selectedBjets1.momentum.Pt() - bjet1_visible.visibleMomentum.Pt());
                double deltaPt2 = std::abs(selectedBjets2.momentum.Pt() - bjet2_visible.visibleMomentum.Pt());
                deltaPtMax_vis =
                        std::max(deltaPt1/bjet1_visible.visibleMomentum.Pt(),deltaPt2/bjet2_visible.visibleMomentum.Pt());
            }
            anaData.deltaPtMax_vis().Fill(deltaPtMax_vis);
        }
    }

private:
    McTruthStudyData anaData;
};
