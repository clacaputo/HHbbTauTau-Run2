/*!
 * \file McTruthStudy_flat.C
 * \brief Study of base analysis object at level of MC truth.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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

#include "Analysis/include/KinFit_CMSSW.h"
#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"

class McTruthStudyData : public analysis::LightFlatAnalyzerData {
public:
    McTruthStudyData(std::shared_ptr<TFile> outputFile) : LightFlatAnalyzerData(outputFile) {}

    TH1D_ENTRY(Delta_Pt, 40, -200, 200)
    TH1D_ENTRY(Pt_Delta, 20, 0, 200)
    TH1D_ENTRY(Delta_Phi, 32, -3.2, 3.2)
};

class McTruthStudy_flat : public analysis::LightBaseFlatTreeAnalyzer {
public:
    McTruthStudy_flat(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName), anaData(GetOutputFile()) {}

protected:
    virtual analysis::LightFlatAnalyzerData& GetAnaData() override { return anaData; }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;

        if(!PassSelection(eventInfo, category)) return;
        if (!analysis::TwoJetsEventCategories.count(category)) return;

        TLorentzVector Htt_sv;
        Htt_sv.SetPtEtaPhiM(eventInfo.event->pt_sv_MC, eventInfo.Htt.Eta(), eventInfo.Htt.Phi(), eventInfo.event->m_sv_MC);

        unsigned n_leptonic_decays = 0;
        for(size_t n = 0; n < 2; ++n) {
            if(!eventInfo.event->isBjet_MC_Bjet.at(n))
                return;
            if(eventInfo.event->isBjet_MC_Bjet.at(n) && eventInfo.event->isBjet_MC_Bjet_withLeptonicDecay.at(n))
                ++n_leptonic_decays;
        }

        std::ostringstream ss_name;
        ss_name << category << "_n_lept_b_" << n_leptonic_decays;

        const TLorentzVector H_true = analysis::MakeLorentzVectorPtEtaPhiM(eventInfo.event->pt_resonance_MC,
                                                                           eventInfo.event->eta_resonance_MC,
                                                                           eventInfo.event->phi_resonance_MC,
                                                                           eventInfo.event->mass_resonance_MC);

        const TLorentzVector sum = eventInfo.Htt + eventInfo.Hbb + eventInfo.MET;
        const TLorentzVector delta = H_true - sum;
        anaData.Delta_Pt(ss_name.str()).Fill(H_true.Pt() - sum.Pt());
        anaData.Delta_Phi(ss_name.str()).Fill(delta.Phi());
        anaData.Pt_Delta(ss_name.str()).Fill(delta.Pt());
    }

private:
    McTruthStudyData anaData;
};
