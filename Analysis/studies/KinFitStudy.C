/*!
 * \file KinFitStudy.C
 * \brief Study of kinematic fit performance.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \author Maria Agnese Ciocci (Siena University, INFN Pisa)
 * \date 2014-10-15 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>,
 *                Maria Agnese Ciocci <mariaagnese.ciocci@pi.infn.it>
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

//#include "Analysis/include/KinFit_CMSSW.h"
#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"

class KinFitStudyData : public root_ext::AnalyzerData {
public:
    KinFitStudyData(TFile& outputFile) : AnalyzerData(outputFile) {}

    TH1D_ENTRY(HHKinFit_convergence, 10, -3.5, 6.5)
    TH1D_ENTRY(HHKinFit_M_bbtt, 30, 0, 600)
    TH1D_ENTRY(HHKinFit_chi2, 10, 0, 50)
    TH1D_ENTRY(HHKinFit_fitprob, 20, 0, 1)
    TH1D_ENTRY(HHKinFit_pull_balance, 20, 0, 10)
    TH1D_ENTRY(HHKinFit_pull_balance_1, 20, 0, 10)
    TH1D_ENTRY(HHKinFit_pull_balance_2, 20, 0, 10)
};

class KinFitStudy : public analysis::LightBaseFlatTreeAnalyzer {
public:
    KinFitStudy(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName), anaData(GetOutputFile())
    {
        recalc_kinfit = true;
        anaData.getOutputFile().cd();
    }

protected:
    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;
        using namespace analysis::kinematic_fit;

        if(DetermineEventRegion(eventInfo,category) != analysis::EventRegion::OS_Isolated) return;
        if (!analysis::TwoJetsEventCategories_MediumBjets.count(category)) return;

        const four_body::FitResults& four_body_result_HHKinFit = eventInfo.fitResults;
        anaData.HHKinFit_convergence(category).Fill(four_body_result_HHKinFit.convergence);
        if(four_body_result_HHKinFit.convergence > 0) {
            anaData.HHKinFit_M_bbtt(category).Fill(four_body_result_HHKinFit.mass);
            anaData.HHKinFit_chi2(category).Fill(four_body_result_HHKinFit.chi2);
            anaData.HHKinFit_fitprob(category).Fill(four_body_result_HHKinFit.fit_probability);
            anaData.HHKinFit_pull_balance(category).Fill(four_body_result_HHKinFit.pull_balance);
            anaData.HHKinFit_pull_balance_1(category).Fill(four_body_result_HHKinFit.pull_balance_1);
            anaData.HHKinFit_pull_balance_2(category).Fill(four_body_result_HHKinFit.pull_balance_2);
        }

//        const two_body::FitInput two_body_input(eventInfo.bjet_momentums.at(0), eventInfo.bjet_momentums.at(1));
//        const two_body::FitResults two_body_result_KinFitter = two_body::Fit_KinFitter(two_body_input);
//        anaData.KinFitter_convergence(category).Fill(two_body_result_KinFitter.convergence);
//        if(!two_body_result_KinFitter.convergence) {
//            const TLorentzVector bbtt = two_body_result_KinFitter.bjet_momentums.at(0)
//                    + two_body_result_KinFitter.bjet_momentums.at(1)
//                    + eventInfo.lepton_momentums.at(0)
//                    + eventInfo.lepton_momentums.at(1);
//            anaData.KinFitter_M_bbtt(category).Fill(bbtt.M());
//        }
    }

private:
    KinFitStudyData anaData;
};

