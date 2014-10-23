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

#include "Analysis/include/KinFit_CMSSW.h"
#include "Analysis/include/LightBaseFlatTreeAnalyzer.h"

class KinFitStudyData : public analysis::LightFlatAnalyzerData {
public:
    KinFitStudyData(TFile& outputFile) : LightFlatAnalyzerData(outputFile) {}

    TH1D_ENTRY(HHKinFit_M_bbtt, 30, 0, 600)
    TH1D_ENTRY(HHKinFit_chi2, 10, 0, 50)
    TH1D_ENTRY(HHKinFit_convergence, 10, 0, 10)
    TH1D_ENTRY(HHKinFit_Failed, 100, -1000, 1000)
    TH1D_ENTRY(HHKinFit_Not_Failed, 100, -1000, 1000)
    TH1D_ENTRY(HHKinFit_Failed_cov, 100, -1000, 1000)
    TH1D_ENTRY(HHKinFit_pulb_99_chi2, 100, -1000, 1000)
    TH1D_ENTRY(HHKinFit_pulb_99_convergence, 10, 0, 10)
    TH1D_ENTRY(KinFitter_M_bbtt, 30, 0, 600)
    TH1D_ENTRY(KinFitter_convergence, 10, 0, 10)
    TH1D_ENTRY(Delta_ptHtt_2bjet_met, 100, 0, 200)
    TH1D_ENTRY(Delta_ptHSVFIT_2bjet_met, 100, 0, 200)   
    TH1D_ENTRY(Delta_ptHSVFIT_met, 100, 0, 200)


};

class KinFitStudy : public analysis::LightBaseFlatTreeAnalyzer {
public:
    KinFitStudy(const std::string& inputFileName, const std::string& outputFileName)
         : LightBaseFlatTreeAnalyzer(inputFileName, outputFileName), anaData(GetOutputFile())
    {
        anaData.getOutputFile().cd();
    }

protected:
    virtual analysis::LightFlatAnalyzerData& GetAnaData() override { return anaData; }

    virtual void AnalyzeEvent(const analysis::FlatEventInfo& eventInfo, analysis::EventCategory category) override
    {
        using analysis::EventCategory;
        using namespace analysis::kinematic_fit;

        if(!PassSelection(eventInfo, category)) return;
        if (category != EventCategory::TwoJets_ZeroBtag && category != EventCategory::TwoJets_OneBtag
                && category != EventCategory::TwoJets_TwoBtag) return;
//mac here check on visible+missing pt 
        TLorentzVector Htt_MET_plus_bjets = eventInfo.bjet_momentums.at(0) + eventInfo.bjet_momentums.at(1) +
                                            eventInfo.Htt_MET;
/*        std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$" << std::endl;
        std::cout << "bjet 1: " << eventInfo.bjet_momentums.at(0) << std::endl;
        std::cout << "bjet 2: " << eventInfo.bjet_momentums.at(1) << std::endl;
        std::cout << "Htt_MET: " << eventInfo.Htt_MET << std::endl;
        std::cout << "Htt_MET_plus_bjets: " << Htt_MET_plus_bjets << std::endl;
*/        
        anaData.Delta_ptHtt_2bjet_met(category).Fill(Htt_MET_plus_bjets.Pt());

        TLorentzVector Htt_SVFIT;
        Htt_SVFIT.SetPtEtaPhiM(eventInfo.event->pt_sv_MC, eventInfo.Htt.Eta(), eventInfo.Htt.Phi(), eventInfo.event->m_sv_MC);
        TLorentzVector MET;
        MET.SetPtEtaPhiM(eventInfo.event->met,0,eventInfo.event->metphi,0);
        TLorentzVector Htt_SVFIT_plus_bjets = eventInfo.Htt + MET;
	for(size_t n = 0; n < eventInfo.bjet_momentums.size(); ++n )
        	Htt_SVFIT_plus_bjets += eventInfo.bjet_momentums.at(n);
         TLorentzVector Htt_SVFIT_plus_MET= Htt_SVFIT +  
                                              eventInfo.MET;
/*        std::cout << "bjet 1: " << eventInfo.bjet_momentums.at(0) << std::endl;
        std::cout << "bjet 2: " << eventInfo.bjet_momentums.at(1) << std::endl;

        std::cout << "Htt_SVFIT: " << Htt_SVFIT << std::endl;
        std::cout << "pt_sv_MC: " << eventInfo.event->pt_sv_MC << std::endl;
        std::cout << "Htt_SVFIT_plus_bjets: " << Htt_SVFIT_plus_bjets << std::endl;
*/
        anaData.Delta_ptHSVFIT_2bjet_met(category).Fill(Htt_SVFIT_plus_bjets.Pt());
         anaData.Delta_ptHSVFIT_met(category).Fill(Htt_SVFIT_plus_MET.Pt());

// mac here end check
        const four_body::FitInput four_body_input(eventInfo.bjet_momentums.at(0), eventInfo.bjet_momentums.at(1),
                                                  eventInfo.lepton_momentums.at(0), eventInfo.lepton_momentums.at(1),
                                                  eventInfo.MET, eventInfo.MET_covariance);
        const four_body::FitResults four_body_result_HHKinFit = four_body::Fit(four_body_input);
        anaData.HHKinFit_convergence(category).Fill(four_body_result_HHKinFit.convergence);
        anaData.HHKinFit_chi2(category).Fill(four_body_result_HHKinFit.chi2);
//failed & not failed fit
// # failed (quality cuts +cov or x2 negative

        if(!four_body_result_HHKinFit.has_valid_mass )anaData.HHKinFit_Failed(category).Fill(999);
//#     of valid
        if(four_body_result_HHKinFit.has_valid_mass )anaData.HHKinFit_Not_Failed(category).Fill(999);

//# of filed because cov or x2 
        if(four_body_result_HHKinFit.pull_balance == -99 ) anaData.HHKinFit_Failed_cov(category).Fill(999);       
//         std::cerr << "four body mass with kin Fit cannot be calculated" << std::endl;
	if(four_body_result_HHKinFit.pull_balance == -99 ){
       		anaData.HHKinFit_pulb_99_convergence(category).Fill(four_body_result_HHKinFit.convergence);
		if (std::isnan(four_body_result_HHKinFit.chi2)){
			anaData.HHKinFit_pulb_99_chi2(category).Fill(999);		
		}
		else anaData.HHKinFit_pulb_99_chi2(category).Fill(four_body_result_HHKinFit.chi2);
	}
        if(four_body_result_HHKinFit.convergence > 0 && four_body_result_HHKinFit.chi2 < 25)
            anaData.HHKinFit_M_bbtt(category).Fill(four_body_result_HHKinFit.mass);

        const two_body::FitInput two_body_input(eventInfo.bjet_momentums.at(0), eventInfo.bjet_momentums.at(1));
        const two_body::FitResults two_body_result_KinFitter = two_body::Fit_KinFitter(two_body_input);
        anaData.KinFitter_convergence(category).Fill(two_body_result_KinFitter.convergence);
        if(!two_body_result_KinFitter.convergence) {
            const TLorentzVector bbtt = two_body_result_KinFitter.bjet_momentums.at(0)
                    + two_body_result_KinFitter.bjet_momentums.at(1)
                    + eventInfo.lepton_momentums.at(0)
                    + eventInfo.lepton_momentums.at(1);
            anaData.KinFitter_M_bbtt(category).Fill(bbtt.M());
        }
    }

private:
    KinFitStudyData anaData;
};

