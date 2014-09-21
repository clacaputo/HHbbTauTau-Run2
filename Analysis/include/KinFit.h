/*!
 * \file SVfit.h
 * \brief Definition of wrappers for SVfit_standalone code.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \author Francesco Brivio (Milano Bicocca University, INFN Milano)
 * \date 2014-09-17 created
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

#include "HHKinFit/source/HHEventRecord.cpp"
#include "HHKinFit/source/HHKinFit.cpp"
#include "HHKinFit/source/HHKinFitMaster.cpp"
#include "HHKinFit/source/HHParticle.cpp"
#include "HHKinFit/source/HHParticleList.cpp"
#include "HHKinFit/source/HHV4Vector.cpp"
#include "HHKinFit/source/PSMath.cpp"
#include "HHKinFit/source/PSTools.cpp"

#include "TreeProduction/interface/MET.h"

#include "AnalysisBase/include/Candidate.h"

namespace analysis {

namespace kinematic_fit {

struct FourBodyFitInput {
    std::vector<TLorentzVector> bjet_momentums;
    std::vector<TLorentzVector> tau_momentums;
    TLorentzVector mvaMET;
    TMatrixD metCov;

    FourBodyFitInput(const TLorentzVector& bjet1, const TLorentzVector& bjet2,
                     const TLorentzVector& tau1, const TLorentzVector& tau2,
                     const TLorentzVector& _mvaMET, const TMatrixD& _metCov)
        : mvaMET(_mvaMET), metCov(_metCov)
    {
        bjet_momentums.push_back(bjet1);
        bjet_momentums.push_back(bjet2);
        tau_momentums.push_back(tau1);
        tau_momentums.push_back(tau2);
    }

    FourBodyFitInput(const FourBodyFitInput& input, double bjet_scale_factor, double tau_scale_factor)
        : mvaMET(input.mvaMET), metCov(input.metCov)
    {
        for(const TLorentzVector& momentum : input.bjet_momentums)
            bjet_momentums.push_back(momentum * bjet_scale_factor);

        for(const TLorentzVector& momentum : input.tau_momentums)
            tau_momentums.push_back(momentum * tau_scale_factor);
    }
};

struct FitResults {
    static constexpr double default_value = std::numeric_limits<double>::lowest();

    bool has_valid_mass;
    double mass;

    FitResults() : has_valid_mass(false), mass(default_value) {}
};

struct FitResultsWithUncertainties {
    FitResults fit_bb_tt;
    FitResults fit_bb_down_tt_down;
    FitResults fit_bb_down_tt_up;
    FitResults fit_bb_up_tt_down;
    FitResults fit_bb_up_tt_up;

    FitResults fit_bb;
    FitResults fit_bb_down;
    FitResults fit_bb_up;
};

inline FitResults FourBodyFit(const FourBodyFitInput& input)
{
    static const bool debug = false;
    static const Int_t higgs_mass_hypotesis = 125;
    static const Int_t convergence_cut = 0;
    static const Double_t chi2_cut = 25;
    static const Double_t pull_balance_cut = 0;

    FitResults result;

    const std::vector<Int_t> hypo_mh1 = { higgs_mass_hypotesis };
    const std::vector<Int_t> hypo_mh2 = { higgs_mass_hypotesis };

    if(debug) {
        input.metCov.Print();
        std::cout << "metDet = " << input.metCov.Determinant() << std::endl;
    }

    //intance of fitter master class
    HHKinFitMaster kinFits = HHKinFitMaster(&input.bjet_momentums.at(0), &input.bjet_momentums.at(1),
                                            &input.tau_momentums.at(0), &input.tau_momentums.at(1));
    kinFits.setAdvancedBalance(&input.mvaMET, input.metCov);
    //kinFits.setSimpleBalance(ptmiss.Pt(),10); //alternative which uses only the absolute value of ptmiss in the fit
    kinFits.addMh1Hypothesis(hypo_mh1);
    kinFits.addMh2Hypothesis(hypo_mh2);
    kinFits.doFullFit();

    if(debug) {
        const Double_t chi2_best = kinFits.getBestChi2FullFit();
        const Double_t mh_best = kinFits.getBestMHFullFit();
        const std::pair<Int_t, Int_t> bestHypo = kinFits.getBestHypoFullFit();

        std::cout << "best chi2:       " << chi2_best << std::endl;
        std::cout << "best hypothesis: " << bestHypo.first << " " << bestHypo.second << std::endl;
        std::cout << "best mH=         " << mh_best << std::endl;
    }

    const std::pair<Int_t, Int_t> hypo(hypo_mh1.at(0), hypo_mh2.at(0));
    const Int_t convergence = kinFits.getConvergenceFullFit().at(hypo);
    const Double_t chi2 = kinFits.getChi2FullFit().at(hypo);
    const Double_t pull_balance = kinFits.getPullBalanceFullFit().at(hypo);
    const Double_t mH = kinFits.getMHFullFit().at(hypo);

    if (debug) {
        const Double_t fitprob = kinFits.getFitProbFullFit().at(hypo);
        const Double_t pull_b1 = kinFits.getPullB1FullFit().at(hypo);
        const Double_t pull_b2 = kinFits.getPullB2FullFit().at(hypo);

        std::cout << "fit convergence =  " << convergence << std::endl;
        std::cout << "fit chi2 =         " << chi2 << std::endl;
        std::cout << "fit fitprob =      " << fitprob << std::endl;
        std::cout << "fit pull_b1 =      " << pull_b1 << std::endl;
        std::cout << "fit pull_b2 =      " << pull_b2 << std::endl;
        std::cout << "fit pull_balance = " << pull_balance << std::endl;
        std::cout << "fit mH =           " << mH << std::endl;
    }

    if (convergence > convergence_cut && chi2 < chi2_cut && pull_balance > pull_balance_cut) {
        result.mass = mH;
        result.has_valid_mass = true;
    }
    else if(debug)
        std::cerr << "four body mass with kin Fit cannot be calculated" << std::endl;

    return result;
}

inline FitResultsWithUncertainties FitWithUncertainties(const FourBodyFitInput& input, double bjet_energy_uncertainty,
                                                        double tau_energy_uncertainty, bool fit_two_bjets,
                                                        bool fit_four_body)
{
    FitResultsWithUncertainties result;
    if(fit_four_body) {
        result.fit_bb_tt = FourBodyFit(input);
        result.fit_bb_down_tt_down = FourBodyFit(FourBodyFitInput(input, 1 - bjet_energy_uncertainty,
                                                                  1 - tau_energy_uncertainty));
        result.fit_bb_down_tt_up = FourBodyFit(FourBodyFitInput(input, 1 - bjet_energy_uncertainty,
                                                                1 + tau_energy_uncertainty));
        result.fit_bb_up_tt_down = FourBodyFit(FourBodyFitInput(input, 1 + bjet_energy_uncertainty,
                                                                1 - tau_energy_uncertainty));
        result.fit_bb_up_tt_up = FourBodyFit(FourBodyFitInput(input, 1 + bjet_energy_uncertainty,
                                                              1 + tau_energy_uncertainty));
    }
    return result;
}

} // namespace kinematic_fit
} // namespace analysis
