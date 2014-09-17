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
double CorrectMassByKinfit(const TLorentzVector& bjet1, const TLorentzVector& bjet2, const TLorentzVector& tau1,
                           const TLorentzVector& tau2, const TLorentzVector& MVAmet, const TMatrixD& metcov)
{
    static const bool debug = true;

    std::vector<Int_t> hypo_mh1;
    hypo_mh1.push_back(125);
    std::vector<Int_t> hypo_mh2;
    hypo_mh2.push_back(125);

//    std::cout << "metcov00 = " << metcov[0][0] << ", metcov01 = " << metcov[0][1] << ",metcov10 = " << metcov[1][0]
//            << ", metcov11 = " << metcov[1][1] << std::endl;
//    std::cout << "metDet = " << metcov.Determinant() << std::endl;
    //intance of fitter master class
    HHKinFitMaster kinFits = HHKinFitMaster(&bjet1,&bjet2,&tau1,&tau2);
    kinFits.setAdvancedBalance(&MVAmet,metcov);
    //kinFits.setSimpleBalance(ptmiss.Pt(),10); //alternative which uses only the absolute value of ptmiss in the fit
    kinFits.addMh1Hypothesis(hypo_mh1);
    kinFits.addMh2Hypothesis(hypo_mh2);
    kinFits.doFullFit();

    Double_t chi2_best = kinFits.getBestChi2FullFit();
    Double_t mh_best = kinFits.getBestMHFullFit();
    std::pair<Int_t, Int_t> bestHypo = kinFits.getBestHypoFullFit();
    std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_chi2 = kinFits.getChi2FullFit();
//    std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_fitprob = kinFits.getFitProbFullFit();
    std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_mH = kinFits.getMHFullFit();
//    std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b1 = kinFits.getPullB1FullFit();
//    std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_b2 = kinFits.getPullB2FullFit();
    std::map< std::pair<Int_t, Int_t>, Double_t> fit_results_pull_balance = kinFits.getPullBalanceFullFit();
    std::map< std::pair<Int_t, Int_t>, Int_t> fit_convergence = kinFits.getConvergenceFullFit();



    std::pair< Int_t, Int_t > hypo(hypo_mh1.at(0),hypo_mh2.at(0));

    if (debug){
        std::cout << "best chi2:       " << chi2_best << std::endl;
        std::cout << "best hypothesis: " << bestHypo.first << " " << bestHypo.second << std::endl;
        std::cout << "best mH=         " << mh_best << std::endl;
        std::cout << "fit_convergence=         " << fit_convergence.at(hypo) << std::endl;
        std::cout << "fit_results_pull_balance=         " << fit_results_pull_balance.at(hypo) << std::endl;
        std::cout << "fit_results_mH=         " << fit_results_mH.at(hypo) << std::endl;
    }

    if (fit_convergence.at(hypo)>0 && fit_results_chi2.at(hypo)<25 && fit_results_pull_balance.at(hypo)>0)
        return fit_results_mH.at(hypo);
    //std::cerr << "mass with kin Fit cannot be calculated" << std::endl;
    return -100;

}

} // analysis
