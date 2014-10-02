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

#include "HHKinFit/src/HHDiJetKinFit.cpp"
#include "HHKinFit/src/HHDiJetKinFitMaster.cpp"
#include "HHKinFit/src/HHEventRecord.cpp"
#include "HHKinFit/src/HHKinFit.cpp"
#include "HHKinFit/src/HHKinFitMaster.cpp"
#include "HHKinFit/src/HHParticle.cpp"
#include "HHKinFit/src/HHParticleList.cpp"
#include "HHKinFit/src/HHTauInputKinFit.cpp"
#include "HHKinFit/src/HHTauInputKinFitMaster.cpp"
#include "HHKinFit/src/HHV4Vector.cpp"
#include "HHKinFit/src/PSMath.cpp"
#include "HHKinFit/src/PSTools.cpp"




#include "FWCore/Utilities/src/EDMException.cc"
#include "FWCore/Utilities/src/Exception.cc"

#include "FWCore/MessageLogger/src/MessageLogger.cc"
#include "FWCore/MessageLogger/src/MessageLoggerQ.cc"
#include "FWCore/MessageLogger/src/MessageDrop.cc"
#include "FWCore/MessageLogger/src/ELseverityLevel.cc"
#include "FWCore/MessageLogger/src/MessageSender.cc"
#include "FWCore/MessageLogger/src/AbstractMLscribe.cc"
#include "FWCore/MessageLogger/src/ErrorObj.cc"
#include "FWCore/MessageLogger/src/ELextendedID.cc"
#include "FWCore/MessageLogger/src/ELstring.cc"

using namespace std;
#include "PhysicsTools/KinFitter/src/TAbsFitConstraint.cc"
#include "PhysicsTools/KinFitter/src/TAbsFitParticle.cc"
#include "PhysicsTools/KinFitter/src/TFitConstraintEp.cc"
#include "PhysicsTools/KinFitter/src/TFitConstraintM.cc"
#include "PhysicsTools/KinFitter/src/TFitConstraintMGaus.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleCart.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleECart.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEScaledMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleESpher.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEtEtaPhi.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleEtThetaPhi.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCCart.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCPInvSpher.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMCSpher.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleMomDev.cc"
#include "PhysicsTools/KinFitter/src/TFitParticleSpher.cc"
#include "PhysicsTools/KinFitter/src/TKinFitter.cc"
#include "PhysicsTools/KinFitter/src/TSLToyGen.cc"

#include "TreeProduction/interface/MET.h"

#include "AnalysisBase/include/Candidate.h"

namespace analysis {

namespace kinematic_fit {

namespace four_body {
struct FitInput {
    std::vector<TLorentzVector> bjet_momentums;
    std::vector<TLorentzVector> tau_momentums;
    TLorentzVector mvaMET;
    TMatrixD metCov;

    FitInput(const TLorentzVector& bjet1, const TLorentzVector& bjet2,
                     const TLorentzVector& tau1, const TLorentzVector& tau2,
                     const TLorentzVector& _mvaMET, const TMatrixD& _metCov)
        : mvaMET(_mvaMET), metCov(_metCov)
    {
        bjet_momentums.push_back(bjet1);
        bjet_momentums.push_back(bjet2);
        tau_momentums.push_back(tau1);
        tau_momentums.push_back(tau2);
    }

    FitInput(const FitInput& input, double bjet_scale_factor, double tau_scale_factor)
        : mvaMET(input.mvaMET), metCov(input.metCov)
    {
        for(const TLorentzVector& momentum : input.bjet_momentums)
            bjet_momentums.push_back(momentum * bjet_scale_factor);

        for(const TLorentzVector& momentum : input.tau_momentums)
            tau_momentums.push_back(momentum * tau_scale_factor);
    }
};

struct FitResults {
    static constexpr double default_double_value = std::numeric_limits<double>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    bool has_valid_mass;
    double mass;
    int convergence;
    double chi2;
    double pull_balance;

    FitResults()
        : has_valid_mass(false), mass(default_double_value), convergence(default_int_value), chi2(default_double_value),
          pull_balance(default_double_value)
    {}
};

inline FitResults Fit(const FitInput& input)
{
    static const bool debug = false;
    static const Int_t higgs_mass_hypotesis = 125;
    static const Int_t convergence_cut = 0;
    static const Double_t chi2_cut = 25;
    static const Double_t pull_balance_cut = 0;

    const std::vector<Int_t> hypo_mh1 = { higgs_mass_hypotesis };
    const std::vector<Int_t> hypo_mh2 = { higgs_mass_hypotesis };

    if(debug) {
        std::cout << "Format: (Pt, eta, phi, E)\n";
        std::cout << "b1 momentum: " << input.bjet_momentums.at(0) << std::endl;
        std::cout << "b2 momentum: " << input.bjet_momentums.at(1) << std::endl;
        std::cout << "tau1 momentum: " << input.tau_momentums.at(0) << std::endl;
        std::cout << "tau2 momentum: " << input.tau_momentums.at(1) << std::endl;
        std::cout << "MET: " << input.mvaMET << std::endl;
        std::cout << "MET covariance:\n";
        input.metCov.Print();
        std::cout << "metDet = " << input.metCov.Determinant() << std::endl;
    }

    //intance of fitter master class
    TLorentzVector b1(input.bjet_momentums.at(0)), b2(input.bjet_momentums.at(1)),
            tau1(input.tau_momentums.at(0)), tau2(input.tau_momentums.at(1)),
            mvaMET(input.mvaMET);
    HHKinFitMaster kinFit(&b1, &b2, &tau1, &tau2);
    kinFit.setAdvancedBalance(&mvaMET, input.metCov);
    //kinFits.setSimpleBalance(ptmiss.Pt(),10); //alternative which uses only the absolute value of ptmiss in the fit
    kinFit.addMh1Hypothesis(hypo_mh1);
    kinFit.addMh2Hypothesis(hypo_mh2);
    kinFit.doFullFit();

    if(debug) {
        const Double_t chi2_best = kinFit.getBestChi2FullFit();
        const Double_t mh_best = kinFit.getBestMHFullFit();
        const std::pair<Int_t, Int_t> bestHypo = kinFit.getBestHypoFullFit();

        std::cout << "best chi2:       " << chi2_best << std::endl;
        std::cout << "best hypothesis: " << bestHypo.first << " " << bestHypo.second << std::endl;
        std::cout << "best mH=         " << mh_best << std::endl;
    }

    const std::pair<Int_t, Int_t> hypo(hypo_mh1.at(0), hypo_mh2.at(0));
    const Int_t convergence = kinFit.getConvergenceFullFit().at(hypo);
    const Double_t chi2 = kinFit.getChi2FullFit().at(hypo);
    const Double_t pull_balance = kinFit.getPullBalanceFullFit().at(hypo);
    const Double_t mH = kinFit.getMHFullFit().at(hypo);

    if (debug) {
        const Double_t fitprob = kinFit.getFitProbFullFit().at(hypo);
        const Double_t pull_b1 = kinFit.getPullB1FullFit().at(hypo);
        const Double_t pull_b2 = kinFit.getPullB2FullFit().at(hypo);

        std::cout << "fit convergence =  " << convergence << std::endl;
        std::cout << "fit chi2 =         " << chi2 << std::endl;
        std::cout << "fit fitprob =      " << fitprob << std::endl;
        std::cout << "fit pull_b1 =      " << pull_b1 << std::endl;
        std::cout << "fit pull_b2 =      " << pull_b2 << std::endl;
        std::cout << "fit pull_balance = " << pull_balance << std::endl;
        std::cout << "fit mH =           " << mH << std::endl;
    }

    FitResults result;
    result.convergence = convergence;
    result.chi2 = chi2;
    result.mass = mH;
    result.pull_balance = pull_balance;
    result.has_valid_mass = convergence > convergence_cut && chi2 < chi2_cut && pull_balance > pull_balance_cut;
    if(!result.has_valid_mass && debug)
        std::cerr << "four body mass with kin Fit cannot be calculated" << std::endl;

    return result;
}

} // namespace four_body

namespace two_body {

struct FitInput {
    std::vector<TLorentzVector> bjet_momentums;

    FitInput(const TLorentzVector& bjet1, const TLorentzVector& bjet2)
    {
        bjet_momentums.push_back(bjet1);
        bjet_momentums.push_back(bjet2);
    }

    FitInput(const FitInput& input, double bjet_scale_factor)
    {
        for(const TLorentzVector& momentum : input.bjet_momentums)
            bjet_momentums.push_back(momentum * bjet_scale_factor);
    }
};

struct FitResults {
    static constexpr double default_double_value = std::numeric_limits<double>::lowest();
    static constexpr int default_int_value = std::numeric_limits<int>::lowest();

    int convergence;
    double chi2;
    std::vector<TLorentzVector> bjet_momentums;

    FitResults() : convergence(default_int_value), chi2(default_double_value) {}
};


FitResults Fit(const FitInput& input)
{
    static const bool debug = false;
    static const Int_t higgs_mass_hypotesis = 125;

    const std::vector<Int_t> hypo_mh = { higgs_mass_hypotesis };

    if(debug) {
        std::cout << "Format: (Pt, eta, phi, E)\n";
        std::cout << "b1 momentum: " << input.bjet_momentums.at(0) << std::endl;
        std::cout << "b2 momentum: " << input.bjet_momentums.at(1) << std::endl;
    }

    //intance of fitter master class
    TLorentzVector b1(input.bjet_momentums.at(0)), b2(input.bjet_momentums.at(1));
    HHDiJetKinFitMaster kinFit(&b1, &b2);
    kinFit.addMhHypothesis(hypo_mh);
    kinFit.doFullFit();

    if(debug) {
        const Double_t chi2_best = kinFit.getBestChi2FullFit();
        const std::pair<Int_t, Int_t> bestHypo = kinFit.getBestHypoFullFit();

        std::cout << "best chi2:       " << chi2_best << std::endl;
        std::cout << "best hypothesis: " << bestHypo.first << " " << bestHypo.second << std::endl;
    }

    const std::pair<Int_t, Int_t> hypo(hypo_mh.at(0), -1);
    const Int_t convergence = kinFit.getConvergenceFullFit().at(hypo);
    const Double_t chi2 = kinFit.getBestChi2FullFit();

    if (debug) {
        const Double_t fitprob = kinFit.getFitProbFullFit().at(hypo);
        const Double_t pull_b1 = kinFit.getPullB1FullFit().at(hypo);
        const Double_t pull_b2 = kinFit.getPullB2FullFit().at(hypo);

        std::cout << "fit convergence =  " << convergence << std::endl;
        std::cout << "fit chi2 =         " << chi2 << std::endl;
        std::cout << "fit fitprob =      " << fitprob << std::endl;
        std::cout << "fit pull_b1 =      " << pull_b1 << std::endl;
        std::cout << "fit pull_b2 =      " << pull_b2 << std::endl;
    }

    FitResults result;
    result.convergence = convergence;
    result.chi2 = chi2;
    result.bjet_momentums.push_back(kinFit.getFitJet1());
    result.bjet_momentums.push_back(kinFit.getFitJet2());

    return result;
}

} // namespace two_body

struct FitResultsWithUncertainties {
    four_body::FitResults fit_bb_tt;
    four_body::FitResults fit_bb_down_tt_down;
    four_body::FitResults fit_bb_down_tt_up;
    four_body::FitResults fit_bb_up_tt_down;
    four_body::FitResults fit_bb_up_tt_up;

    two_body::FitResults fit_bb;
    two_body::FitResults fit_bb_down;
    two_body::FitResults fit_bb_up;
};

/*
inline Double_t ErrEta(Float_t Et, Float_t Eta)
{
    Double_t InvPerr2, a, b, c;
    if(fabs(Eta) < 1.4){
        a = 1.215;
        b = 0.037;
        c = 7.941 * 0.0001;
    } else {
        a = 1.773;
        b = 0.034;
        c = 3.56 * 0.0001;
    }
    InvPerr2 = a/(Et * Et) + b/Et + c;
    return InvPerr2;
}

inline Double_t ErrEt(Float_t Et, Float_t Eta)
{
    Double_t InvPerr2, a, b, c;
      if(fabs(Eta) < 1.4){
        a = 5.6;
    b = 1.25;
        c = 0.033;
      }else{
        a = 4.8;
        b = 0.89;
        c = 0.043;
    }
    InvPerr2 = (a * a) + (b * b) * Et + (c * c) * Et * Et;
    return InvPerr2;
}

inline Double_t ErrPhi(Float_t Et, Float_t Eta)
{
Double_t InvPerr2, a, b, c;
  if(fabs(Eta) < 1.4){
    a = 6.65;
    b = 0.04;
    c = 8.49 * 0.00001;
  }else{
      a = 2.908;
      b = 0.021;
      c = 2.59 * 0.0001;
  }
InvPerr2 = a/(Et * Et) + b/Et + c;
  return InvPerr2;
}


inline FitResults TwoBodyFit(const TLorentzVector& momentum_1, const TLorentzVector& momentum_2)
{
    FitResults result;
    TMatrixD m1(3, 3);
    TMatrixD m2(3, 3);
    m1.Zero();
    m2.Zero();

    m1(0, 0) = ErrEt(momentum_1.Et(), momentum_1.Eta());
    m1(1, 1) = ErrEta(momentum_1.Et(), momentum_1.Eta());
    m1(2, 2) = ErrPhi(momentum_1.Et(), momentum_1.Eta());

    m2(0, 0) = ErrEt(momentum_2.Et(), momentum_2.Eta());
    m2(1, 1) = ErrEta(momentum_2.Et(), momentum_2.Eta());
    m2(2, 2) = ErrPhi(momentum_2.Et(), momentum_2.Eta());

    TLorentzVector mom1(momentum_1), mom2(momentum_2);
    TFitParticleEtEtaPhi _jet1(&mom1, &m1);
    TFitParticleEtEtaPhi _jet2(&mom2, &m2);

    TFitConstraintM m_bb;
    m_bb.addParticle1(&_jet1);
    m_bb.addParticle1(&_jet2);
    m_bb.setMassConstraint(125.0);

    TKinFitter _fitter;
    _fitter.addMeasParticle(&_jet1);
    _fitter.addMeasParticle(&_jet2);
    _fitter.addConstraint(&m_bb);

    _fitter.setMaxNbIter(30);
    _fitter.setMaxDeltaS(1e-2);
    _fitter.setMaxF(1e-1);
    _fitter.setVerbosity(0);
    const Int_t fit_result = _fitter.fit();

    std::cout << "fit result: " << fit_result << std::endl;
    if(fit_result == 0) {
        std::cout << "Old 1: " << momentum_1 << std::endl;
        std::cout << "Old 2: " << momentum_2 << std::endl;
        std::cout << "Old mass 1: " << momentum_1.M() << std::endl;
        std::cout << "Old mass 2: " << momentum_2.M() << std::endl;
        std::cout << "New 1: " << *_jet1.getCurr4Vec() << std::endl;
        std::cout << "New 2: " << *_jet2.getCurr4Vec() << std::endl;
        std::cout << "New mass 1: " << _jet1.getCurr4Vec()->M() << std::endl;
        std::cout << "New mass 2: " << _jet2.getCurr4Vec()->M() << std::endl;
        const double m_old = (momentum_1 + momentum_2).M();
        const double m_new = ((*_jet1.getCurr4Vec()) + (*_jet2.getCurr4Vec())).M();
        std::cout << "M old: " << m_old << std::endl;
        std::cout << "M new: " << m_new << std::endl;
    }
    return result;
}
*/
inline FitResultsWithUncertainties FitWithUncertainties(const four_body::FitInput& input, double bjet_energy_uncertainty,
                                                        double tau_energy_uncertainty, bool fit_two_bjets,
                                                        bool fit_four_body)
{
    FitResultsWithUncertainties result;
    if(fit_four_body) {
        result.fit_bb_tt = four_body::Fit(input);
        result.fit_bb_down_tt_down = four_body::Fit(four_body::FitInput(input, 1 - bjet_energy_uncertainty,
                                                                        1 - tau_energy_uncertainty));
        result.fit_bb_down_tt_up = four_body::Fit(four_body::FitInput(input, 1 - bjet_energy_uncertainty,
                                                                      1 + tau_energy_uncertainty));
        result.fit_bb_up_tt_down = four_body::Fit(four_body::FitInput(input, 1 + bjet_energy_uncertainty,
                                                                      1 - tau_energy_uncertainty));
        result.fit_bb_up_tt_up = four_body::Fit(four_body::FitInput(input, 1 + bjet_energy_uncertainty,
                                                                    1 + tau_energy_uncertainty));
    }
//    if(fit_two_bjets) {
//        const two_body::FitInput two_body_input(input.bjet_momentums.at(0), input.bjet_momentums.at(1));
//        result.fit_bb = two_body::Fit(two_body_input);
//        result.fit_bb_down = two_body::Fit(two_body::FitInput(two_body_input, 1 - bjet_energy_uncertainty));
//        result.fit_bb_up = two_body::Fit(two_body::FitInput(two_body_input, 1 + bjet_energy_uncertainty));
//    }
    return result;
}

typedef std::pair<size_t, size_t> FitId;
typedef std::map<FitId, FitResultsWithUncertainties> FitResultsMap;

} // namespace kinematic_fit
} // namespace analysis
