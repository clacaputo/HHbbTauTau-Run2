/*!
 * \file JetEnergyUncertainty.h
 * \brief Definition of wrapper to estimate jet energy uncertainties.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-11-16 created
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

#include <TLorentzVector.h>

#include "FWCore/Utilities/src/EDMException.cc"
#include "FWCore/Utilities/src/Exception.cc"
#include "FWCore/Utilities/src/typelookup.cc"

#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc"

#include "TreeProduction/interface/Jet.h"
#include "AnalysisBase/include/AnalysisMath.h"

namespace analysis {

class JetEnergyUncertaintyCorrector {
public:
    JetEnergyUncertaintyCorrector(const std::string& input_file_name, const std::string& section_name)
        : parameters(input_file_name, section_name), jetCorrector(parameters) {}

    void ApplyCorrection(ntuple::Jet& jet, bool scale_up)
    {
        try {
            jetCorrector.setJetPt(jet.pt);
            jetCorrector.setJetEta(jet.eta);
            const double uncertainty = jetCorrector.getUncertainty(scale_up);
            const double sign = scale_up ? +1 : -1;
            const double sf = 1.0 + sign * uncertainty;
            const TLorentzVector original_momentum = MakeLorentzVectorPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass);
            const TLorentzVector corrected_momentum = original_momentum * sf;
            jet.pt = corrected_momentum.Pt();
            jet.eta = corrected_momentum.Eta();
            jet.phi = corrected_momentum.Phi();
            jet.mass = corrected_momentum.M();
        } catch(cms::Exception&) {
            std::cerr << "WARNING: Jet uncertainty is not calculated for jet with pt = " << jet.pt
                      << ", eta = " << jet.eta << std::endl;
        }
    }

    void ApplyCorrection(ntuple::JetVector& jets, bool scale_up)
    {
        for(ntuple::Jet& jet : jets)
            ApplyCorrection(jet, scale_up);
    }

private:
    void CmsswWarningsWorkaround() const
    {
        float fx[0], fy[0];
        quadraticInterpolation(0, fx, fy);
    }

private:
    JetCorrectorParameters parameters;
    JetCorrectionUncertainty jetCorrector;
};

} // namespace analysis
