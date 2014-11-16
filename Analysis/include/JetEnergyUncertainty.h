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

#include "CondFormats/JetMETObjects/src/Utilities.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectorParameters.cc"
#include "CondFormats/JetMETObjects/src/SimpleJetCorrectionUncertainty.cc"
#include "CondFormats/JetMETObjects/src/JetCorrectionUncertainty.cc"

#include "TreeProduction/interface/Jet.h"

namespace analysis {

class JetEnergyUncertaintyCorrector {
public:
    JetEnergyUncertaintyCorrector(const std::string& input_file_name)
        : parameters(input_file_name, "Total"), jetCorrector(parameters)
    {

    }

    ntuple::Jet ApplyCorrection(const ntuple::Jet& original_jet, bool scale_up)
    {
        ntuple::Jet corrected_jet = original_jet;
        jetCorrector.setJetPt(original_jet.pt);
        jetCorrector.setJetEta(original_jet.eta);
        const double uncer
    }

private:
    JetCorrectorParameters parameters;
    JetCorrectionUncertainty jetCorrector;
};

} // namespace analysis
