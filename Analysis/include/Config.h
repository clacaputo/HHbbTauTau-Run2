/*!
 * \file Config.h
 * \brief Definition of Config class.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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

#include "AnalysisBase/include/BaseConfig.h"

namespace analysis {

class Config : public BaseConfig {
public:
    Config(const std::string& fileName) { Read(fileName); }

    ANA_CONFIG_PARAMETER(unsigned, ReportInterval, 10)
    ANA_CONFIG_PARAMETER(bool, RunSingleEvent, false)
    ANA_CONFIG_PARAMETER(unsigned, SingleEventId, 0)
    ANA_CONFIG_PARAMETER(unsigned, MaxTreeVersion, 1)

    ANA_CONFIG_PARAMETER(bool, isMC, false)
    ANA_CONFIG_PARAMETER(bool, ApplyTauESCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrectionForW, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrectionForZ, false)
    ANA_CONFIG_PARAMETER(bool, ExpectedOneNonSMResonance, false)
    ANA_CONFIG_PARAMETER(bool, ExpectedAtLeastOneSMResonanceToTauTauOrToBB, false)
    ANA_CONFIG_PARAMETER(bool, RequireSpecificFinalState, false)
    ANA_CONFIG_PARAMETER(bool, DoZEventCategorization, false)
    ANA_CONFIG_PARAMETER(bool, isDYEmbeddedSample, false)
    ANA_CONFIG_PARAMETER(bool, isTTEmbeddedSample, false)
    ANA_CONFIG_PARAMETER(bool, ApplyEtoTauFakeRate, false)
    ANA_CONFIG_PARAMETER(bool, ApplyJetToTauFakeRate, false)
    ANA_CONFIG_PARAMETER(bool, ApplyDMweight, false)

    ANA_CONFIG_PARAMETER(bool, EstimateTauEnergyUncertainties, false)
    ANA_CONFIG_PARAMETER(bool, EstimateJetEnergyUncertainties, false)
    ANA_CONFIG_PARAMETER(bool, EstimateBtagEfficiencyUncertainties, false)
    ANA_CONFIG_PARAMETER(bool, EstimateBtagFakeUncertainties, false)
    ANA_CONFIG_PARAMETER(std::string, JetEnergyUncertainties_inputFile, "")
    ANA_CONFIG_PARAMETER(std::string, JetEnergyUncertainties_inputSection, "")

    ANA_CONFIG_PARAMETER(bool, ApplyPUreweight, false)
    ANA_CONFIG_PARAMETER(std::string, PUreweight_fileName, "")
    ANA_CONFIG_PARAMETER(double, PUreweight_maxAvailablePU, 60.0)
    ANA_CONFIG_PARAMETER(double, PUreweight_defaultWeight, 0.0)

    ANA_CONFIG_PARAMETER(double, MvaMet_dZcut, 0.1)
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameU, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameDPhi, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameCovU1, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameCovU2, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo_MuTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData_MuTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC_MuTau, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo_ETau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData_ETau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC_ETau, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo_TauTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData_TauTau, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC_TauTau, "")

    bool extractMCtruth()
    {
        return ApplyTauESCorrection() || ApplyRecoilCorrection() || RequireSpecificFinalState()
                || ExpectedOneNonSMResonance() || ExpectedAtLeastOneSMResonanceToTauTauOrToBB()
                || DoZEventCategorization() || ApplyEtoTauFakeRate() || IsEmbeddedSample() || ApplyDMweight();
    }

    bool IsEmbeddedSample() { return isDYEmbeddedSample() || isTTEmbeddedSample(); }

};

} // analysis
