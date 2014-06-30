/*!
 * \file Config.h
 * \brief Definition of Config class.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
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

#include "BaseConfig.h"

namespace analysis {

class Config : public BaseConfig {
public:
    Config(const std::string& fileName) { Read(fileName); }

    ANA_CONFIG_PARAMETER(unsigned, ReportInterval, 10)
    ANA_CONFIG_PARAMETER(bool, RunSingleEvent, false)
    ANA_CONFIG_PARAMETER(unsigned, SingleEventId, 0)

    ANA_CONFIG_PARAMETER(bool, isMC, false)
    ANA_CONFIG_PARAMETER(bool, ApplyTauESCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ApplyRecoilCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ExpectedOneResonanceToTauTau, false)
    ANA_CONFIG_PARAMETER(bool, RequireSpecificFinalState, false)

    ANA_CONFIG_PARAMETER(bool, ApplyPUreweight, false)
    ANA_CONFIG_PARAMETER(std::string, PUreweight_fileName, "")
    ANA_CONFIG_PARAMETER(double, PUreweight_maxAvailablePU, 60.0)
    ANA_CONFIG_PARAMETER(double, PUreweight_defaultWeight, 0.0)

    ANA_CONFIG_PARAMETER(double, MvaMet_dZcut, 0.1)
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameU, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameDPhi, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameCovU1, "")
    ANA_CONFIG_PARAMETER(std::string, MvaMet_inputFileNameCovU2, "")

    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileCorrectTo, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmData, "")
    ANA_CONFIG_PARAMETER(std::string, RecoilCorrection_fileZmmMC, "")

    bool extractMCtruth()
    {
        return ApplyTauESCorrection() || ApplyRecoilCorrection() ||
                RequireSpecificFinalState() || ExpectedOneResonanceToTauTau();
    }

};

} // analysis