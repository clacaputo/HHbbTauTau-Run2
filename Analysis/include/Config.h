/*!
 * \file Config.h
 * \brief Definition of Config class.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>
 *
 * This file is part of X->HH->bbTauTau.
 *
 * X->HH->bbTauTau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * VoltageSourceControl is distributed in the hope that it will be useful,
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

    ANA_CONFIG_PARAMETER(std::string, ReweightFileName, "none")
    ANA_CONFIG_PARAMETER(bool, UseMCtruth, false)
    ANA_CONFIG_PARAMETER(bool, ApplyTauESCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ApplyPostRecoilCorrection, false)
    ANA_CONFIG_PARAMETER(bool, ExpectedOneResonanceToTauTau, false)
    ANA_CONFIG_PARAMETER(bool, RequireSpecificFinalState, false)

};

} // analysis
