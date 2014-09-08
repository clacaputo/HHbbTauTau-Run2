/*!
 * \file FlatTreeAnalyzer_etau.C
 * \brief Analyze flat-tree for etau channel for Htautau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \author Rosamaria Venditti (INFN Bari, Bari University)
 * \author Claudio Caputo (INFN Bari, Bari University)
 * \date 2014-09-03 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
 *                Maria Teresa Grippo <grippomariateresa@gmail.com>,
 *                Rosamaria Venditti,
 *                Claudio Caputo
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

#include "../include/BaseFlatTreeAnalyzer.h"
#include "../include/Htautau_Summer13.h"

class FlatTreeAnalyzer_etau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_etau(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                          const std::string& outputFileName, const std::string& _signalName,
                          const std::string& _dataName, bool _WjetsData = false, bool _isBlind=false)
        : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, _signalName, _dataName, _WjetsData,
                               _isBlind)
    {
    }

protected:
    virtual analysis::EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) override
    {
        using analysis::EventType_QCD;

        if (event.againstElectronLooseMVA_2 > 0.5){
            // OS - isolated
            if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < 1.5 &&
                    event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < 1.5 &&
                    (event.q_1 * event.q_2 == -1))
                return EventType_QCD::OS_Isolated;
            else if (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < 1.5 &&
                    event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < 1.5 &&
                    (event.q_1 * event.q_2 == +1)) // SS - isolated
                return EventType_QCD::SS_Isolated;
            else if ((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= 1.5 ||
                     event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= 1.5) &&
                     (event.q_1 * event.q_2 == -1)) // OS - not isolated
                return EventType_QCD::OS_NotIsolated;
            else if ((event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= 1.5 ||
                     event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= 1.5) &&
                     (event.q_1 * event.q_2 == +1)) // SS - not isolated
                return EventType_QCD::SS_NotIsolated;
            else
                return EventType_QCD::Unknown;
        }
        else
            return EventType_QCD::Unknown;
    }

    virtual analysis::EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) override
    {
        using analysis::EventType_Wjets;
        return EventType_Wjets::Unknown;
    }

    virtual void EstimateWjets() override
    {

    }

};
