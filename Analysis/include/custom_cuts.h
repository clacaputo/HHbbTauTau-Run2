/*!
 * \file custom_cuts.h
 * \brief Cuts which are customized for HHbbtautau analysis
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-04-01 created
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

namespace cuts {
    const double minDeltaR_betweenHiggses = 1; // >

    const double DeltaR_MC_Match = 0.3; // <

    namespace skim {
        namespace TauTau {

            namespace tauID {
                const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 10; // <
            }
        }

        namespace ETau {

            const double pFRelIso = 0.5; // <
            namespace tauID {
                const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 10; // <
            }
        }

        namespace MuTau {

            const double pFRelIso = 0.5; // <
            namespace tauID {
                const double againstMuonLoose = 0.5; // >
                const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 10; // <
            }
        }
    }

    namespace btag {
        const double CSVL = 0.244; // loose twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
        const double CSVM = 0.679; //medium twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
        const double CSVT = 0.898; //medium twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5

        namespace signal {
            const double pt = 20; // > AN-2013/075 HH->bb gamma gamma
            const double eta = 2.4; // < twiki HiggsToTauTauWorkingSummer2013#bjets
            const double CSV = CSVM; // recommended twiki HiggsToTauTauWorkingSummer2013#bjets
        }
        namespace veto {
            const double pt = 20; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
            const double eta = 2.4; // < twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
            const double CSV = CSVT; // twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
            const double passLooseID = true; // = twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
            const double deltaR_signalObjects = 0.4; // > twiki HiggsToTauTauWorkingSummer2013#Cuts_for_VH_e_mu_tau
        }
    }

}
