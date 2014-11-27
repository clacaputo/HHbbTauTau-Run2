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
const double DeltaR_MC_Match = 0.3; // <

namespace massWindow{
    const double m_tautau_low = 90;
    const double m_tautau_high = 150;
    const double m_bb_low = 70;
    const double m_bb_high = 150;
}

namespace WjetsBackgroundEstimation {
    const double HighMtRegion = 70; // > For W-jets data driven estimation
    const double HighMtRegion_low = 60; // > For W-jets data driven estimation in 2jet2tag for ltau channels
    const double HighMtRegion_high = 120; // < For W-jets data driven estimation in 2jet2tag for ltau channels
}

namespace IsolationRegionForLeptonicChannel {
    const double pfRelIso = 0.1;
    const double isolation_low = 0.2; // > For QCD data driven estimation in 2jet*tag for ltau channels
    const double isolation_high = 0.5; // < For QCD data driven estimation in 2jet*tag for ltau channels
}

namespace jetCorrections {
    const double energyUncertainty = 0.05;
}

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
} // namespace skim
} // namespace cuts
