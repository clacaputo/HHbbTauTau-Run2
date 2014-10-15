/*!
 * \file FlatEventInfo.h
 * \brief Definiton of analysis::FlatEventInfo class.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-07 created
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
#include <TMatrixD.h>

#include "FlatTree.h"
#include "AnalysisTypes.h"
#include "exception.h"

namespace analysis {

struct FlatEventInfo {
    typedef std::pair<size_t, size_t> BjetPair;

    static size_t CombinationPairToIndex(const BjetPair& pair, size_t n_bjets)
    {
        const size_t min = std::min(pair.first, pair.second);
        const size_t max = std::max(pair.first, pair.second);
        if(n_bjets < 2 || min == max || min >= n_bjets || max >= n_bjets)
            throw exception("bad combination pair (") << pair.first << ", " << pair.second
                                                      << ") for n b-jets = " << n_bjets << ".";
        return max - 1 + min * (2 * n_bjets - 3 - min) / 2;
    }

    static BjetPair CombinationIndexToPair(size_t index, size_t n_bjets)
    {
        if(n_bjets < 2 || index >= n_bjets * (n_bjets - 1) / 2)
            throw exception("bad combination index = ") << index << " for n b-jets = " << n_bjets << ".";

        for(size_t min = 0;; ++min) {
            const size_t l = CombinationPairToIndex(BjetPair(min, n_bjets - 1), n_bjets);
            if(l >= index) {
                const size_t max = index + n_bjets - 1 - l;
                return BjetPair(min, max);
            }
        }
    }

    const ntuple::Flat* event;
    ntuple::EventType eventType;
    analysis::Channel channel;
    std::vector<TLorentzVector> lepton_momentums;
    std::vector<TLorentzVector> bjet_momentums;
    BjetPair selected_bjets;
    bool has_bjet_pair;
    TLorentzVector MET, Htt, Htt_MET, Hbb, resonance;
    TMatrixD MET_covariance;
    size_t kinfit_data_index;
    double mva_BDT, mva_BDTD, mva_BDTMitFisher;

    FlatEventInfo(const ntuple::Flat& _event, const BjetPair& _selected_bjets)
        : event(&_event), eventType(static_cast<ntuple::EventType>(_event.eventType)),
          channel(static_cast<analysis::Channel>(_event.channel)),
          lepton_momentums(2), bjet_momentums(_event.pt_Bjets.size()), selected_bjets(_selected_bjets),
          has_bjet_pair(false), MET_covariance(2, 2), kinfit_data_index(std::numeric_limits<size_t>::max()),
          mva_BDT(-1), mva_BDTD(-1), mva_BDTMitFisher(-1)
    {
        lepton_momentums.at(0).SetPtEtaPhiM(event->pt_1, event->eta_1, event->phi_1, event->m_1);
        lepton_momentums.at(1).SetPtEtaPhiM(event->pt_2, event->eta_2, event->phi_2, event->m_2);

        for(size_t n = 0; n < bjet_momentums.size(); ++n)
            bjet_momentums.at(n).SetPtEtaPhiE(event->pt_Bjets.at(n), event->eta_Bjets.at(n), event->phi_Bjets.at(n),
                                              event->energy_Bjets.at(n));

        MET.SetPtEtaPhiM(event->mvamet, 0, event->mvametphi, 0);
        MET_covariance(0, 0) = event->mvacov00;
        MET_covariance(1, 0) = event->mvacov10;
        MET_covariance(0, 1) = event->mvacov01;
        MET_covariance(1, 1) = event->mvacov11;

        Htt = lepton_momentums.at(0) + lepton_momentums.at(1);
        Htt_MET = Htt + MET;
        has_bjet_pair = selected_bjets.first < bjet_momentums.size() && selected_bjets.second < bjet_momentums.size();
        if(has_bjet_pair) {
            Hbb = bjet_momentums.at(selected_bjets.first) + bjet_momentums.at(selected_bjets.second);
            resonance = Htt_MET + Hbb;
            kinfit_data_index = CombinationPairToIndex(selected_bjets, bjet_momentums.size());
        }
    }
};

} // namespace analysis
