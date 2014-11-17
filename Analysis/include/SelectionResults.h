/*!
 * \file SelectionResults.h
 * \brief Definition of selection results container used in HH->bbTauTau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-11-17 created
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

#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/Candidate.h"

#include "SVfit.h"
#include "KinFit.h"

#define SELECTION_ENTRY(name) \
    ENTRY_1D(cuts::ObjectSelector, name) \
    /**/

#define X(name) \
    selectionManager.FillHistogram(object.name, #name)

#define Y(name) \
    selectionManager.FillHistogram(name, #name)


namespace analysis {
class SelectionManager {
public:
    SelectionManager(root_ext::AnalyzerData& _anaData, const std::string& _selection_label, double _weight)
        : anaData(&_anaData), selection_label(_selection_label), weight(_weight) {}

    template<typename ValueType>
    ValueType FillHistogram(ValueType value, const std::string& histogram_name)
    {
        auto& histogram = anaData->Get(static_cast<TH1D*>(nullptr), histogram_name, selection_label);
        return cuts::fill_histogram(value, histogram, weight);
    }

private:
    root_ext::AnalyzerData* anaData;
    std::string selection_label;
    double weight;
};

struct SelectionResults {
    static constexpr size_t NumberOfLegs = 2;

    virtual ~SelectionResults() {}
    CandidatePtr higgs;
    sv_fit::CombinedFitResults svfitResults;
    kinematic_fit::four_body::FitResults kinfitResults;
    CandidatePtrVector jets;
    CandidatePtrVector jetsPt20;
    CandidatePtrVector bjets_all;
    CandidatePtrVector retagged_bjets;
    VertexPtrVector vertices;
    ntuple::MET pfMET;
    ntuple::MET MET_with_recoil_corrections;
    ntuple::EventType eventType;

    VertexPtr GetPrimaryVertex() const { return vertices.front(); }
    virtual CandidatePtr GetLeg(size_t leg_id) const = 0;
    virtual const finalState::bbTauTau& GetFinalStateMC() const = 0;
};
} // namespace analysis
