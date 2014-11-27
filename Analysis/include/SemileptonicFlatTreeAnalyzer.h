/*!
 * \file SemileptonicFlatTreeAnalyzer.h
 * \brief Definition of the base class for semi-leptonic flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-11-20 created
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

#include "BaseFlatTreeAnalyzer.h"

namespace analysis {

class SemileptonicFlatTreeAnalyzer : public BaseFlatTreeAnalyzer {
public:
    SemileptonicFlatTreeAnalyzer(const DataCategoryCollection& _dataCategoryCollection, const std::string& _inputPath,
                                 const std::string& _outputFileName)
         : BaseFlatTreeAnalyzer(_dataCategoryCollection, _inputPath, _outputFileName)
    {
    }

protected:
    bool IsHighMtRegion(const ntuple::Flat& event, analysis::EventCategory eventCategory)
    {
        using namespace cuts;
        if (eventCategory == analysis::EventCategory::TwoJets_TwoBtag)
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion_low &&
                    event.mt_1 < WjetsBackgroundEstimation::HighMtRegion_high;
        else
            return event.mt_1 > WjetsBackgroundEstimation::HighMtRegion;
    }

    virtual PhysicalValue CalculateQCDScaleFactor(const std::string& hist_name, EventCategory eventCategory) override
    {
        static const PhysicalValue sf(1.06, 0.001);
        return sf;
    }

    virtual void EstimateQCD(const std::string& hist_name, EventCategory eventCategory,
                             const PhysicalValue& scale_factor) override
    {
        return EstimateQCDEx(hist_name, eventCategory, eventCategory, EventRegion::SS_Isolated, scale_factor);
    }

    virtual PhysicalValueMap CalculateWjetsScaleFactors(EventCategory eventCategory,
                                                        const std::string& hist_name) override
    {
        PhysicalValueMap valueMap;
        static const analysis::PhysicalValue one(1, 0.001);
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);
        analysis::EventCategory refEventCategory = eventCategory;
        if(CategoriesToRelax().count(eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(eventCategory);

        for (const auto& eventRegion : HighMt_LowMt_RegionMap){
            auto hist_data = GetHistogram(eventCategory, data.name, eventRegion.first, hist_name);
            if(!hist_data)
                throw exception("Unable to find data histograms for Wjet scale factors estimation");
            TH1D& hist_HighMt = CloneHistogram(eventCategory, wjets.name, eventRegion.first, *hist_data);
            SubtractBackgroundHistograms(hist_HighMt, eventCategory, eventRegion.first, wjets.name, true);
            const PhysicalValue n_HighMt = Integral(hist_HighMt, false);

            PhysicalValue n_HighMt_mc;
            bool hist_mc_found = false;
            PhysicalValue ratio_HighToLowMt;
            for(const analysis::DataCategory* wjets_category : wjets_mc_categories){
                if(auto hist_mc = GetHistogram(eventCategory, wjets_category->name, eventRegion.first, hist_name)) {
                    hist_mc_found = true;
                    n_HighMt_mc += Integral(*hist_mc, false);
                }
            }
            try {
                if (!hist_mc_found)
                    throw exception("Unable to find mc histograms for Wjet scale factors estimation.");
                if(n_HighMt.value < 0)
                    throw exception("Negative number of estimated events in Wjets SF estimation for ")
                            << eventCategory << " " << eventRegion.second << ".";
                ratio_HighToLowMt = n_HighMt / n_HighMt_mc;
                //valueMap[eventRegion.second] = n_HighMt / n_HighMt_mc;
            } catch(exception& ex) {
                static const PhysicalValue defaultValue(1, 0.0001);
                std::cerr << ex.message() << " Using default value " << defaultValue << "." << std::endl;
                //valueMap[eventRegion.second] = defaultValue;
                ratio_HighToLowMt = defaultValue;
            }

            analysis::PhysicalValue sf_MediumToLow = one;
            if(eventCategory != refEventCategory) {
                try {
                    const auto n_medium = CalculateFullIntegral(eventCategory, eventRegion.first, hist_name, wjets_mc_categories);
                    const auto n_ref = CalculateFullIntegral(refEventCategory, eventRegion.first, hist_name, wjets_mc_categories);
                    sf_MediumToLow = n_medium / n_ref;
                } catch(analysis::exception& ex) {
                    std::cerr << ex.message() << " Using default value " << one << "." << std::endl;
                }
            }
            valueMap[eventRegion.second] = ratio_HighToLowMt*sf_MediumToLow;
        }

        return valueMap;
    }
};

} // namespace analysis
