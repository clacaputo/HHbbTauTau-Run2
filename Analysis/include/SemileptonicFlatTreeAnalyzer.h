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
    SemileptonicFlatTreeAnalyzer(const std::string& source_cfg, const std::string& hist_cfg,
                                 const std::string& _inputPath, const std::string& outputFileName,
                                 Channel channel_id, const std::string& signal_list)
         : BaseFlatTreeAnalyzer(source_cfg, hist_cfg, _inputPath, outputFileName, channel_id, signal_list)
    {
    }

protected:
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
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        for (const auto& eventRegion : HighMt_LowMt_RegionMap){
            auto hist_data = GetHistogram(eventCategory, data.name, eventRegion.first, hist_name);
            if(!hist_data)
                throw exception("Unable to find data histograms for Wjet scale factors estimation");
            TH1D& hist_HighMt = CloneHistogram(eventCategory, wjets.name, eventRegion.first, *hist_data);
            SubtractBackgroundHistograms(hist_HighMt, eventCategory, eventRegion.first, wjets.name, true);
            const PhysicalValue n_HighMt = Integral(hist_HighMt, false);

            PhysicalValue n_HighMt_mc;
            bool hist_mc_found = false;
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
                valueMap[eventRegion.second] = n_HighMt / n_HighMt_mc;
            } catch(exception& ex) {
                static const PhysicalValue defaultValue(1, 0.0001);
                std::cerr << ex.message() << " Using default value " << defaultValue << "." << std::endl;
                valueMap[eventRegion.second] = defaultValue;
            }
        }

        return valueMap;
    }
};

} // namespace analysis
