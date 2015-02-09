/*!
 * \file FlatTreeAnalyzer_tautau.C
 * \brief Analyze flat-tree for tau-tau channel for HHbbtautau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-03 created
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

#include "Analysis/include/BaseFlatTreeAnalyzer.h"

class FlatTreeAnalyzer_tautau : public analysis::BaseFlatTreeAnalyzer {
public:
    FlatTreeAnalyzer_tautau(const std::string& source_cfg, const std::string& _inputPath,
                            const std::string& outputFileName, const std::string& signal_list)
          : BaseFlatTreeAnalyzer(analysis::DataCategoryCollection(source_cfg, signal_list, ChannelId()), _inputPath,
                                 outputFileName)
    {
    }

protected:
    virtual analysis::Channel ChannelId() const override { return analysis::Channel::TauTau; }

    virtual analysis::EventRegion DetermineEventRegion(const ntuple::Flat& event,
                                                       analysis::EventCategory /*eventCategory*/) override
    {
        using analysis::EventRegion;
        using namespace cuts::Htautau_Summer13::TauTau::tauID;

        if(!event.againstElectronLooseMVA_2
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= BackgroundEstimation::Isolation_upperLimit
                || event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= BackgroundEstimation::Isolation_upperLimit
                || (event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 >= byCombinedIsolationDeltaBetaCorrRaw3Hits
                    && event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 >= byCombinedIsolationDeltaBetaCorrRaw3Hits)
                /*|| event.kinfit_bb_tt_chi2 > 20*/)
            return EventRegion::Unknown;

        const bool os = event.q_1 * event.q_2 == -1;
        const bool iso = event.byCombinedIsolationDeltaBetaCorrRaw3Hits_1 < byCombinedIsolationDeltaBetaCorrRaw3Hits &&
                         event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < byCombinedIsolationDeltaBetaCorrRaw3Hits;

        if(iso) return os ? EventRegion::OS_Isolated : EventRegion::SS_Isolated;
        return os ? EventRegion::OS_AntiIsolated : EventRegion::SS_AntiIsolated;
    }

    virtual analysis::PhysicalValue CalculateQCDYield(
            const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, const std::string& hist_name,
            analysis::DataCategoryType dataCategoryType, std::ostream& s_out) override
    {
        if (dataCategoryType == analysis::DataCategoryType::QCD)
            return CalculateQCDYield_base(anaDataMetaId, hist_name, s_out);
        else
            return CalculateQCDYield_alternative(anaDataMetaId, hist_name, s_out);
    }

    virtual void EstimateQCD(const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                             const std::string& hist_name, const analysis::PhysicalValue& yield,
                             analysis::DataCategoryType dataCategoryType) override
    {
        using analysis::EventCategory;
        using analysis::DataCategory;
        using analysis::DataCategoryType;
        using analysis::EventRegion;

        static const analysis::EventCategorySet categories = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag
        };
        static const std::map<DataCategoryType, EventRegion> dataCategoryType_eventRegions= {
            { DataCategoryType::QCD, EventRegion::OS_AntiIsolated },
            { DataCategoryType::QCD_alternative, EventRegion::SS_Isolated }
        };

        const DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(dataCategoryType);
        const DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        EventCategory refEventCategory = anaDataMetaId.eventCategory;
        if(categories.count(anaDataMetaId.eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory);

        EventRegion eventRegion = EventRegion::OS_AntiIsolated;
        if(dataCategoryType_eventRegions.count(dataCategoryType))
            eventRegion = dataCategoryType_eventRegions.at(dataCategoryType);

        const analysis::FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_ref(refEventCategory,
                                                                                 anaDataMetaId.eventSubCategory,
                                                                                 anaDataMetaId.eventEnergyScale);

        auto hist_shape_data = GetHistogram(anaDataMetaId_ref, eventRegion, data.name, hist_name);
        if(!hist_shape_data) return;

        TH1D& histogram = CloneHistogram(anaDataMetaId, EventRegion::OS_Isolated, qcd.name, *hist_shape_data);
        std::string debug_info, negative_bins_info;
        SubtractBackgroundHistograms(anaDataMetaId_ref, eventRegion, histogram, qcd.name, debug_info,
                                     negative_bins_info);
        if(negative_bins_info.size())
            std::cerr << negative_bins_info;
        analysis::RenormalizeHistogram(histogram, yield, true);
    }

    virtual analysis::PhysicalValue CalculateQCDYield_base(
            const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, const std::string& hist_name,
            std::ostream& s_out)
    {
        using analysis::EventCategory;
        using analysis::PhysicalValue;
        using analysis::EventRegion;

        static const analysis::EventCategorySet categories = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag
        };

        EventCategory refEventCategory = anaDataMetaId.eventCategory;
        if(categories.count(anaDataMetaId.eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory);

        const analysis::FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_ref(refEventCategory,
                                                                                 anaDataMetaId.eventSubCategory,
                                                                                 anaDataMetaId.eventEnergyScale);

        const PhysicalValue yield_OSAntiIso =
                CalculateYieldsForQCD(anaDataMetaId, EventRegion::OS_AntiIsolated, hist_name, s_out);

        const PhysicalValue yield_SSIso =
                CalculateYieldsForQCD(anaDataMetaId_ref, EventRegion::SS_Isolated, hist_name, s_out);

        const analysis::PhysicalValue yield_SSAntiIso =
                CalculateYieldsForQCD(anaDataMetaId_ref, EventRegion::SS_AntiIsolated, hist_name, s_out);

        const auto iso_antiIso_sf = yield_SSIso / yield_SSAntiIso;
        s_out << anaDataMetaId.eventCategory << ": QCD SF Iso / AntiIso = " << iso_antiIso_sf << ". \n";
        return yield_OSAntiIso * iso_antiIso_sf;
    }

    virtual analysis::PhysicalValue CalculateQCDYield_alternative(
            const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId, const std::string& hist_name,
            std::ostream& s_out)
    {
        using analysis::EventCategory;
        using analysis::PhysicalValue;
        using analysis::EventRegion;

        static const analysis::EventCategorySet categories = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag
        };

        EventCategory refEventCategory = anaDataMetaId.eventCategory;
        if(categories.count(anaDataMetaId.eventCategory))
            refEventCategory = analysis::MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory);

        const analysis::FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_ref(refEventCategory,
                                                                                 anaDataMetaId.eventSubCategory,
                                                                                 anaDataMetaId.eventEnergyScale);

        const PhysicalValue yield_SSIso =
                CalculateYieldsForQCD(anaDataMetaId, EventRegion::SS_Isolated, hist_name, s_out);

        const PhysicalValue yield_OSAntiIso =
                CalculateYieldsForQCD(anaDataMetaId_ref, EventRegion::OS_AntiIsolated, hist_name, s_out);

        const PhysicalValue yield_SSAntiIso =
                CalculateYieldsForQCD(anaDataMetaId_ref, EventRegion::SS_AntiIsolated, hist_name, s_out);

        const auto os_ss_antiIso_sf = yield_OSAntiIso / yield_SSAntiIso;
        s_out << anaDataMetaId.eventCategory << ": QCD alternative SF OS_AntiIso / SS_AntiIso = "
              << os_ss_antiIso_sf << ". \n";

        return yield_SSIso * os_ss_antiIso_sf;
    }

    virtual PhysicalValueMap CalculateWjetsYields(const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                                  const std::string& hist_name, bool fullEstimate) override
    {
        using analysis::EventCategory;

        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(analysis::DataCategoryType::WJets_MC);

        static const analysis::EventCategorySet NotEstimated_categories = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_OneLooseBtag,
            EventCategory::TwoJets_TwoLooseBtag, EventCategory::TwoJets_AtLeastOneBtag,
            EventCategory::TwoJets_AtLeastOneLooseBtag
        };

        PhysicalValueMap valueMap;
        if (fullEstimate || !NotEstimated_categories.count(anaDataMetaId.eventCategory)) {
            for (analysis::EventRegion eventRegion : analysis::QcdRegions)
                valueMap[eventRegion] = CalculateFullIntegral(anaDataMetaId, eventRegion, hist_name,
                                                              wjets_mc_categories);
        }
        return valueMap;
    }

    virtual void CreateHistogramForZTT(const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                       const std::string& hist_name, const PhysicalValueMap& ztt_yield_map,
                                       bool useEmbedded) override
    {
        using analysis::DataCategory;
        using analysis::DataCategoryType;
        using analysis::EventCategory;

        const DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const DataCategory& ZTT = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT);
        const DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);
        const DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::TT_Embedded);

        static const analysis::EventCategorySet categoriesToRelax =
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };

        EventCategory shapeEventCategory = anaDataMetaId.eventCategory;
        if (categoriesToRelax.count(anaDataMetaId.eventCategory))
            shapeEventCategory = analysis::MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory);

        const analysis::FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_shape(shapeEventCategory,
                                                                                   anaDataMetaId.eventSubCategory,
                                                                                   anaDataMetaId.eventEnergyScale);

        for(const auto& eventRegionKey : ztt_yield_map) {
            const analysis::EventRegion eventRegion = eventRegionKey.first;
            const analysis::PhysicalValue& ztt_yield = eventRegionKey.second;
            auto ztt_l_hist = GetHistogram(anaDataMetaId, eventRegion, ZTT_L.name, hist_name);
            const std::string embeddedName = useEmbedded && eventRegion == analysis::EventRegion::OS_Isolated ?
                                             embedded.name : ZTT_MC.name;

            auto embedded_hist = GetHistogram(anaDataMetaId_shape, eventRegion, embeddedName, hist_name);
            auto TTembedded_hist = GetHistogram(anaDataMetaId_shape, eventRegion, TTembedded.name, hist_name);

            if (embedded_hist){
                TH1D& ztt_hist = CloneHistogram(anaDataMetaId, eventRegion, ZTT.name, *embedded_hist);
                if (TTembedded_hist && useEmbedded)
                    ztt_hist.Add(TTembedded_hist, -1);
                analysis::RenormalizeHistogram(ztt_hist, ztt_yield, true);
                if (ztt_l_hist)
                    ztt_hist.Add(ztt_l_hist);                
            }
            if (!embedded_hist && ztt_l_hist)
                CloneHistogram(anaDataMetaId, eventRegion, ZTT.name, *ztt_l_hist);
        }
    }

    virtual void CreateHistogramForVVcategory(const analysis::FlatAnalyzerDataMetaId_noRegion_noName& anaDataMetaId,
                                              const std::string& hist_name) override
    {
        using analysis::DataCategoryType;
        using analysis::DataCategory;
        using analysis::EventRegion;
        using analysis::EventCategory;

        static const std::map<DataCategoryType, DataCategoryType> diboson_category_map = {
            { DataCategoryType::DiBoson_MC, DataCategoryType::DiBoson }
        };

        static const analysis::EventCategorySet categoriesToRelax = {
            EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag
        };

        for (const auto& diboson_category : diboson_category_map) {
            const DataCategory& originalVVcategory = dataCategoryCollection.GetUniqueCategory(diboson_category.first);
            const DataCategory& newVVcategory = dataCategoryCollection.GetUniqueCategory(diboson_category.second);

            PhysicalValueMap valueMap;

            for(EventRegion eventRegion : analysis::AllEventRegions) {
                auto vv_hist_yield = GetHistogram(anaDataMetaId, eventRegion, originalVVcategory.name, hist_name);
                if (vv_hist_yield)
                    valueMap[eventRegion] = analysis::Integral(*vv_hist_yield, true);
            }

            const EventCategory shapeEventCategory = categoriesToRelax.count(anaDataMetaId.eventCategory)
                    ? analysis::MediumToLoose_EventCategoryMap.at(anaDataMetaId.eventCategory)
                    : anaDataMetaId.eventCategory;
            const analysis::FlatAnalyzerDataMetaId_noRegion_noName anaDataMetaId_shape(shapeEventCategory,
                                                                                       anaDataMetaId.eventSubCategory,
                                                                                       anaDataMetaId.eventEnergyScale);

            for(const auto& yield_iter : valueMap) {
                const analysis::EventRegion eventRegion = yield_iter.first;
                const analysis::PhysicalValue& yield = yield_iter.second;
                auto vv_hist_shape = GetHistogram(anaDataMetaId_shape, eventRegion, originalVVcategory.name, hist_name);
                if (vv_hist_shape){
                    TH1D& vv_hist = CloneHistogram(anaDataMetaId, eventRegion, newVVcategory.name, *vv_hist_shape);
                    RenormalizeHistogram(vv_hist, yield, true);
                }
            }
        }
    }
};
