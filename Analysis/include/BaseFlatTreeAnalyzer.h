/*!
 * \file BaseFlatTreeAnalyzer.h
 * \brief Definition of BaseFlatTreeAnalyzer class, the base class for flat tree analyzers.
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

#pragma once

#ifndef __APPLE__
#define override
#endif

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>

#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalysisBase/include/FlatEventInfo.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/AnalysisTypes.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"
#include "FlatAnalyzerData.h"

namespace analysis {

class BaseFlatTreeAnalyzer {
public:
    typedef std::map<EventRegion, std::shared_ptr<FlatAnalyzerData> > RegionAnaData;

    typedef std::map<std::string, RegionAnaData> AnaDataForEventCategory;
    typedef std::map<EventCategory, AnaDataForEventCategory> FullAnaData;
    typedef std::map<EventRegion, PhysicalValue> PhysicalValueMap;

    static const EventCategorySet& EventCategoriesToProcess() { return AllEventCategories; }

    BaseFlatTreeAnalyzer(const DataCategoryCollection& _dataCategoryCollection, const std::string& _inputPath,
                         const std::string& _outputFileName)
        : inputPath(_inputPath), outputFileName(_outputFileName), dataCategoryCollection(_dataCategoryCollection)
    {
        TH1::SetDefaultSumw2();
    }

    void Run()
    {
        std::cout << "Processing data categories... " << std::endl;
        for(const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            if(!dataCategory->sources_sf.size()) continue;
            std::cout << *dataCategory << std::endl;
            for(const auto& source_entry : dataCategory->sources_sf) {
                const std::string fullFileName = inputPath + "/" + source_entry.first;
                std::shared_ptr<TFile> file(new TFile(fullFileName.c_str(), "READ"));
                if(file->IsZombie())
                    throw exception("Input file '") << source_entry.first << "' not found.";
                std::shared_ptr<ntuple::FlatTree> tree(new ntuple::FlatTree(*file, "flatTree"));
                ProcessDataSource(*dataCategory, tree, source_entry.second);
            }
        }

        static const std::set<std::string> interesting_histograms = {
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_sv_Name(), EventSubCategory::MassWindow,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                                EventSubCategory::KinematicFitConverged,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                                EventSubCategory::KinematicFitConvergedWithMassWindow,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_bb_slice_Name(), EventSubCategory::NoCuts,
                                                EventEnergyScale::Central)
        };

        for (const auto& hist_name : FlatAnalyzerData::GetAllHistogramNames()) {
            std::cout << "Processing '" << hist_name << "'..." << std::endl;

            std::ostringstream ss_debug;
            std::ostream& s_out = interesting_histograms.count(hist_name) ? std::cout : ss_debug;

            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                const auto ZTT_matched_yield = CalculateZTTmatchedYield(hist_name,eventCategory,true);
                for (const auto yield_entry : ZTT_matched_yield)
                    s_out << eventCategory << ": ZTT MC yield in " << yield_entry.first << " = " << yield_entry.second
                          << ".\n";

                CreateHistogramForZTT(eventCategory, hist_name, ZTT_matched_yield, true);
                CreateHistogramForZcategory(eventCategory,hist_name);
                CreateHistogramForVVcategory(eventCategory,hist_name);
            }

            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                const auto wjets_yields = CalculateWjetsYields(eventCategory, hist_name,false);
                for (const auto yield_entry : wjets_yields){
                    s_out << eventCategory << ": W+jets yield in " << yield_entry.first << " = " << yield_entry.second
                          << ".\n";
                }
                EstimateWjets(eventCategory, hist_name, wjets_yields);
            }

            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                const auto qcd_yield = CalculateQCDYield(hist_name, eventCategory, s_out);
                s_out << eventCategory << ": QCD yield = " << qcd_yield << ".\n";
                EstimateQCD(hist_name, eventCategory, qcd_yield);
                ProcessCompositDataCategories(eventCategory, hist_name);
            }
            s_out << std::endl;
        }

        std::cout << "\nSaving tables... " << std::endl;
        PrintTables("comma", L",");
        PrintTables("semicolon", L";");

        std::cout << "Saving datacards... " << std::endl;
        ProduceFileForLimitsCalculation(FlatAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts);
        ProduceFileForLimitsCalculation(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                        EventSubCategory::KinematicFitConverged);
        ProduceFileForLimitsCalculation(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                        EventSubCategory::KinematicFitConvergedWithMassWindow);
        ProduceFileForLimitsCalculation(FlatAnalyzerData::m_bb_slice_Name(), EventSubCategory::NoCuts);

        std::cout << "Printing stacked plots... " << std::endl;
        PrintStackedPlots(true, true, EventRegion::OS_Isolated);
    }

protected:
    virtual Channel ChannelId() const = 0;
    virtual EventRegion DetermineEventRegion(const ntuple::Flat& event, EventCategory eventCategory) = 0;
    virtual bool PassMvaCut(const FlatEventInfo& eventInfo, EventCategory eventCategory) { return true; }

    virtual PhysicalValue CalculateQCDYield(const std::string& hist_name, EventCategory eventCategory,
                                            std::ostream& s_out) = 0;
    virtual void EstimateQCD(const std::string& hist_name, EventCategory eventCategory,
                             const PhysicalValue& scale_factor) = 0;
    virtual PhysicalValueMap CalculateWjetsYields(EventCategory eventCategory, const std::string& hist_name,
                                                  bool fullEstimate) = 0;
    virtual void CreateHistogramForZTT(EventCategory eventCategory, const std::string& hist_name,
                               const PhysicalValueMap& ztt_yield, bool useEmbedded) = 0;
    virtual void CreateHistogramForVVcategory(EventCategory eventCategory,const std::string& hist_name) = 0;

    PhysicalValueMap CalculateZTTmatchedYield(const std::string& hist_name,
                                                             analysis::EventCategory eventCategory,
                                                             bool useEmbedded)
    {

        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::ZTT_MC);

        PhysicalValueMap zttYield;
        if (useEmbedded){
            const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Embedded);
            const analysis::DataCategory& TTembedded = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::TT_Embedded);
            TH1D* hist_embedded_inclusive = GetSignalHistogram(analysis::EventCategory::Inclusive, embedded.name, hist_name);
            TH1D* hist_embedded_category = GetSignalHistogram(eventCategory, embedded.name, hist_name);
            TH1D* hist_TTembedded_inclusive = GetSignalHistogram(analysis::EventCategory::Inclusive, TTembedded.name, hist_name);
            TH1D* hist_TTembedded_category = GetSignalHistogram(eventCategory, TTembedded.name, hist_name);
            TH1D* hist_ztautau_inclusive = GetSignalHistogram(analysis::EventCategory::Inclusive, ZTT_MC.name, hist_name);
            if(!hist_embedded_inclusive)
                throw analysis::exception("embedded hist in inclusive category not found");
            if(!hist_embedded_category)
                throw analysis::exception("embedded hist not found in event category: ") << eventCategory << "\n";
            if(!hist_TTembedded_inclusive)
                throw analysis::exception("TTembedded hist in inclusive category not found");
            if(!hist_TTembedded_category)
                throw analysis::exception("TTembedded hist not found in event category: ") << eventCategory << "\n";
            if(!hist_ztautau_inclusive )
                throw analysis::exception("ztt hist in inclusive category not found");

            const analysis::PhysicalValue n_emb_inclusive =
                    analysis::Integral(*hist_embedded_inclusive, true) - analysis::Integral(*hist_TTembedded_inclusive, true);
            const analysis::PhysicalValue n_emb_category =
                    analysis::Integral(*hist_embedded_category, true) - analysis::Integral(*hist_TTembedded_category, true);
            const analysis::PhysicalValue n_ztautau_inclusive = analysis::Integral(*hist_ztautau_inclusive, true);
            const analysis::PhysicalValue embedded_eff = n_emb_category/n_emb_inclusive;
            zttYield[analysis::EventRegion::OS_Isolated] = n_ztautau_inclusive * embedded_eff;
        }

        for (analysis::EventRegion eventRegion : analysis::AllEventRegions){
            if (zttYield.count(eventRegion)) continue;
            TH1D* hist_ztautau = GetHistogram(eventCategory, ZTT_MC.name,eventRegion, hist_name);

            if(!hist_ztautau ){
                if (eventRegion == analysis::EventRegion::OS_Isolated)
                    throw analysis::exception("ztt hist not found in event category") << eventCategory;
                continue;
            }

            zttYield[eventRegion] = analysis::Integral(*hist_ztautau, true);
        }
        return zttYield;
    }

    virtual void CreateHistogramForZcategory(EventCategory eventCategory, const std::string& hist_name)
    {
        const std::map<DataCategoryType, DataCategoryType> z_type_category_map = {
            { DataCategoryType::ZL_MC, DataCategoryType::ZL }, { DataCategoryType::ZJ_MC, DataCategoryType::ZJ }
        };

        for (const auto& z_category : z_type_category_map){
            const analysis::DataCategory& originalZcategory = dataCategoryCollection.GetUniqueCategory(z_category.first);
            const analysis::DataCategory& newZcategory = dataCategoryCollection.GetUniqueCategory(z_category.second);

            PhysicalValueMap valueMap;

            for(EventRegion eventRegion : AllEventRegions) {
                auto z_hist_yield = GetHistogram(eventCategory, originalZcategory.name, eventRegion, hist_name);
                if (z_hist_yield)
                    valueMap[eventRegion] = Integral(*z_hist_yield,true);
            }

            static const EventCategorySet categoriesToRelax =
                { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
            const EventCategory shapeEventCategory = categoriesToRelax.count(eventCategory)
                    ? analysis::MediumToLoose_EventCategoryMap.at(eventCategory) : eventCategory;

            for(const auto& yield_iter : valueMap) {
                const EventRegion eventRegion = yield_iter.first;
                const PhysicalValue& yield = yield_iter.second;
                auto z_hist_shape = GetHistogram(shapeEventCategory, originalZcategory.name, eventRegion, hist_name);
                if (z_hist_shape){
                    TH1D& z_hist = CloneHistogram(eventCategory, newZcategory.name, eventRegion, *z_hist_shape);
                    RenormalizeHistogram(z_hist,yield,true);
                }
            }

        }
    }

    PhysicalValue CalculateYieldsForQCD(const std::string& hist_name, EventCategory eventCategory,
                                        EventRegion eventRegion, std::ostream& s_out)
    {
        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(analysis::DataCategoryType::Data);

        std::string bkg_yield_debug;
        const analysis::PhysicalValue bkg_yield =
                CalculateBackgroundIntegral(hist_name, eventCategory, eventRegion, qcd.name, false, bkg_yield_debug);

        s_out << bkg_yield_debug;

        auto hist_data = GetHistogram(eventCategory, data.name, eventRegion, hist_name);
        if(!hist_data)
            throw analysis::exception("Unable to find data histograms for QCD yield estimation.");
        const auto data_yield = analysis::Integral(*hist_data, true);
        const analysis::PhysicalValue yield = data_yield - bkg_yield;
        s_out << "Data yield = " << data_yield << "\nData-MC yield = " << yield << std::endl;
        if(yield.value < 0) {
            std::cout << bkg_yield_debug << "\nData yield = " << data_yield << std::endl;
            throw exception("Negative QCD yield for histogram '") << hist_name << "' in " << eventCategory << " "
                                                                  << eventRegion << ".";
        }
        return yield;
    }

    virtual void EstimateWjets(EventCategory eventCategory, const std::string& hist_name,
                         const PhysicalValueMap& yield_map)
    {
        static const EventCategorySet categoriesToRelax =
            { EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag, EventCategory::TwoJets_AtLeastOneBtag };
        const EventCategory shapeEventCategory = categoriesToRelax.count(eventCategory)
                ? analysis::MediumToLoose_EventCategoryMap.at(eventCategory) : eventCategory;
        return EstimateWjetsEx(eventCategory,shapeEventCategory,hist_name,yield_map);
    }

    void EstimateWjetsEx(EventCategory eventCategory, EventCategory shapeEventCategory, const std::string& hist_name,
                         const PhysicalValueMap& yield_map)
    {
        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategoryPtrSet& wjets_mc_categories =
                dataCategoryCollection.GetCategories(DataCategoryType::WJets_MC);

        for(const auto& yield_entry : yield_map) {
            const EventRegion eventRegion = yield_entry.first;
            const PhysicalValue& yield = yield_entry.second;
            TH1D* hist = nullptr;
            for (const analysis::DataCategory* wjets_category : wjets_mc_categories){
                if(auto hist_mc = GetHistogram(shapeEventCategory, wjets_category->name, eventRegion, hist_name)) {
                    if (!hist)
                        hist = &CloneHistogram(eventCategory, wjets.name, eventRegion, *hist_mc);
                    else
                        hist->Add(hist_mc);
                }
            }
            if (hist)
                RenormalizeHistogram(*hist,yield,true);
        }
    }

    const std::string& ChannelName() const { return detail::ChannelNameMap.at(ChannelId()); }
    const std::string& ChannelNameLatex() const { return detail::ChannelNameMapLatex.at(ChannelId()); }

    virtual std::shared_ptr<FlatAnalyzerData> MakeAnaData(EventCategory eventCategory)
    {
        return std::shared_ptr<FlatAnalyzerData>(new FlatAnalyzerData());
    }

    FlatAnalyzerData& GetAnaData(EventCategory eventCategory, const std::string& dataCategoryName,
                                 EventRegion eventRegion)
    {
        std::shared_ptr<FlatAnalyzerData>& anaData = fullAnaData[eventCategory][dataCategoryName][eventRegion];
        if(!anaData)
            anaData = MakeAnaData(eventCategory);
        return *anaData;
    }

    root_ext::SmartHistogram<TH1D>* GetHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                                                 EventRegion eventRegion, const std::string& histogramName)
    {
        return GetAnaData(eventCategory, dataCategoryName, eventRegion).GetPtr<TH1D>(histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetSignalHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                                                       const std::string& histogramName)
    {
        return GetHistogram(eventCategory, dataCategoryName, EventRegion::OS_Isolated, histogramName);
    }

    TH1D& CloneHistogram(EventCategory eventCategory, const std::string& dataCategoryName, EventRegion eventRegion,
                         const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return GetAnaData(eventCategory, dataCategoryName, eventRegion).Clone(originalHistogram);
    }

    TH1D& CloneSignalHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                         const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(eventCategory, dataCategoryName, EventRegion::OS_Isolated, originalHistogram);
    }


    void ProcessDataSource(const DataCategory& dataCategory, std::shared_ptr<ntuple::FlatTree> tree, double scale_factor)
    {
        static const bool applyMVAcut = false;

        const analysis::DataCategory& DYJets_incl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_incl);
        const analysis::DataCategory& DYJets_excl = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets_excl);
        const analysis::DataCategory& DY_Embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);

        for(Long64_t current_entry = 0; current_entry < tree->GetEntries(); ++current_entry) {
            tree->GetEntry(current_entry);
            const ntuple::Flat& event = tree->data;

            const bool useRetag = dataCategory.IsData() || dataCategory.name == DY_Embedded.name  ? false : true;
            const EventCategoryVector eventCategories = DetermineEventCategories(event.csv_Bjets,
                                                                                 event.nBjets_retagged,
                                                                                 cuts::Htautau_Summer13::btag::CSVL,
                                                                                 cuts::Htautau_Summer13::btag::CSVM,
                                                                                 useRetag);
            const bool useCustomSF = dataCategory.exclusive_sf.count(event.n_extraJets_MC);
            const double corrected_sf = useCustomSF
                    ? dataCategory.exclusive_sf.at(event.n_extraJets_MC) : scale_factor;

            if(std::isnan(event.weight)) {
                std::cerr << "ERROR: event " << event.evt << " will not be processed because event weight is NaN."
                          << std::endl;
                continue;
            }

            double weight;
            if (dataCategory.IsData())
                weight = 1;
            else if (!dataCategory.IsSignal() && static_cast<ntuple::EventType>(event.eventType) != ntuple::EventType::ZTT
                     && event.decayMode_2 == ntuple::tau_id::kOneProng0PiZero)
                weight = event.weight * corrected_sf / cuts::Htautau_Summer13::tauCorrections::DecayModeWeight;
            else
                weight = event.weight * corrected_sf;

            std::shared_ptr<FlatEventInfo> eventInfo;
            for(auto eventCategory : eventCategories) {
                if (!EventCategoriesToProcess().count(eventCategory)) continue;
                const EventRegion eventRegion = DetermineEventRegion(event, eventCategory);
                if(eventRegion == EventRegion::Unknown) continue;

                if(!eventInfo)
                    eventInfo = std::shared_ptr<FlatEventInfo>(new FlatEventInfo(event, FlatEventInfo::BjetPair(0, 1),
                                                                                 false));

                UpdateMvaInfo(*eventInfo, eventCategory, false, false, false);
                if(applyMVAcut && !PassMvaCut(*eventInfo, eventCategory)) continue;


                if(dataCategory.name == DYJets_excl.name || dataCategory.name == DYJets_incl.name)
                    FillDYjetHistograms(*eventInfo, eventCategory, eventRegion, weight);
                auto& anaData = GetAnaData(eventCategory, dataCategory.name, eventRegion);
                if (dataCategory.IsData())
                    anaData.FillAllEnergyScales(*eventInfo, weight);
                else if(dataCategory.name == DY_Embedded.name && eventInfo->eventEnergyScale == EventEnergyScale::Central)
                    anaData.FillCentralAndJetEnergyScales(*eventInfo, weight);
                else
                    anaData.FillSubCategories(*eventInfo, weight, eventInfo->eventEnergyScale);
            }
        }
    }

    void FillDYjetHistograms(const FlatEventInfo& eventInfo, EventCategory eventCategory, EventRegion eventRegion,
                             double weight)
    {
        const analysis::DataCategory& ZL_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZL_MC);
        const analysis::DataCategory& ZJ_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZJ_MC);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const analysis::DataCategory& ZTT_L = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_L);

        const std::map<ntuple::EventType, std::string> type_category_map = {
            { ntuple::EventType::ZL, ZL_MC.name }, { ntuple::EventType::ZJ, ZJ_MC.name },
            { ntuple::EventType::ZTT, ZTT_MC.name }, {ntuple::EventType::ZTT_L, ZTT_L.name}
        };

        if(type_category_map.count(eventInfo.eventType)) {
            const std::string& name = type_category_map.at(eventInfo.eventType);
            GetAnaData(eventCategory, name, eventRegion).FillSubCategories(eventInfo, weight,
                                                                           eventInfo.eventEnergyScale);
        }
    }

    void UpdateMvaInfo(FlatEventInfo& eventInfo, EventCategory eventCategory, bool calc_BDT, bool calc_BDTD,
                       bool calc_BDTMitFisher)
    {
        static double const default_value = std::numeric_limits<double>::lowest();
        std::ostringstream category_name;
        category_name << eventCategory;

        auto getMVA = [&](bool calc_MVA, MVA_Selections::MvaMethod method) -> double {
            if(calc_MVA) {
                auto mvaReader = MVA_Selections::MvaReader::Get(ChannelName(), category_name.str(), method);
                if(mvaReader)
                    return mvaReader->GetMva(eventInfo.lepton_momentums.at(0), eventInfo.lepton_momentums.at(1),
                                             eventInfo.bjet_momentums.at(0), eventInfo.bjet_momentums.at(1),
                                             eventInfo.MET);
            }
            return default_value;
        };

        eventInfo.mva_BDT = getMVA(calc_BDT, MVA_Selections::BDT);
        eventInfo.mva_BDTD = getMVA(calc_BDTD, MVA_Selections::BDTD);
        eventInfo.mva_BDTMitFisher = getMVA(calc_BDTMitFisher, MVA_Selections::BDTMitFisher);
    }

    void PrintStackedPlots(bool isBlind, bool drawRatio, EventRegion eventRegion)
    {
        const std::string blindCondition = isBlind ? "_blind" : "_noBlind";
        const std::string ratioCondition = drawRatio ? "_ratio" : "_noRatio";
        std::ostringstream eventRegionName;
        eventRegionName << outputFileName << blindCondition << ratioCondition << "_" << eventRegion << ".pdf";
        root_ext::PdfPrinter printer(eventRegionName.str());

        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            for (const auto& hist_name : FlatAnalyzerData::GetAllHistogramNames()) {
                //root_ext::PdfPrinter printer(outputFileName + blindCondition + "_" + hist.name + ratioCondition +".pdf");
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist_name;
                StackedPlotDescriptor stackDescriptor(ss_title.str(), false, ChannelNameLatex(), drawRatio);

                for(const DataCategory* category : dataCategoryCollection.GetAllCategories()) {
                    if(!category->draw) continue;

//                    const auto histogram = GetSignalHistogram(eventCategory, category->name, hist_name);
                    const auto histogram = GetHistogram(eventCategory,category->name,eventRegion,hist_name);
                    if(!histogram) continue;

                    if(category->IsSignal())
                        stackDescriptor.AddSignalHistogram(*histogram, category->title, category->color, category->draw_sf);
                    else if(category->IsBackground())
                        stackDescriptor.AddBackgroundHistogram(*histogram, category->title, category->color);
                    else if(category->IsData())
                        stackDescriptor.AddDataHistogram(*histogram, category->title, isBlind, GetBlindRegion(hist_name));
                }

                printer.PrintStack(stackDescriptor);
            }
        }
    }

    std::string FullDataCardName(const std::string& datacard_name, EventEnergyScale eventEnergyScale) const
    {
        if(eventEnergyScale == EventEnergyScale::Central)
            return datacard_name;

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
        std::ostringstream full_name;
        full_name << datacard_name << "_CMS_scale_";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::TauDown)
            full_name << "t_" << channel_name;
        else if(eventEnergyScale == EventEnergyScale::JetUp || eventEnergyScale == EventEnergyScale::JetDown)
            full_name << "j";
        else if(eventEnergyScale == EventEnergyScale::BtagEfficiencyUp || eventEnergyScale == EventEnergyScale::BtagEfficiencyDown)
            full_name << "btag";
        else
            throw exception("Unsupported event energy scale ") << eventEnergyScale;
        full_name << "_8TeV";
        if(eventEnergyScale == EventEnergyScale::TauUp || eventEnergyScale == EventEnergyScale::JetUp)
            full_name << "Up";
        else
            full_name << "Down";
        return full_name.str();
    }


    void ProduceFileForLimitsCalculation(const std::string& hist_name, EventSubCategory eventSubCategory)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" },
            { EventCategory::TwoJets_ZeroLooseBtag, "2jetloose0tag" },
            { EventCategory::TwoJets_OneLooseBtag, "2jetloose1tag" },
            { EventCategory::TwoJets_TwoLooseBtag, "2jetloose2tag" },
            { EventCategory::TwoJets_AtLeastOneBtag, "2jet_at_least_1tag" },
            { EventCategory::TwoJets_AtLeastOneLooseBtag, "2jet_at_least_loose1tag" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        const std::string file_name = outputFileName + "_"
                + FlatAnalyzerData::FullHistogramName(hist_name, eventSubCategory, EventEnergyScale::Central)
                + ".root";

        std::shared_ptr<TFile> outputFile(new TFile(file_name.c_str(), "RECREATE"));
        outputFile->cd();
        for(EventCategory eventCategory : EventCategoriesToProcess()) {
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '") << dataCategory->name << "'.";
                for(const EventEnergyScale& eventEnergyScale : AllEventEnergyScales) {
                    const std::string full_hist_name = FlatAnalyzerData::FullHistogramName(hist_name, eventSubCategory,
                                                                                           eventEnergyScale);
                    std::shared_ptr<TH1D> hist;
                    if(auto hist_orig = GetSignalHistogram(eventCategory, dataCategory->name, full_hist_name))
                        hist = std::shared_ptr<TH1D>(static_cast<TH1D*>(hist_orig->Clone()));
                    else {
                        std::cout << "Warning - Datacard histogram '" << full_hist_name
                                  << "' not found for data category '" << dataCategory->name << "' for eventCategory '"
                                  << categoryToDirectoryNameSuffix.at(eventCategory) << ".\n";
                        continue;
                    }
                    const std::string full_datacard_name = FullDataCardName(dataCategory->datacard, eventEnergyScale);
                    hist->Scale(dataCategory->limits_sf);
                    hist->Write(full_datacard_name.c_str());

                    if(eventEnergyScale == EventEnergyScale::Central && dataCategory->datacard == "ZL") {
                        std::string channel_name = ChannelName();
                        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);
                        const std::string name_syst_prefix = dataCategory->datacard + "_CMS_htt_" + dataCategory->datacard
                                + "Scale_" + channel_name + "_8TeV";
                        const std::string name_syst_up = name_syst_prefix + "Up";
                        const std::string name_syst_down = name_syst_prefix + "Down";
                        std::shared_ptr<TH1D> hist_syst_up(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_up->Scale(1.02);
                        hist_syst_up->Write(name_syst_up.c_str());
                        std::shared_ptr<TH1D> hist_syst_down(static_cast<TH1D*>(hist->Clone()));
                        hist_syst_down->Scale(0.98);
                        hist_syst_down->Write(name_syst_down.c_str());
                    }
                }
            }
        }
        outputFile->Close();
    }

    void SubtractBackgroundHistograms(TH1D& histogram, EventCategory eventCategory, EventRegion eventRegion,
                                      const std::string& current_category, std::string& debug_info,
                                      std::string& negative_bins_info)
    {
        static const double correction_factor = 0.0000001;

        std::ostringstream ss_debug;

        ss_debug << "\nSubtracting background for '" << histogram.GetName() << "' in region " << eventRegion
                 << " for Event category '" << eventCategory << "' for data category '" << current_category
                 << "'.\nInitial integral: " << Integral(histogram, true) << ".\n";
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract ) continue;

            ss_debug << "Sample '" << category->name << "': ";
            if(auto other_histogram = GetHistogram(eventCategory, category->name, eventRegion, histogram.GetName())) {
                histogram.Add(other_histogram, -1);
                ss_debug << Integral(*other_histogram, true) << ".\n";
            } else
                ss_debug << "not found.\n";
        }

        const PhysicalValue original_Integral = Integral(histogram, true);
        ss_debug << "Integral after bkg subtraction: " << original_Integral.value << ".\n";
        debug_info = ss_debug.str();
        if (original_Integral.value < 0) {
            std::cout << debug_info << std::endl;
            throw exception("Integral after bkg subtraction is negative for histogram '")
                << histogram.GetName() << "' in event category " << eventCategory << " for event region " << eventRegion
                << ".";
        }

        std::ostringstream ss_negative;

        for (Int_t n = 1; n <= histogram.GetNbinsX(); ++n){
            if (histogram.GetBinContent(n) >= 0) continue;
            const std::string prefix = histogram.GetBinContent(n) + histogram.GetBinError(n) >= 0 ? "WARNING" : "ERROR";

            ss_negative << prefix << ": " << histogram.GetName() << " Bin " << n << ", content = "
                        << histogram.GetBinContent(n) << ", error = " << histogram.GetBinError(n)
                        << ", bin limits=[" << histogram.GetBinLowEdge(n) << "," << histogram.GetBinLowEdge(n+1)
                        << "].\n";
            const double error = correction_factor - histogram.GetBinContent(n);
            const double new_error = std::sqrt(error*error + histogram.GetBinError(n)*histogram.GetBinError(n));
            histogram.SetBinContent(n,correction_factor);
            histogram.SetBinError(n,new_error);
        }
        const PhysicalValue new_Integral = Integral(histogram, true);
        const PhysicalValue correctionSF = original_Integral/new_Integral;
        histogram.Scale(correctionSF.value);
        negative_bins_info = ss_negative.str();
    }

    PhysicalValue CalculateBackgroundIntegral(const std::string& hist_name, EventCategory eventCategory,
                                              EventRegion eventRegion, const std::string& current_category,
                                              bool expect_at_least_one_contribution = false)
    {
        std::string debug_info;
        return CalculateBackgroundIntegral(hist_name, eventCategory, eventRegion, current_category,
                                           expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateBackgroundIntegral(const std::string& hist_name, EventCategory eventCategory,
                                              EventRegion eventRegion, const std::string& current_category,
                                              bool expect_at_least_one_contribution, std::string& debug_info)
    {
        DataCategoryPtrSet bkg_dataCategories;
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category || !category->isCategoryToSubtract ) continue;
            bkg_dataCategories.insert(category);
        }

        return CalculateFullIntegral(eventCategory, eventRegion, hist_name, bkg_dataCategories,
                                     expect_at_least_one_contribution, debug_info);
    }

    PhysicalValue CalculateFullIntegral(analysis::EventCategory eventCategory, analysis::EventRegion eventRegion,
                                        const std::string& hist_name, const DataCategoryPtrSet& dataCategories,
                                        bool expect_at_least_one_contribution = false)
    {
        std::string debug_info;
        return CalculateFullIntegral(eventCategory, eventRegion, hist_name, dataCategories,
                                     expect_at_least_one_contribution, debug_info);
    }


    PhysicalValue CalculateFullIntegral(analysis::EventCategory eventCategory, analysis::EventRegion eventRegion,
                                        const std::string& hist_name, const DataCategoryPtrSet& dataCategories,
                                        bool expect_at_least_one_contribution, std::string& debug_info)
    {
        PhysicalValue integral;
        bool hist_found = false;

        std::ostringstream ss_debug;

        ss_debug << "\nCalculating full integral for '" << hist_name << "' in region " << eventRegion
                 << " for Event category '" << eventCategory << "' considering data categories (" << dataCategories
                 << ").\n";

        for(const auto& dataCategory : dataCategories) {
            ss_debug << "Sample '" << dataCategory->name << "': ";

            if(auto hist = GetHistogram(eventCategory, dataCategory->name, eventRegion, hist_name)) {
                hist_found = true;
                const auto hist_integral = Integral(*hist, true);
                integral += hist_integral;
                ss_debug << hist_integral << ".\n";
            } else
                ss_debug << "not found.\n";
        }

        ss_debug << "Full integral: " << integral << ".\n";
        debug_info = ss_debug.str();

        if(!hist_found && expect_at_least_one_contribution) {
            std::cout << debug_info << std::endl;
            throw exception("No histogram with name '")
                << hist_name << "' was found in the given data category set (" << dataCategories
                << "), eventCategory: " << eventCategory << ", eventRegion: " << eventRegion
                << "to calculate full integral.";
        }

        return integral;
    }

private:

    static const std::vector< std::pair<double, double> >& GetBlindRegion(const std::string& hist_name)
    {
        static const std::vector< std::vector< std::pair<double, double> > > blindingRegions = {
            { { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() } },
            { { 100, 150 } },
            { { 200, 400 } },
            { { 100, 150 }, { 450, 500 }, { 800, 850 }, { 1150, 1200 }, { 1500, 1550 } }
        };

        static const std::map<std::string, size_t> histogramsToBlind = {
            { FlatAnalyzerData::m_sv_Name(), 1 }, { FlatAnalyzerData::m_vis_Name(), 1 },
            { FlatAnalyzerData::m_bb_Name(), 1 }, { FlatAnalyzerData::m_ttbb_Name(), 2 },
            { FlatAnalyzerData::m_ttbb_kinfit_Name(), 2 }, { FlatAnalyzerData::m_bb_slice_Name(), 3 }
        };

        static const std::set<EventSubCategory> sidebandSubCategories = { EventSubCategory::OutsideMassWindow };

        const auto findRegionId = [&]() -> size_t {
            if(histogramsToBlind.count(hist_name))
                return histogramsToBlind.at(hist_name);
            for(EventEnergyScale eventEnergyScale : AllEventEnergyScales) {
                for(EventSubCategory subCategory : AllEventSubCategories) {
                    if(sidebandSubCategories.count(subCategory)) continue;
                    for(const auto& entry : histogramsToBlind) {
                        const auto full_hist_name = FlatAnalyzerData::FullHistogramName(entry.first, subCategory,
                                                                                        eventEnergyScale);
                        if(full_hist_name == hist_name)
                            return entry.second;
                    }
                }
            }
            return 0;
        };

        const size_t regionId = findRegionId();
        if(regionId >= blindingRegions.size())
            throw exception("Bad blinding region index = ") << regionId;
        return blindingRegions.at(regionId);
    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(outputFileName + "_" + name_suffix + ".csv");

        static const std::set<std::string> interesting_histograms = {
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_sv_Name(), EventSubCategory::NoCuts,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_sv_Name(), EventSubCategory::MassWindow,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                                EventSubCategory::KinematicFitConverged,
                                                EventEnergyScale::Central),
            FlatAnalyzerData::FullHistogramName(FlatAnalyzerData::m_ttbb_kinfit_Name(),
                                                EventSubCategory::KinematicFitConvergedWithMassWindow,
                                                EventEnergyScale::Central)
        };

        for(const auto& hist_name : interesting_histograms)
            PrintTables(of, sep, hist_name, false, true);

        of.flush();
        of.close();
    }

    void PrintTables(std::wostream& of, const std::wstring& sep, const std::string& hist_name, bool includeOverflow,
                        bool includeError)
    {
        of << std::wstring(hist_name.begin(), hist_name.end());

        std::wstring table_name_suffix = L"";
        if(includeOverflow && includeError)
            table_name_suffix = L" with overflow and error";
        else if(includeOverflow && !includeError)
            table_name_suffix = L" with overflow";
        else if(!includeOverflow && includeError)
            table_name_suffix = L" with error";
        of << table_name_suffix << sep;

        for (EventCategory eventCategory : EventCategoriesToProcess())
            of << eventCategory << sep;
        of << std::endl;

        for (const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            of << std::wstring(dataCategory->title.begin(), dataCategory->title.end()) << sep;
            for (EventCategory eventCategory : EventCategoriesToProcess()) {
                if( TH1D* histogram = GetSignalHistogram(eventCategory, dataCategory->name, hist_name) ) {
                    const PhysicalValue integral = Integral(*histogram, includeOverflow);
                    of << integral.ToString<wchar_t>(includeError) << sep;
                }
                else
                    of << "not found" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;
    }

    void ProcessCompositDataCategories(EventCategory eventCategory, const std::string& hist_name)
    {
        for (analysis::EventRegion eventRegion : analysis::AllEventRegions){
            for(const DataCategory* composit : dataCategoryCollection.GetCategories(DataCategoryType::Composit)) {
                for(const std::string& sub_name : composit->sub_categories) {
                    const DataCategory& sub_category = dataCategoryCollection.FindCategory(sub_name);
                    //auto sub_hist = GetSignalHistogram(eventCategory, sub_category.name, hist_name);
                    auto sub_hist = GetHistogram(eventCategory,sub_category.name,eventRegion,hist_name);
                    if(!sub_hist) continue;

//                    if(auto composit_hist = GetSignalHistogram(eventCategory, composit->name, hist_name))
                    if(auto composit_hist = GetHistogram(eventCategory, composit->name, eventRegion, hist_name))
                        composit_hist->Add(sub_hist);
                    else
                        //CloneSignalHistogram(eventCategory, composit->name, *sub_hist);
                        CloneHistogram(eventCategory, composit->name,eventRegion, *sub_hist);
                }
            }
        }
    }

protected:
    std::string inputPath;
    std::string outputFileName;
    DataCategoryCollection dataCategoryCollection;
    FullAnaData fullAnaData;
};

} // namespace analysis
