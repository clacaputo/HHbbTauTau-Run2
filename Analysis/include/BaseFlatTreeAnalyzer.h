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


#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatEventInfo.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/AnalysisTypes.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"
//#include "KinFit.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

static const std::vector<double> mass_bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140., 150,
                                               160, 170, 180, 190, 200, 225, 250, 275, 300, 325, 350 };

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:
    TH1D_ENTRY(pt_1, 50, 0, 200)
    TH1D_ENTRY(eta_1, 60, -3, 3)
    TH1D_ENTRY(pt_2, 50, 0, 200)
    TH1D_ENTRY(eta_2, 60, -3, 3)
    TH1D_ENTRY(pt_b1, 50, 0, 200)
    TH1D_ENTRY(eta_b1, 60, -3, 3)
    TH1D_ENTRY(pt_b2, 50, 0, 200)
    TH1D_ENTRY(eta_b2, 60, -3, 3)
    TH1D_ENTRY(pt_H_tt, 35, 0, 350)
    TH1D_ENTRY(pt_H_bb, 35, 0, 350)
    TH1D_ENTRY(pt_H_hh, 35, 0, 350)
    TH1D_ENTRY_CUSTOM(m_sv, mass_bins)
    TH1D_ENTRY_CUSTOM(m_sv_up, mass_bins)
    TH1D_ENTRY_CUSTOM(m_sv_down, mass_bins)
    TH1D_ENTRY_CUSTOM(m_vis, mass_bins)
    TH1D_ENTRY(m_bb, 30, 0, 600)
    TH1D_ENTRY(m_ttbb, 50, 0, 1000)
    TH1D_ENTRY(m_ttbb_kinfit, 50, 0, 1000)
    TH1D_ENTRY(m_ttbb_kinfit_up, 50, 0, 1000)
    TH1D_ENTRY(m_ttbb_kinfit_down, 50, 0, 1000)
    TH1D_ENTRY(DeltaPhi_tt, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_bb, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_bb_MET, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_tt_MET, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_hh, 80, -4, 4)
    TH1D_ENTRY(DeltaR_tt, 60, 0, 6)
    TH1D_ENTRY(DeltaR_bb, 60, 0, 6)
    TH1D_ENTRY(DeltaR_hh, 60, 0, 6)
    TH1D_ENTRY(MVA_BDT, 40, -1, 1)
    TH1D_ENTRY(MVA_BDTD, 40, -1, 1)
    TH1D_ENTRY(MVA_BDTMitFisher, 40, -1, 1)
    TH1D_ENTRY(mt_1, 50, 0, 50)
    TH1D_ENTRY(mt_2, 50, 0, 300)
    TH1D_ENTRY(pt_H_tt_MET, 35, 0, 350)
};

class BaseFlatTreeAnalyzer {
public:
    typedef std::map<EventRegion, FlatAnalyzerData> RegionAnaData;

    typedef std::map<std::string, RegionAnaData> AnaDataForEventCategory;
    typedef std::map<EventCategory, AnaDataForEventCategory> FullAnaData;

    static const std::string ReferenceHistogramName() { return "m_sv"; }

    BaseFlatTreeAnalyzer(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& _outputFileName, Channel channel_id,
                         const std::string& signal_list)
        : inputPath(_inputPath), outputFileName(_outputFileName),
          dataCategoryCollection(source_cfg, signal_list, channel_id)
    {
        TH1::SetDefaultSumw2();

        histograms = HistogramDescriptor::ReadFromFile(hist_cfg);
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

        std::cout << "Calculating embedded scale factor... " << std::endl;
        //const double embeddedSF = CalculateEmbeddedScaleFactor(ReferenceHistogramName());
        const double embeddedSF = 1;

        std::cout << "Estimating QCD, Wjets and composit data categories... " << std::endl;
        for (auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForEventCategory& anaData = fullAnaDataEntry.second;

            CreateHistogramForZTT(anaData, ReferenceHistogramName(), embeddedSF);

            const auto wjets_scale_factors = CalculateWjetsScaleFactors(eventCategory, ReferenceHistogramName());
            std::cout << eventCategory << " OS_HighMt_data / OS_HighMt_mc = "
                      << wjets_scale_factors.first << "\n"
                      << eventCategory << " SS_HighMt_data / SS_HighMt_mc = "
                      << wjets_scale_factors.second << std::endl;
            EstimateWjets(eventCategory, ReferenceHistogramName(), wjets_scale_factors);

            const auto qcd_scale_factor = CalculateQCDScaleFactor(eventCategory, ReferenceHistogramName());
            std::cout << eventCategory << " OS_NotIso / SS_NotIso = " << qcd_scale_factor << std::endl;

            for (const auto& hist : histograms) {
                if(hist.name != ReferenceHistogramName()) {
                    CreateHistogramForZTT(anaData, hist.name, embeddedSF);
                    EstimateWjets(eventCategory, hist.name, wjets_scale_factors);
                }
                EstimateQCD(eventCategory, hist.name, qcd_scale_factor);
                ProcessCompositDataCategories(eventCategory, hist.name);
            }
        }

        std::cout << "Saving tables... " << std::endl;
        PrintTables("comma", L",");
        PrintTables("semicolon", L";");

        std::cout << "Saving datacards... " << std::endl;
        ProduceFileForLimitsCalculation("m_sv", "m_sv_up", "m_sv_down", false);
        ProduceFileForLimitsCalculation("m_ttbb_kinfit", "m_ttbb_kinfit_up", "m_ttbb_kinfit_down", false);

        std::cout << "Printing stacked plots... " << std::endl;
        PrintStackedPlots(false);
        PrintStackedPlots(true);
    }

protected:
    virtual Channel ChannelId() const = 0;
    virtual EventRegion DetermineEventRegion(const ntuple::Flat& event) = 0;
    virtual bool PassMvaCut(const FlatEventInfo& eventInfo, EventCategory eventCategory) { return true; }

    virtual PhysicalValue CalculateQCDScaleFactor(EventCategory eventCategory, const std::string& hist_name)
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_OSnotIso_data = GetHistogram(eventCategory, data.name, EventRegion::OS_NotIsolated, hist_name);
        auto hist_SSnotIso_data = GetHistogram(eventCategory, data.name, EventRegion::SS_NotIsolated, hist_name);
        if(!hist_OSnotIso_data || !hist_SSnotIso_data)
            throw analysis::exception("Unable to find histograms for QCD scale factor estimation");

        TH1D& hist_OSnotIso = CloneHistogram(eventCategory, qcd.name, EventRegion::OS_NotIsolated, *hist_OSnotIso_data);
        TH1D& hist_SSnotIso = CloneHistogram(eventCategory, qcd.name, EventRegion::SS_NotIsolated, *hist_SSnotIso_data);

        SubtractBackgroundHistograms(hist_OSnotIso, eventCategory, EventRegion::OS_NotIsolated, qcd.name, true);
        SubtractBackgroundHistograms(hist_SSnotIso, eventCategory, EventRegion::SS_NotIsolated, qcd.name, true);

        const PhysicalValue n_OSnotIso = Integral(hist_OSnotIso, false);
        const PhysicalValue n_SSnotIso = Integral(hist_SSnotIso, false);
        return n_OSnotIso / n_SSnotIso;
    }

    virtual void EstimateQCD(EventCategory eventCategory, const std::string& hist_name,
                             const PhysicalValue& scale_factor)
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& qcd = dataCategoryCollection.GetUniqueCategory(DataCategoryType::QCD);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_SSIso_data = GetHistogram(eventCategory, data.name, EventRegion::SS_Isolated, hist_name);
        if(!hist_SSIso_data) return;

        TH1D& histogram = CloneHistogram(eventCategory, qcd.name, EventRegion::OS_Isolated, *hist_SSIso_data);
        SubtractBackgroundHistograms(histogram, eventCategory, EventRegion::SS_Isolated, qcd.name);
        histogram.Scale(scale_factor.value);
    }

    virtual PhysicalValuePair CalculateWjetsScaleFactors(EventCategory eventCategory, const std::string& hist_name)
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategory& wjets_mc = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets_MC);
        const analysis::DataCategory& data = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Data);

        auto hist_OS_HighMt_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::OS_HighMt, hist_name);
        auto hist_SS_HighMt_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::SS_HighMt, hist_name);
        auto hist_OS_HighMt_data = GetHistogram(eventCategory, data.name, EventRegion::OS_HighMt, hist_name);
        auto hist_SS_HighMt_data = GetHistogram(eventCategory, data.name, EventRegion::SS_HighMt, hist_name);
        if(!hist_OS_HighMt_mc || !hist_SS_HighMt_mc || !hist_OS_HighMt_data || !hist_SS_HighMt_data)
            throw analysis::exception("Unable to find histograms for Wjet scale factors estimation");

        TH1D& hist_OS_HighMt = CloneHistogram(eventCategory, wjets.name, EventRegion::OS_HighMt, *hist_OS_HighMt_data);
        TH1D& hist_SS_HighMt = CloneHistogram(eventCategory, wjets.name, EventRegion::SS_HighMt, *hist_SS_HighMt_data);

        SubtractBackgroundHistograms(hist_OS_HighMt, eventCategory, EventRegion::OS_HighMt, wjets.name, true);
        SubtractBackgroundHistograms(hist_SS_HighMt, eventCategory, EventRegion::SS_HighMt, wjets.name, true);

        const PhysicalValue n_OS_HighMt = Integral(hist_OS_HighMt, false);
        const PhysicalValue n_OS_HighMt_mc = Integral(*hist_OS_HighMt_mc, false);
        const PhysicalValue n_SS_HighMt = Integral(hist_SS_HighMt, false);
        const PhysicalValue n_SS_HighMt_mc = Integral(*hist_SS_HighMt_mc, false);

        const PhysicalValue OS_ratio = n_OS_HighMt / n_OS_HighMt_mc;
        const PhysicalValue SS_ratio = n_SS_HighMt / n_SS_HighMt_mc;
        return PhysicalValuePair(OS_ratio, SS_ratio);
    }

    virtual void EstimateWjets(EventCategory eventCategory, const std::string& hist_name,
                               const PhysicalValuePair& scale_factors)
    {
        using analysis::EventRegion;
        using analysis::DataCategoryType;

        const analysis::DataCategory& wjets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets);
        const analysis::DataCategory& wjets_mc = dataCategoryCollection.GetUniqueCategory(DataCategoryType::WJets_MC);

        if(auto hist_OS_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::OS_Isolated, hist_name)) {
            TH1D& hist_OS = CloneHistogram(eventCategory, wjets.name, EventRegion::OS_Isolated, *hist_OS_mc);
            hist_OS.Scale(scale_factors.first.value);
        }

        if(auto hist_SS_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::SS_Isolated, hist_name)) {
            TH1D& hist_SS = CloneHistogram(eventCategory, wjets.name, EventRegion::SS_Isolated, *hist_SS_mc);
            hist_SS.Scale(scale_factors.second.value);
        }

        if(auto hist_OS_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::OS_NotIsolated, hist_name)) {
            TH1D& hist_OS = CloneHistogram(eventCategory, wjets.name, EventRegion::OS_NotIsolated, *hist_OS_mc);
            hist_OS.Scale(scale_factors.first.value);
        }

        if(auto hist_SS_mc = GetHistogram(eventCategory, wjets_mc.name, EventRegion::SS_NotIsolated, hist_name)) {
            TH1D& hist_SS = CloneHistogram(eventCategory, wjets.name, EventRegion::SS_NotIsolated, *hist_SS_mc);
            hist_SS.Scale(scale_factors.second.value);
        }
    }

    const std::string& ChannelName() const { return detail::ChannelNameMap.at(ChannelId()); }

    root_ext::SmartHistogram<TH1D>* GetHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                                                 EventRegion eventRegion, const std::string& histogramName)
    {
        return fullAnaData[eventCategory][dataCategoryName][eventRegion].GetPtr<TH1D>(histogramName);
    }

    root_ext::SmartHistogram<TH1D>* GetSignalHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                                                       const std::string& histogramName)
    {
        return GetHistogram(eventCategory, dataCategoryName, EventRegion::OS_Isolated, histogramName);
    }

    TH1D& CloneHistogram(EventCategory eventCategory, const std::string& dataCategoryName, EventRegion eventRegion,
                         const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return fullAnaData[eventCategory][dataCategoryName][eventRegion].Clone(originalHistogram);
    }

    TH1D& CloneSignalHistogram(EventCategory eventCategory, const std::string& dataCategoryName,
                         const root_ext::SmartHistogram<TH1D>& originalHistogram)
    {
        return CloneHistogram(eventCategory, dataCategoryName, EventRegion::OS_Isolated, originalHistogram);
    }


    void ProcessDataSource(const DataCategory& dataCategory, std::shared_ptr<ntuple::FlatTree> tree, double scale_factor)
    {
        static const bool applyMVAcut = false;

        const analysis::DataCategory& DYJets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets);

        for(Long64_t current_entry = 0; current_entry < tree->GetEntries(); ++current_entry) {
            tree->GetEntry(current_entry);
            const ntuple::Flat& event = tree->data;

            const EventRegion eventRegion = DetermineEventRegion(event);
            if(eventRegion == EventRegion::Unknown) continue;

            const EventCategoryVector eventCategories = DetermineEventCategories(event.csv_Bjets,
                                                                                 cuts::Htautau_Summer13::btag::CSVM,
                                                                                 cuts::Htautau_Summer13::btag::CSVT);
            FlatEventInfo eventInfo(event, FlatEventInfo::BjetPair(0, 1));

            const double weight = dataCategory.IsData() ? 1 : event.weight * scale_factor;

            for(auto eventCategory : eventCategories) {
                UpdateMvaInfo(eventInfo, eventCategory, false, false, false);
                if(applyMVAcut && !PassMvaCut(eventInfo, eventCategory)) continue;
                if (dataCategory.name == DYJets.name)
                    FillDYjetHistograms(eventInfo, eventCategory, eventRegion, weight);

                FillHistograms(fullAnaData[eventCategory][dataCategory.name][eventRegion], eventInfo, weight);
            }
        }
    }

    void FillDYjetHistograms(const FlatEventInfo& eventInfo, EventCategory eventCategory, EventRegion eventRegion,
                             double weight)
    {
        const analysis::DataCategory& ZL = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZL);
        const analysis::DataCategory& ZJ = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZJ);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        const std::map<ntuple::EventType, std::string> type_category_map = {
            { ntuple::EventType::ZL, ZL.name }, { ntuple::EventType::ZJ, ZJ.name },
            { ntuple::EventType::ZTT, ZTT_MC.name }, { ntuple::EventType::ZTT_no_match, ZTT_MC.name }
        };

        if(type_category_map.count(eventInfo.eventType)) {
            const std::string& name = type_category_map.at(eventInfo.eventType);
            FillHistograms(fullAnaData[eventCategory][name][eventRegion], eventInfo, weight);
        }
    }

    double CalculateEmbeddedScaleFactor(const std::string& hist_name)
    {
        const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        TH1D* hist_embedded = GetSignalHistogram(EventCategory::Inclusive, embedded.name, hist_name);
        TH1D* hist_ztautau = GetSignalHistogram(EventCategory::Inclusive, ZTT_MC.name, hist_name);
        if(!hist_embedded || !hist_ztautau )
            throw std::runtime_error("embedded or ztt hist not found");

        const double n_ztautau = hist_ztautau->Integral(0, hist_ztautau->GetNbinsX() + 1);
        const double n_embedded = hist_embedded->Integral(0, hist_embedded->GetNbinsX() + 1);
        return n_ztautau / n_embedded;
    }

    void CreateHistogramForZTT(AnaDataForEventCategory& anaData, const std::string& hist_name, double scale_factor)
    {
        //const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);
        const analysis::DataCategory& ZTT = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT);

        for(auto& map_entry : anaData[embedded.name]) {
            if(auto embedded_hist = map_entry.second.GetPtr<TH1D>(hist_name)) {
                TH1D& ztt_hist = anaData[ZTT.name][map_entry.first].Clone(*embedded_hist);
                ztt_hist.Scale(scale_factor);
            }
        }
    }

    void FillHistograms(FlatAnalyzerData& anaData, const FlatEventInfo& eventInfo, double weight)
    {
        const ntuple::Flat& event = *eventInfo.event;
        anaData.m_sv().Fill(event.m_sv_vegas, weight);
        anaData.m_sv_up().Fill(event.m_sv_up_vegas, weight);
        anaData.m_sv_down().Fill(event.m_sv_down_vegas, weight);
        anaData.pt_1().Fill(event.pt_1, weight);
        anaData.eta_1().Fill(event.eta_1, weight);
        anaData.pt_2().Fill(event.pt_2, weight);
        anaData.eta_2().Fill(event.eta_2, weight);
        anaData.DeltaPhi_tt().Fill(eventInfo.lepton_momentums.at(0).DeltaPhi(eventInfo.lepton_momentums.at(1)), weight);
        anaData.DeltaR_tt().Fill(eventInfo.lepton_momentums.at(0).DeltaR(eventInfo.lepton_momentums.at(1)), weight);
        anaData.pt_H_tt().Fill(eventInfo.Htt.Pt(),weight);
        anaData.m_vis().Fill(eventInfo.Htt.M(),weight);
        anaData.pt_H_tt_MET().Fill(eventInfo.Htt_MET.Pt(), weight);
        anaData.DeltaPhi_tt_MET().Fill(eventInfo.Htt.DeltaPhi(eventInfo.MET), weight);
        anaData.mt_1().Fill(event.mt_1, weight);
        anaData.mt_2().Fill(event.mt_2, weight);
        if(eventInfo.has_bjet_pair) {
            anaData.pt_b1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Pt(), weight);
            anaData.eta_b1().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).Eta(), weight);
            anaData.pt_b2().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Pt(), weight);
            anaData.eta_b2().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second).Eta(), weight);
            anaData.DeltaPhi_bb().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaPhi(
                                           eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second)), weight);
            anaData.DeltaR_bb().Fill(eventInfo.bjet_momentums.at(eventInfo.selected_bjets.first).DeltaR(
                                         eventInfo.bjet_momentums.at(eventInfo.selected_bjets.second)), weight);
            anaData.pt_H_bb().Fill(eventInfo.Hbb.Pt(),weight);
            anaData.m_bb().Fill(eventInfo.Hbb.M(), weight);
            anaData.DeltaPhi_bb_MET().Fill(eventInfo.Hbb.DeltaPhi(eventInfo.MET), weight);
            anaData.DeltaPhi_hh().Fill(eventInfo.Htt.DeltaPhi(eventInfo.Hbb), weight);
            anaData.DeltaR_hh().Fill(eventInfo.Htt.DeltaR(eventInfo.Hbb), weight);
            anaData.m_ttbb().Fill(eventInfo.resonance.M(), weight);
            anaData.pt_H_hh().Fill(eventInfo.resonance.Pt(), weight);
            const double m_ttbb_kinFit = event.kinfit_bb_tt_mass.at(eventInfo.kinfit_data_index);
            anaData.m_ttbb_kinfit().Fill(m_ttbb_kinFit,weight);
            anaData.m_ttbb_kinfit_up().Fill(1.04*m_ttbb_kinFit,weight);
            anaData.m_ttbb_kinfit_down().Fill(0.96*m_ttbb_kinFit,weight);

            anaData.MVA_BDT().Fill(eventInfo.mva_BDT, weight);
            anaData.MVA_BDTD().Fill(eventInfo.mva_BDTD, weight);
            anaData.MVA_BDTMitFisher().Fill(eventInfo.mva_BDTMitFisher, weight);
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

    void PrintStackedPlots(bool isBlind)
    {
        const std::string blindCondition = isBlind ? "_blind" : "_noBlind";
        root_ext::PdfPrinter printer(outputFileName + blindCondition + ".pdf");

        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory eventCategory = fullAnaDataEntry.first;

            for (const HistogramDescriptor& hist : histograms) {
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist.title;
                StackedPlotDescriptor stackDescriptor(hist, ss_title.str());

                for(const DataCategory* category : dataCategoryCollection.GetAllCategories()) {
                    if(!category->draw) continue;

                    TH1D* histogram = GetSignalHistogram(eventCategory, category->name, hist.name);
                    if(!histogram) continue;

                    if(category->IsSignal())
                        stackDescriptor.AddSignalHistogram(histogram, category->title, category->color, category->draw_sf);
                    else if(category->IsBackground())
                        stackDescriptor.AddBackgroundHistogram(histogram, category->title, category->color);
                    else if(category->IsData())
                        stackDescriptor.AddDataHistogram(histogram, category->title, isBlind, GetBlindRegion(hist.name));
                }

                printer.PrintStack(stackDescriptor);
            }
        }
    }

    void ProduceFileForLimitsCalculation(const std::string& hist_name, const std::string& hist_name_up,
                                         const std::string& hist_name_down, bool include_one_jet_categories)
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::OneJet_ZeroBtag, "1jet0tag" },
            { EventCategory::OneJet_OneBtag, "1jet1tag" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);

        std::shared_ptr<TFile> outputFile(new TFile((outputFileName + hist_name + ".root").c_str(), "RECREATE"));
        outputFile->cd();
        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            if(!categoryToDirectoryNameSuffix.count(eventCategory)
                    || (!include_one_jet_categories && OneJetEventCategories.count(eventCategory))) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '") << dataCategory->name << "'.";
                TH1D* hist_orig = GetSignalHistogram(eventCategory, dataCategory->name, hist_name);
                if(!hist_orig)
                    throw exception("Datacard histogram '") << hist_name << "' not found for data category '"
                                                            << dataCategory->name << "'.";
                std::shared_ptr<TH1D> hist(static_cast<TH1D*>(hist_orig->Clone()));
                hist->Scale(dataCategory->limits_sf);
                hist->Write(dataCategory->datacard.c_str());
                const std::string namePrefix = dataCategory->datacard + "_CMS_scale_t_" + channel_name + "_8TeV";
                const std::string nameDown = namePrefix + "Down";
                const std::string nameUp = namePrefix + "Up";

                TH1D* hist_up_orig = GetSignalHistogram(eventCategory, dataCategory->name, hist_name_up);
                if(!hist_up_orig)
                    throw exception("Datacard histogram '") << hist_name_up << "' not found for data category '"
                                                            << dataCategory->name << "'.";
                std::shared_ptr<TH1D> hist_up(static_cast<TH1D*>(hist_up_orig->Clone()));
                hist_up->Scale(dataCategory->limits_sf);
                hist_up->Write(nameUp.c_str());
                TH1D* hist_down_orig = GetSignalHistogram(eventCategory, dataCategory->name, hist_name_down);
                if(!hist_down_orig)
                    throw exception("Datacard histogram '") << hist_name_down << "' not found for data category '"
                                                            << dataCategory->name << "'.";
                std::shared_ptr<TH1D> hist_down(static_cast<TH1D*>(hist_down_orig->Clone()));
                hist_down->Scale(dataCategory->limits_sf);
                hist_down->Write(nameDown.c_str());
            }
        }
        outputFile->Close();
    }

    void SubtractBackgroundHistograms(TH1D& histogram, EventCategory eventCategory, EventRegion eventRegion,
                                      const std::string& current_category, bool verbose = false)
    {
        if(verbose)
            std::cout << "\nSubtracting background for '" << histogram.GetName() << "' in region " << eventRegion
                      << " for data category '" << current_category << "'.\n"
                      << "Initial integral: " << Integral(histogram, false) << ".\n";
        for (auto category : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            if(category->IsComposit() || category->name == current_category) continue;

            if(verbose)
                std::cout << "Sample '" << category->name << "': ";
            if(auto other_histogram = GetHistogram(eventCategory, category->name, eventRegion, histogram.GetName())) {
                histogram.Add(other_histogram, -1);
                if(verbose)
                    std::cout << Integral(*other_histogram, false) << ".\n";
            } else if(verbose)
                std::cout << "not found.\n";
        }
        if(verbose)
            std::cout << "Integral after bkg subtraction: " << Integral(histogram, false) << ".\n" << std::endl;
    }

private:

    static const std::pair<double, double>& GetBlindRegion(const std::string& hist_name)
    {
        static const std::vector< std::pair<double, double> > blindingRegions = {
            { std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest() }, { 100, 150 }, { 250, 350 }
        };
        static const std::map<std::string, size_t> histogramsToBlind = {
            { "m_sv", 1 }, { "m_sv_up", 1 }, { "m_sv_down", 1 }, { "m_vis", 1 }, { "m_bb", 1 },
            { "m_ttbb", 2 }, { "m_ttbb_kinfit", 2 }, { "m_ttbb_kinfit_up", 2 }, { "m_ttbb_kinfit_down", 2 }
        };

        if(!histogramsToBlind.count(hist_name)) return blindingRegions.at(0);
        const size_t regionId = histogramsToBlind.at(hist_name);
        if(regionId >= blindingRegions.size())
            throw analysis::exception("Bad blinding region index = ") << regionId;
        return blindingRegions.at(regionId);
    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(outputFileName + "_" + name_suffix + ".csv");

        for(const HistogramDescriptor& hist : histograms) {
            if (hist.name != "m_sv") continue;
            PrintTables(of, sep, hist, false,false);
            PrintTables(of, sep, hist, true, false);
            PrintTables(of, sep, hist, false,true);
            PrintTables(of, sep, hist, true, true);
        }
        of.flush();
        of.close();
    }

    void PrintTables(std::wostream& of, const std::wstring& sep, const HistogramDescriptor& hist, bool includeOverflow,
                        bool includeError)
    {
        of << std::wstring(hist.title.begin(), hist.title.end());

        std::wstring table_name_suffix = L"";
        if(includeOverflow && includeError)
            table_name_suffix = L" with overflow and error";
        else if(includeOverflow && !includeError)
            table_name_suffix = L" with overflow";
        else if(!includeOverflow && includeError)
            table_name_suffix = L" with error";
        of << table_name_suffix << sep;

        for (const auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            of << eventCategory << sep;
        }
        of << std::endl;

        for (const DataCategory* dataCategory : dataCategoryCollection.GetAllCategories()) {
            of << std::wstring(dataCategory->title.begin(), dataCategory->title.end()) << sep;
            for (auto& fullAnaDataEntry : fullAnaData) {
                const EventCategory& eventCategory = fullAnaDataEntry.first;
                if( TH1D* histogram = GetSignalHistogram(eventCategory, dataCategory->name, hist.name) ) {
                    const PhysicalValue integral = Integral(*histogram, includeOverflow);
                    if(includeError) of << integral;
                    else of << integral.value;
                    of << sep;
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
        for(const DataCategory* composit : dataCategoryCollection.GetCategories(DataCategoryType::Composit)) {
            for(const std::string& sub_name : composit->sub_categories) {
                const DataCategory& sub_category = dataCategoryCollection.FindCategory(sub_name);
                auto sub_hist = GetSignalHistogram(eventCategory, sub_category.name, hist_name);
                if(!sub_hist) continue;

                if(auto composit_hist = GetSignalHistogram(eventCategory, composit->name, hist_name))
                    composit_hist->Add(sub_hist);
                else
                    CloneSignalHistogram(eventCategory, composit->name, *sub_hist);
            }
        }
    }

protected:
    std::string inputPath;
    std::string outputFileName;
    DataCategoryCollection dataCategoryCollection;
    std::vector<HistogramDescriptor> histograms;
    FullAnaData fullAnaData;
};

} // namespace analysis
