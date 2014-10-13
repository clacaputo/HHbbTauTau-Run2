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
    typedef std::map<EventType_QCD, FlatAnalyzerData> AnaDataQCD;
    typedef std::map<EventType_Wjets, FlatAnalyzerData> AnaDataWjets;

    struct AnaDataEventTypes {
        AnaDataQCD QCD;
        AnaDataWjets Wjets;

        FlatAnalyzerData& Signal() { return QCD[EventType_QCD::OS_Isolated]; }
    };

    typedef std::map<std::string, AnaDataEventTypes> AnaDataForDataCategory;
    typedef std::map<EventCategory, AnaDataForDataCategory> FullAnaData;

    BaseFlatTreeAnalyzer(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& _outputFileName, Channel channel_id,
                         const std::string& signal_list, bool _WjetsData = false)
        : inputPath(_inputPath), outputFileName(_outputFileName),
          dataCategoryCollection(source_cfg, signal_list, channel_id), WjetsData(_WjetsData)
    {
        TH1::SetDefaultSumw2();

        histograms = HistogramDescriptor::ReadFromFile(hist_cfg);
    }

    void Run()
    {
        std::cout << "Processing input source files... " << std::endl;
        for(const DataSource* source : dataCategoryCollection.GetSources()) {
            std::cout << *source << std::endl;
            const std::string fullFileName = inputPath + "/" + source->file_name;
            std::shared_ptr<TFile> file(new TFile(fullFileName.c_str(), "READ"));
            if(file->IsZombie())
                throw exception("Input file '") << source->file_name << "' not found.";
            std::shared_ptr<ntuple::FlatTree> tree(new ntuple::FlatTree(*file, "flatTree"));
            ProcessDataSource(*source, tree);
        }

        std::cout << "Calculating embedded scale factor... " << std::endl;
        const double embeddedSF = CalculateEmbeddedScaleFactor();

        std::cout << "Estimating QCD, Wjets and bkg sum... " << std::endl;
        for (auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const auto& hist : histograms) {
                CreateHistogramForZTT(anaData,hist,embeddedSF);
                if (WjetsData) EstimateWjets(eventCategory, anaData, hist);
                EstimateQCD(eventCategory, anaData, hist);
                EstimateSumBkg(eventCategory, anaData, hist);
            }
        }

        std::cout << "Saving tables and printing stacked plots... " << std::endl;
        PrintTables("comma", L",");
        PrintTables("semicolon", L";");
        ProduceFileForLimitsCalculation("m_sv","m_sv_up","m_sv_down");
        ProduceFileForLimitsCalculation("m_ttbb_kinfit","m_ttbb_kinfit_up","m_ttbb_kinfit_down");
        std::cout << "plots for limits done" << std::endl;
        PrintStackedPlots(false);
        PrintStackedPlots(true);
    }

protected:
    void ProcessDataSource(const DataSource& dataSource, std::shared_ptr<ntuple::FlatTree> tree)
    {
        static const bool applyMVAcut = false;

        const analysis::DataCategory& DYJets = dataCategoryCollection.GetUniqueCategory(DataCategoryType::DYJets);

        for(size_t current_entry = 0; current_entry < tree->GetEntries(); ++current_entry) {
            tree->GetEntry(current_entry);
            const ntuple::Flat& event = tree->data;

            const EventType_QCD eventTypeQCD = DetermineEventTypeForQCD(event);
            const EventType_Wjets eventTypeWjets = DetermineEventTypeForWjets(event);
            if(eventTypeQCD == EventType_QCD::Unknown && eventTypeWjets == EventType_Wjets::Unknown)
                continue;

            const EventCategoryVector eventCategories = DetermineEventCategories(event);
            FlatEventInfo eventInfo(event, FlatEventInfo::BjetPair(0, 1));

            for(const DataCategory* dataCategory : dataSource.data_categories) {
                const double weight = dataCategory->IsData() ? 1 : event.weight * dataSource.scale_factor(dataCategory);

                for(auto eventCategory : eventCategories) {
                    UpdateMvaInfo(eventInfo, eventCategory, false, false, false);
                    if(applyMVAcut && !PassMvaCut(eventInfo, eventCategory)) continue;
                    if (dataCategory->name == DYJets.name)
                        FillDYjetHistograms(eventInfo, eventCategory, eventTypeQCD, eventTypeWjets, weight);

                    FillHistograms(fullAnaData[eventCategory][dataCategory->name].QCD[eventTypeQCD], eventInfo, weight);
                    FillHistograms(fullAnaData[eventCategory][dataCategory->name].Wjets[eventTypeWjets], eventInfo, weight);
                }
            }
        }
    }

    void FillDYjetHistograms(const FlatEventInfo& eventInfo, EventCategory eventCategory, EventType_QCD eventTypeQCD,
                             EventType_Wjets eventTypeWjets, double weight)
    {
        const analysis::DataCategory& ZL = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZL);
        const analysis::DataCategory& ZJ = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZJ);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        if (eventInfo.eventType == ntuple::EventType::ZL) {
            FillHistograms(fullAnaData[eventCategory][ZL.name].QCD[eventTypeQCD], eventInfo, weight);
            FillHistograms(fullAnaData[eventCategory][ZL.name].Wjets[eventTypeWjets], eventInfo, weight);
        } else if (eventInfo.eventType == ntuple::EventType::ZJ) {
            FillHistograms(fullAnaData[eventCategory][ZJ.name].QCD[eventTypeQCD], eventInfo, weight);
            FillHistograms(fullAnaData[eventCategory][ZJ.name].Wjets[eventTypeWjets], eventInfo, weight);
        } else if (eventInfo.eventType == ntuple::EventType::ZTT
                   && std::abs(eventInfo.event->pdgId_1_MC) != particles::NONEXISTENT.RawCode()
                   && std::abs(eventInfo.event->pdgId_2_MC) != particles::NONEXISTENT.RawCode()) {
            FillHistograms(fullAnaData[eventCategory][ZTT_MC.name].QCD[eventTypeQCD], eventInfo, weight);
            FillHistograms(fullAnaData[eventCategory][ZTT_MC.name].Wjets[eventTypeWjets], eventInfo, weight);
        }
    }

    virtual Channel ChannelId() = 0;

    virtual EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) = 0;
    virtual EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) = 0;
    virtual bool PassMvaCut(const FlatEventInfo& eventInfo, EventCategory eventCategory) = 0;

    const std::string& ChannelName() const { return detail::ChannelNameMap.at(ChannelId()); }

    virtual EventCategoryVector DetermineEventCategories(const ntuple::Flat& event)
    {
        EventCategoryVector categories;
        categories.push_back(EventCategory::Inclusive);

        if(event.csv_Bjets.size() == 1) {
            const EventCategory category = event.csv_Bjets.at(0) > cuts::Htautau_Summer13::btag::CSVT
                    ? EventCategory::OneJet_OneBtag : EventCategory::OneJet_ZeroBtag;
            categories.push_back(category);
        } else if(event.csv_Bjets.size() >= 2) {
            EventCategory category;
            if (event.csv_Bjets.at(0) <= cuts::Htautau_Summer13::btag::CSVM )
                category = EventCategory::TwoJets_ZeroBtag;
            else if ( event.csv_Bjets.at(1) <= cuts::Htautau_Summer13::btag::CSVM )
                category = EventCategory::TwoJets_OneBtag;
            else
                category = EventCategory::TwoJets_TwoBtag;
            categories.push_back(category);
        }

        return categories;
    }

    double CalculateEmbeddedScaleFactor()
    {
        static const std::string hist_name = "m_sv";

        const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const analysis::DataCategory& ZTT_MC = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT_MC);

        TH1D* hist_embedded =
                fullAnaData[EventCategory::Inclusive][embedded.name].Signal().GetPtr<TH1D>(hist_name);
        TH1D* hist_ztautau =
                fullAnaData[EventCategory::Inclusive][ZTT_MC.name].Signal().GetPtr<TH1D>(hist_name);

        if(!hist_embedded || !hist_ztautau )
            throw std::runtime_error("embedded or ztt hist not found");
        const double n_ztautau = hist_ztautau->Integral(0, hist_ztautau->GetNbinsX() + 1);
        const double n_embedded = hist_embedded->Integral(0, hist_embedded->GetNbinsX() + 1);
        return n_ztautau / n_embedded;
    }

    void CreateHistogramForZTT(AnaDataForDataCategory& anaData, const HistogramDescriptor& hist, double scale_factor)
    {
        const analysis::DataCategory& embedded = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Embedded);
        const analysis::DataCategory& ZTT = dataCategoryCollection.GetUniqueCategory(DataCategoryType::ZTT);

        for(auto& map_entry : anaData[embedded.name].QCD) {
            if(auto embedded_hist = map_entry.second.GetPtr<TH1D>(hist.name)) {
                TH1D& ztt_hist = anaData[ZTT.name].QCD[map_entry.first].Clone(*embedded_hist);
                ztt_hist.Scale(scale_factor);
            }
        }

        for(auto& map_entry : anaData[embedded.name].Wjets) {
            if(auto embedded_hist = map_entry.second.GetPtr<TH1D>(hist.name)) {
                TH1D& ztt_hist = anaData[ZTT.name].Wjets[map_entry.first].Clone(*embedded_hist);
                ztt_hist.Scale(scale_factor);
            }
        }
    }

    virtual void EstimateQCD(EventCategory eventCategory, AnaDataForDataCategory& anaData,
                             const HistogramDescriptor& hist) = 0;
    virtual void EstimateWjets(EventCategory eventCategory, AnaDataForDataCategory& anaData,
                               const HistogramDescriptor& hist) = 0;

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
        const std::string category_name = eventCategoryMapName.at(eventCategory);

        auto getMVA = [&](bool calc_MVA, MVA_Selections::MvaMethod method) -> double {
            if(calc_MVA) {
                auto mvaReader = MVA_Selections::MvaReader::Get(ChannelName(), category_name, method);
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
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;

            for (const HistogramDescriptor& hist : histograms) {
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist.title;
                StackedPlotDescriptor stackDescriptor(hist, ss_title.str());

                for(const DataCategory* category : dataCategoryCollection.GetAllCategories()) {
                    if(!category->draw) continue;

                    TH1D* histogram;
                    if(!(histogram = anaData[category->name].Signal().GetPtr<TH1D>(hist.name)))
                        continue;

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
                                         const std::string& hist_name_down)
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
            AnaDataForDataCategory& anaDataForCategory = fullAnaDataEntry.second;
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const DataCategory* dataCategory : dataCategoryCollection.GetCategories(DataCategoryType::Limits)) {
                if(!dataCategory->datacard.size())
                    throw exception("Empty datacard name for data category '") << dataCategory->name << "'.";
                FlatAnalyzerData& anaData = anaDataForCategory[dataCategory->name].Signal();
                TH1D* hist_orig = anaData.GetPtr<TH1D>(hist_name);
                if(!hist_orig)
                    throw exception("Datacard histogram '") << hist_name << "' not found for data category '"
                                                            << dataCategory->name << "'.";
                std::shared_ptr<TH1D> hist(static_cast<TH1D*>(hist_orig->Clone()));
                hist->Scale(dataCategory->limits_sf);
                hist->Write(dataCategory->datacard.c_str());
                const std::string namePrefix = dataCategory->datacard + "_CMS_scale_t_" + channel_name + "_8TeV";
                const std::string nameDown = namePrefix + "Down";
                const std::string nameUp = namePrefix + "Up";

                TH1D* hist_up_orig = anaData.GetPtr<TH1D>(hist_name_up);
                if(!hist_up_orig)
                    throw exception("Datacard histogram '") << hist_name_up << "' not found for data category '"
                                                            << dataCategory->name << "'.";
                std::shared_ptr<TH1D> hist_up(static_cast<TH1D*>(hist_up_orig->Clone()));
                hist_up->Scale(dataCategory->limits_sf);
                hist_up->Write(nameUp.c_str());
                TH1D* hist_down_orig = anaData.GetPtr<TH1D>(hist_name_down);
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
                AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
                if( TH1D* histogram = anaData[dataCategory->name].Signal().GetPtr<TH1D>(hist.name) ) {
                    typedef std::pair<Int_t, Int_t> limit_pair;
                    const limit_pair limits = includeOverflow ? limit_pair(0, histogram->GetNbinsX() + 1)
                                                              : limit_pair(1, histogram->GetNbinsX());
                    double error = 0.;
                    const double integral = histogram->IntegralAndError(limits.first, limits.second, error);
                    of << MakeStringRepresentationForValueWithError(integral, error, includeError) << sep;
                }
                else
                    of << "NaN" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;
    }

    void EstimateSumBkg(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                             const analysis::HistogramDescriptor& hist)
    {
        const analysis::DataCategory& sumBkg = dataCategoryCollection.GetUniqueCategory(DataCategoryType::Sum);
        for(const DataCategory* bkgCategory : dataCategoryCollection.GetCategories(DataCategoryType::Background)) {
            root_ext::SmartHistogram<TH1D>* bkg_hist;
            if(!(bkg_hist = anaData[bkgCategory->name].Signal().GetPtr<TH1D>(hist.name))) continue;

            if(TH1D* sum_hist = anaData[sumBkg.name].Signal().GetPtr<TH1D>(hist.name))
                sum_hist->Add(bkg_hist);
            else
                anaData[sumBkg.name].Signal().Clone(*bkg_hist);
        }
    }

protected:
    std::string inputPath;
    std::string outputFileName;
    DataCategoryCollection dataCategoryCollection;
    std::vector<HistogramDescriptor> histograms;
    FullAnaData fullAnaData;
    bool WjetsData;
};

} // namespace analysis
