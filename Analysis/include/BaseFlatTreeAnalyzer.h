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

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>

#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatTree.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"

#include "Htautau_Summer13.h"

namespace analysis {

struct DataSource {
    std::string file_name;
    double scale_factor;
    TFile* file;
    std::shared_ptr<ntuple::FlatTree> tree;
};

struct HistogramDescriptor {
    std::string name;
    std::string title;
    std::string Xaxis_title;
    std::string Yaxis_title;
    bool useLogY;
};

typedef std::vector<DataSource> DataSourceVector;

struct DataCategory {
    std::string name;
    std::string title;
    EColor color;

    DataSourceVector sources;

public:
    bool IsData() const { return NameContains("DATA"); }
    bool IsSignal() const { return NameContains("SIGNAL"); }
    bool IsReference() const { return NameContains("REFERENCE"); }
    bool IsVirtual() const { return NameContains("VIRTUAL"); }
    bool IsForLimitsOnly() const { return NameContains("LIMITS"); }
    bool IsSumBkg() const { return NameContains("SUM"); }

    bool NameContains(const std::string& substring) const { return name.find(substring) != std::string::npos; }
};

typedef std::list<DataCategory> DataCategoryCollection;

static const std::map<std::string, EColor> colorMapName = {{"white",kWhite}, {"black",kBlack}, {"gray",kGray},
                                                           {"red",kRed}, {"green",kGreen}, {"blue",kBlue},
                                                           {"yellow",kYellow}, {"magenta",kMagenta}, {"cyan",kCyan},
                                                           {"orange",kOrange}, {"spring",kSpring}, {"teal",kTeal},
                                                           {"azure",kAzure}, {"violet",kViolet},{"pink",kPink},
                                                           {"pink_custom", (EColor) TColor::GetColor(250,202,255)},
                                                           {"red_custom", (EColor) TColor::GetColor(222,90,106)},
                                                           {"violet_custom", (EColor) TColor::GetColor(155,152,204)},
                                                           {"yellow_custom", (EColor) TColor::GetColor(248,206,104)}};

std::istream& operator>>(std::istream& s, EColor& color){
    std::string name;
    s >> name;
    if (!colorMapName.count(name)){
        std::ostringstream ss;
        ss << "undefined color: '" << name ;
        throw std::runtime_error(ss.str());
    }
    color = colorMapName.at(name);
    return s;
}

std::ostream& operator<<(std::ostream& s, const EColor& color){
    for(const auto& entry : colorMapName) {
        if(entry.second == color) {
            s << entry.first;
            return s;
        }
    }
    s << "Unknown color " << color;
    return s;
}

std::ostream& operator<<(std::ostream& s, const DataSource& source){
    s << "File: " << source.file_name << ", SF: " << source.scale_factor;
    return s;
}

std::ostream& operator<<(std::ostream& s, const DataCategory& category){
    s << "Name: " << category.name << ", Title: '" << category.title << "', Color: " << category.color << std::endl;
    for(const DataSource& source : category.sources)
        s << source << std::endl;
    return s;
}

std::ostream& operator<<(std::ostream& s, const HistogramDescriptor& hist){
    s << "Name: " << hist.name << ", Title: " << hist.title << ", useLog: " << hist.useLogY  ;
    return s;
}

enum class EventType_QCD { Unknown, OS_Isolated, OS_NotIsolated, SS_Isolated, SS_NotIsolated };
enum class EventType_Wjets { Unknown, Signal, HighMt };
enum class EventCategory { Inclusive, OneJet_ZeroBtag, OneJet_OneBtag, TwoJets_ZeroBtag, TwoJets_OneBtag, TwoJets_TwoBtag };

static const std::map<EventCategory, std::string> eventCategoryMapName =
          { { EventCategory::Inclusive, "Inclusive" }, { EventCategory::OneJet_ZeroBtag, "1jet0btag" },
            { EventCategory::OneJet_OneBtag, "1jet1btag" }, { EventCategory::TwoJets_ZeroBtag, "2jets0btag" },
          { EventCategory::TwoJets_OneBtag, "2jets1btag"}, { EventCategory::TwoJets_TwoBtag, "2jets2btag" } };
typedef std::vector<EventCategory> EventCategoryVector;

std::ostream& operator<<(std::ostream& s, const EventCategory& eventCategory) {
    s << eventCategoryMapName.at(eventCategory);
    return s;
}

std::wostream& operator<<(std::wostream& s, const EventCategory& eventCategory) {
    const std::string str = eventCategoryMapName.at(eventCategory);
    s << std::wstring(str.begin(), str.end());
    return s;
}


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
    TH1D_ENTRY(m_ttbb_nomet, 100, 0, 1000)
    TH1D_ENTRY(DeltaPhi_tt, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_bb, 80, -4, 4)
    TH1D_ENTRY(DeltaPhi_hh, 80, -4, 4)
    TH1D_ENTRY(DeltaR_tt, 60, 0, 6)
    TH1D_ENTRY(DeltaR_bb, 60, 0, 6)
    TH1D_ENTRY(DeltaR_hh, 60, 0, 6)
};

class BaseFlatTreeAnalyzer {
public:
    typedef root_ext::PdfPrinter Printer;


    typedef std::map<EventType_QCD, FlatAnalyzerData> AnaDataQCD;
    typedef std::map<EventType_Wjets, FlatAnalyzerData> AnaDataWjets;

    struct AnaDataEventTypes {
        AnaDataQCD QCD;
        AnaDataWjets Wjets;
    };

    typedef std::map<std::string, AnaDataEventTypes> AnaDataForDataCategory;
    typedef std::map<EventCategory, AnaDataForDataCategory> FullAnaData;

    BaseFlatTreeAnalyzer(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                         const std::string& _outputFileName, const std::string& _signalName,
                         const std::string& _dataName, bool _WjetsData = false, bool _isBlind=false)
        : inputPath(_inputPath), signalName(_signalName), dataName(_dataName), outputFileName(_outputFileName),
          WjetsData(_WjetsData), isBlind(_isBlind)
    {
        TH1::SetDefaultSumw2();

        ReadSourceCfg(source_cfg);
        ReadHistCfg(hist_cfg);
    }

    void Run()
    {
        std::cout << "Processing sources... " << std::endl;
        for(DataCategory& category : categories) {
            std::cout << category << std::endl;
            for (DataSource& source : category.sources){
                const std::string fullFileName = inputPath + "/" + source.file_name;
                source.file = new TFile(fullFileName.c_str(), "READ");
                if(source.file->IsZombie()) {
                    std::ostringstream ss;
                    ss << "Input file '" << source.file_name << "' not found.";
                    throw std::runtime_error(ss.str());
                }
                source.tree = std::shared_ptr<ntuple::FlatTree>(new ntuple::FlatTree(*source.file, "flatTree"));
                ProcessDataSource(category, source);
            }
        }

        std::cout << "Estimating QCD and Wjets... " << std::endl;
        for (auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const auto& hist : histograms) {
                EstimateQCD(eventCategory, anaData, hist);
                if (WjetsData) EstimateWjets(eventCategory, anaData, hist);
            }
        }
        std::cout << "Estimating Sum od backgrounds... " << std::endl;
        for (auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const auto& hist : histograms) {
                EstimateSumBkg(eventCategory, anaData, hist);
            }
        }

        std::cout << "Saving tables and printing stacked plots... " << std::endl;
        static const std::wstring comma = L",", semicolon = L";";
        PrintTables("comma", comma);
        PrintTables("semicolon", semicolon);
        ProduceFileForLimitsCalculation();
        PrintStackedPlots();
    }

protected:
    void ProcessDataSource(const DataCategory& dataCategory, const DataSource& dataSource)
    {
        const analysis::DataCategory& Ztautau = FindCategory("LIMITS Ztautau");
        const analysis::DataCategory& DYJets = FindCategory("DYJets");

        for(size_t current_entry = 0; current_entry < dataSource.tree->GetEntries(); ++current_entry) {
            dataSource.tree->GetEntry(current_entry);
            const ntuple::Flat& event = dataSource.tree->data;
            const EventType_QCD eventTypeQCD = DetermineEventTypeForQCD(event);
            const EventType_Wjets eventTypeWjets = DetermineEventTypeForWjets(event);
            const EventCategoryVector eventCategories = DetermineEventCategories(event);
            const double weight = dataCategory.IsData() ? 1 : event.weight * dataSource.scale_factor;
            for(auto eventCategory : eventCategories) {
                FillHistograms(fullAnaData[eventCategory][dataCategory.name].QCD[eventTypeQCD], event, weight);
                FillHistograms(fullAnaData[eventCategory][dataCategory.name].Wjets[eventTypeWjets], event, weight);
                if(dataCategory.name == DYJets.name && std::abs(event.pdgId_2_MC) == particles::tau.RawCode())
                    FillHistograms(fullAnaData[eventCategory][Ztautau.name].QCD[eventTypeQCD], event, weight);
            }
        }
    }

    virtual const std::string& ChannelName() = 0;

    virtual EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) = 0;
    virtual EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) = 0;

    virtual EventCategoryVector DetermineEventCategories(const ntuple::Flat& event)
    {
        EventCategoryVector categories;
        categories.push_back(EventCategory::Inclusive);

        std::vector<Float_t> goodCVSvalues;
        for (unsigned i = 0; i < event.eta_Bjets.size(); ++i){
            if ( std::abs(event.eta_Bjets.at(i)) >= cuts::Htautau_Summer13::btag::eta) continue;
            goodCVSvalues.push_back(event.csv_Bjets.at(i));
        }

        if (goodCVSvalues.size() == 1) {
            if (goodCVSvalues.at(0) <= cuts::Htautau_Summer13::btag::CSVT )
                categories.push_back(EventCategory::OneJet_ZeroBtag);
            else
                categories.push_back(EventCategory::OneJet_OneBtag);
        }

        if (goodCVSvalues.size() >= 2) {
            if (goodCVSvalues.at(0) <= cuts::Htautau_Summer13::btag::CSVM )
                categories.push_back(EventCategory::TwoJets_ZeroBtag);
            else if ( goodCVSvalues.at(1) <= cuts::Htautau_Summer13::btag::CSVM )
                categories.push_back(EventCategory::TwoJets_OneBtag);
            else
                categories.push_back(EventCategory::TwoJets_TwoBtag);
        }
        return categories;
    }

    virtual void EstimateQCD(EventCategory eventCategory, AnaDataForDataCategory& anaData,
                             const HistogramDescriptor& hist) = 0;
    virtual void EstimateWjets(EventCategory eventCategory, AnaDataForDataCategory& anaData,
                               const HistogramDescriptor& hist) = 0;

    virtual void FillHistograms(FlatAnalyzerData& anaData, const ntuple::Flat& event, double weight)
    {
        anaData.m_sv().Fill(event.m_sv, weight);
        anaData.m_sv_up().Fill(event.m_sv_Up, weight);
        anaData.m_sv_down().Fill(event.m_sv_Down, weight);
        anaData.pt_1().Fill(event.pt_1, weight);
        anaData.eta_1().Fill(event.eta_1, weight);
        anaData.pt_2().Fill(event.pt_2, weight);
        anaData.eta_2().Fill(event.eta_2, weight);
        TLorentzVector first_cand, second_cand;
        first_cand.SetPtEtaPhiE(event.pt_1,event.eta_1,event.phi_1,event.energy_1);
        second_cand.SetPtEtaPhiE(event.pt_2,event.eta_2,event.phi_2,event.energy_2);
        anaData.DeltaPhi_tt().Fill(first_cand.DeltaPhi(second_cand), weight);
        anaData.DeltaR_tt().Fill(first_cand.DeltaR(second_cand), weight);
        TLorentzVector Htt = first_cand + second_cand;
        anaData.pt_H_tt().Fill(Htt.Pt(),weight);
        anaData.m_vis().Fill(Htt.M(),weight);
        TLorentzVector MET;
        MET.SetPtEtaPhiE(event.mvamet,0,event.mvametphi,0);
        if(event.mass_Bjets.size() >= 2) {
            std::vector<TLorentzVector> b_momentums(2);
            for(size_t n = 0; n < b_momentums.size(); ++n)
                b_momentums.at(n).SetPtEtaPhiM(event.pt_Bjets.at(n), event.eta_Bjets.at(n), event.phi_Bjets.at(n),
                                               event.mass_Bjets.at(n));
            anaData.pt_b1().Fill(b_momentums.at(0).Pt(), weight);
            anaData.eta_b1().Fill(b_momentums.at(0).Eta(), weight);
            anaData.pt_b2().Fill(b_momentums.at(1).Pt(), weight);
            anaData.eta_b2().Fill(b_momentums.at(1).Eta(), weight);
            anaData.DeltaPhi_bb().Fill(b_momentums.at(0).DeltaPhi(b_momentums.at(1)), weight);
            anaData.DeltaR_bb().Fill(b_momentums.at(0).DeltaR(b_momentums.at(1)), weight);
            const TLorentzVector Hbb = b_momentums.at(0) + b_momentums.at(1);
            anaData.pt_H_bb().Fill(Hbb.Pt(),weight);
            anaData.m_bb().Fill(Hbb.M(), weight);
            anaData.DeltaPhi_hh().Fill(Htt.DeltaPhi(Hbb), weight);
            anaData.DeltaR_hh().Fill(Htt.DeltaR(Hbb), weight);
            const TLorentzVector Candidate_ttbb = Hbb + Htt + MET;
            anaData.m_ttbb().Fill(Candidate_ttbb.M(), weight);
            anaData.pt_H_hh().Fill(Candidate_ttbb.Pt(), weight);
            const TLorentzVector Candidate_ttbb_noMET = Hbb + Htt;
            anaData.m_ttbb_nomet().Fill(Candidate_ttbb_noMET.M(), weight);
        }
    }

    void PrintStackedPlots()
    {
        Printer printer(outputFileName + ".pdf");
        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const HistogramDescriptor& hist : histograms) {
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist.title;

                page.side.use_log_scaleY = hist.useLogY;
                page.side.fit_range_x = false;
                page.title = ss_title.str();
                page.layout.has_title = false;
                page.side.axis_titleX = hist.Xaxis_title;
                page.side.axis_titleY = hist.Yaxis_title;
                page.side.layout.has_stat_pad = false;
                page.side.layout.main_pad.right_top.x=1;
                page.side.layout.main_pad.right_top.y=1.;
                page.side.layout.main_pad.left_bottom.x=0.;
                page.side.layout.main_pad.left_bottom.y=0;
                page.side.layout.stat_pad.left_bottom.x=0.7;
                page.side.layout.stat_pad.left_bottom.y=0.3;

                std::shared_ptr<THStack> stack = std::shared_ptr<THStack>(new THStack(hist.name.c_str(),hist.title.c_str()));

                TLegend* leg = new TLegend (0.60, 0.65, 0.67, 0.90);
                leg->SetFillColor(0);
                leg->SetTextSize(0.035);
                SetLegendStyle(leg);
                TString lumist="19.7 fb^{-1}";
                TPaveText *ll = new TPaveText(0.15, 0.95, 0.95, 0.99, "NDC");
                ll->SetTextSize(0.03);
                ll->SetTextFont(42);
                ll->SetFillColor(0);
                ll->SetBorderSize(0);
                ll->SetMargin(0.01);
                ll->SetTextAlign(12); // align left
                TString text = "CMS Preliminary";
                ll->AddText(0.01,0.5,text);
                text = "#sqrt{s} = 8 TeV  L = ";
                text = text + lumist;
                ll->AddText(0.25, 0.6, text);

                TH1D* data_histogram = nullptr;
                for(const DataCategory& category : categories) {
                    TH1D* histogram;
                    if(!(histogram = anaData[category.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name))) continue;
                    if(category.IsReference() || (category.IsSignal() && !category.NameContains(signalName)) ||
                            category.IsSumBkg() || category.IsForLimitsOnly()) continue;
                    histogram->SetLineColor(colorMapName.at("black"));
                    histogram->SetLineWidth(1.);
                    ReweightWithBinWidth(*histogram);

                    std::string legend_option = "f";
                    if (!category.IsData()){
                        histogram->SetFillColor(category.color);
//                        if(category.IsSignal()) histogram.Scale(10);
                        stack->Add(histogram);
                    }
                    else {
                        legend_option = "LP";

                        if (isBlind){
                            TH1D* refilledData_histogram = refillData(histogram);
                            data_histogram = refilledData_histogram;
                        }
                        else
                            data_histogram = histogram;
                    }
                    if (category.IsSignal()){
                        leg->AddEntry(histogram , category.title.c_str() , "F" );
                        leg->SetLineColor(category.color);
                    }
                    else
                        leg->AddEntry(histogram, category.title.c_str(), legend_option.c_str());
                }
                if(!data_histogram) {
                    std::cerr << "WARNING: DATA histogram not found for " << hist.name << std::endl;
                    continue;
                }
                TH1D* bkg_sum = nullptr;
                for(const DataCategory& category : categories) {
                    if(!anaData[category.name].QCD[EventType_QCD::OS_Isolated].Contains(hist.name)) continue;
                    if(category.IsReference() || category.IsData() || category.IsSignal() || category.IsSumBkg()
                            || category.IsForLimitsOnly()) continue;
                    TH1D& bkg = anaData[category.name].QCD[EventType_QCD::OS_Isolated].Get<TH1D>(hist.name);
                    if(!bkg_sum)
                        bkg_sum = static_cast<TH1D*>(bkg.Clone());
                    else
                        bkg_sum->Add(&bkg);
                }

                TH1D* ratioData_Bkg = (TH1D*)data_histogram->Clone("ratioData_Bkg");
                if(bkg_sum)
                    ratioData_Bkg->Divide(bkg_sum);
                printer.PrintStack(page, *stack, *data_histogram, *ratioData_Bkg, *leg, *ll);
            }
        }
    }

    const analysis::DataCategory& FindCategory(const std::string& prefix) const
    {
        for (const analysis::DataCategory& category : categories){
            if (category.NameContains(prefix))
                return category;
        }
        throw analysis::exception("category not found : ") << prefix ;
    }

    void ProduceFileForLimitsCalculation()
    {
        static const std::map<EventCategory, std::string> categoryToDirectoryNameSuffix = {
            { EventCategory::Inclusive, "inclusive" }, { EventCategory::OneJet_ZeroBtag, "1jet0tag" },
            { EventCategory::OneJet_OneBtag, "1jet1tag" }, { EventCategory::TwoJets_ZeroBtag, "2jet0tag" },
            { EventCategory::TwoJets_OneBtag, "2jet1tag" }, { EventCategory::TwoJets_TwoBtag, "2jet2tag" }
        };

        static const std::vector< std::pair<std::string, std::string> > dataCategoriesForLimits = {
            { "VIRTUAL QCD", "QCD" }, { "TTbar", "TT" }, { "LIMITS VH125", "VH125" },
            { "LIMITS DiBoson", "VV" }, { "LIMITS Wjets", "W" }, { "ZJ", "ZJ" }, { "ZL", "ZL" },
            { "ZLL", "ZLL" }, { "LIMITS Ztautau", "ZTT" }, { "bbH100", "bbH100" }, { "bbH110", "bbH110" },
            { "bbH120", "bbH120" }, { "bbH130", "bbH130" }, { "bbH140", "bbH140" }, { "bbH160", "bbH160" },
            { "bbH180", "bbH180" }, { "bbH200", "bbH200" }, { "bbH250", "bbH250" }, { "bbH300", "bbH300" },
            { "bbH350", "bbH350" }, { "bbH400", "bbH400" }, { "bbH90", "bbH90" },
            { "DATA Tau", "data_obs" }, { "DATA TauPlusX", "data_obs" },
            { "ggAToZhToLLBB260", "ggAToZhToLLBB260" }, { "ggAToZhToLLBB270", "ggAToZhToLLBB270" },
            { "ggAToZhToLLBB280", "ggAToZhToLLBB280" }, { "ggAToZhToLLBB290", "ggAToZhToLLBB290" },
            { "ggAToZhToLLBB300", "ggAToZhToLLBB300" }, { "ggAToZhToLLBB310", "ggAToZhToLLBB310" },
            { "ggAToZhToLLBB320", "ggAToZhToLLBB320" }, { "ggAToZhToLLBB330", "ggAToZhToLLBB330" },
            { "ggAToZhToLLBB340", "ggAToZhToLLBB340" }, { "ggAToZhToLLBB350", "ggAToZhToLLBB350" },
            { "ggAToZhToLLTauTau260", "ggAToZhToLLTauTau260" }, { "ggAToZhToLLTauTau270", "ggAToZhToLLTauTau270" },
            { "ggAToZhToLLTauTau280", "ggAToZhToLLTauTau280" }, { "ggAToZhToLLTauTau290", "ggAToZhToLLTauTau290" },
            { "ggAToZhToLLTauTau300", "ggAToZhToLLTauTau300" }, { "ggAToZhToLLTauTau310", "ggAToZhToLLTauTau310" },
            { "ggAToZhToLLTauTau320", "ggAToZhToLLTauTau320" }, { "ggAToZhToLLTauTau330", "ggAToZhToLLTauTau330" },
            { "ggAToZhToLLTauTau340", "ggAToZhToLLTauTau340" }, { "ggAToZhToLLTauTau350", "ggAToZhToLLTauTau350" },
            { "LIMITS ggH125", "ggH125" }, { "LIMITS qqH125", "qqH125" },
            { "LIMITS ggHhh260", "ggHTohhTo2Tau2B260" }, { "LIMITS ggHhh270", "ggHTohhTo2Tau2B270" },
            { "LIMITS ggHhh280", "ggHTohhTo2Tau2B280" }, { "LIMITS ggHhh290", "ggHTohhTo2Tau2B290" },
            { "LIMITS ggHhh300", "ggHTohhTo2Tau2B300" }, { "LIMITS ggHhh310", "ggHTohhTo2Tau2B310" },
            { "LIMITS ggHhh320", "ggHTohhTo2Tau2B320" }, { "LIMITS ggHhh330", "ggHTohhTo2Tau2B330" },
            { "LIMITS ggHhh340", "ggHTohhTo2Tau2B340" }, { "LIMITS ggHhh350", "ggHTohhTo2Tau2B350" }
        };

        static const std::map<std::string, std::string> channelNameForFolder = {
            { "eTau", "eleTau" }, { "muTau", "muTau" }, { "tauTau", "tauTau" }
        };

        std::string channel_name = ChannelName();
        std::transform(channel_name.begin(), channel_name.end(), channel_name.begin(), ::tolower);

        std::shared_ptr<TFile> outputFile(new TFile((outputFileName + ".root").c_str(), "RECREATE"));
        outputFile->cd();
        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaDataForCategory = fullAnaDataEntry.second;
            if(!categoryToDirectoryNameSuffix.count(eventCategory)) continue;
            const std::string directoryName = channelNameForFolder.at(ChannelName()) + "_"
                    + categoryToDirectoryNameSuffix.at(eventCategory);
            outputFile->mkdir(directoryName.c_str());
            outputFile->cd(directoryName.c_str());
            for(const auto& limitDataCategory : dataCategoriesForLimits) {
                if(!anaDataForCategory.count(limitDataCategory.first)) {
                    std::cerr << "WARNING: Dataset '" << limitDataCategory.first << "' for event category '"
                              << eventCategory << "' not found." << std::endl;
                    continue;
                }
                FlatAnalyzerData& anaData = anaDataForCategory[limitDataCategory.first].QCD[EventType_QCD::OS_Isolated];
                anaData.m_sv().Write(limitDataCategory.second.c_str());
                const std::string namePrefix = limitDataCategory.second + "_CMS_scale_t_" + channel_name + "_8TeV";
                const std::string nameDown = namePrefix + "Down";
                const std::string nameUp = namePrefix + "Up";

                anaData.m_sv_down().Write(nameDown.c_str());
                anaData.m_sv_up().Write(nameUp.c_str());
            }
        }
        outputFile->Close();
    }

private:
    void ReadSourceCfg(const std::string& cfg_name)
    {
        std::ifstream cfg(cfg_name);
        std::shared_ptr<DataCategory> currentCategory;
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            if (cfgLine.at(0) == '[') {
                if(currentCategory)
                    categories.push_back(*currentCategory);
                currentCategory = std::shared_ptr<DataCategory>(new DataCategory());
                const size_t pos = cfgLine.find(']');
                currentCategory->name = cfgLine.substr(1, pos - 1);
                std::getline(cfg, currentCategory->title);
                std::string colorLine;
                std::getline(cfg,colorLine);
                std::istringstream ss(colorLine);
                ss >> currentCategory->color;
            }
            else if (currentCategory) {
                std::istringstream ss(cfgLine);
                DataSource source;
                ss >> source.file_name;
                ss >> source.scale_factor;
                currentCategory->sources.push_back(source);
            }
            else
                throw std::runtime_error("bad source file format");
          }
        if(currentCategory)
            categories.push_back(*currentCategory);
        DataCategoryCollection filteredCategories;
        for(const DataCategory& category : categories) {
//            if(category.IsSignal()) {
//                const size_t sub_name_pos = category.name.find(' ');
//                const std::string sub_name = category.name.substr(sub_name_pos + 1);
//                if(sub_name != signalName)
//                    continue;
//            }
            if(category.IsData()) {
                const size_t sub_name_pos = category.name.find(' ');
                const std::string sub_name = category.name.substr(sub_name_pos + 1);
                if(sub_name != dataName)
                    continue;
            }
            filteredCategories.push_back(category);
        }
        categories = filteredCategories;
    }

    void ReadHistCfg(const std::string& cfg_name)
    {
        std::ifstream cfg(cfg_name);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            std::istringstream ss(cfgLine);
            HistogramDescriptor hist;
            ss >> hist.name;
            ss >> hist.title;
            ss >> hist.Xaxis_title;
            ss >> hist.Yaxis_title;
            ss >> std::boolalpha >> hist.useLogY;
            histograms.push_back(hist);
          }
    }

    static void ReweightWithBinWidth(TH1D& histogram)
    {
        for(Int_t n = 1; n <= histogram.GetNbinsX(); ++n) {
            const double new_value = histogram.GetBinContent(n) / histogram.GetBinWidth(n);
            const double new_bin_error = histogram.GetBinError(n) / histogram.GetBinWidth(n);
            histogram.SetBinContent(n, new_value);
            histogram.SetBinError(n, new_bin_error);
        }
    }

    TH1D* refillData(TH1D* histogram)
    {
        static const std::vector< std::pair<double, double> > blindingRegions = { { 100, 150 }, { 250, 350 } };
        static const std::map<std::string, size_t> histogramsToBlind = {
            { "m_sv", 0 }, { "m_sv_up", 0 }, { "m_sv_down", 0 }, { "m_vis", 0 }, { "m_bb", 0 },
            { "m_ttbb", 1 }, { "m_ttbb_nomet", 1 }
        };

        static auto need_blind = [&](double mass) -> bool {
            if(!histogramsToBlind.count(histogram->GetName())) return false;
            const size_t regionId = histogramsToBlind.at(histogram->GetName());
            if(regionId >= blindingRegions.size())
                throw analysis::exception("Bad blinding region index = ") << regionId;
            const auto& region = blindingRegions.at(regionId);
            return mass > region.first && mass < region.second;
        };

        TH1D* histData = (TH1D*)histogram->Clone();
        histData->Clear();
        for(int i=1; i<=histData->GetNbinsX(); ++i){

            const bool blinding = need_blind(histData->GetBinCenter(i));
            histData->SetBinContent(i, blinding ? 0. : histogram->GetBinContent(i));
            histData->SetBinError(i, blinding ? 0. : histogram->GetBinError(i));
        }
        return histData;

    }

    void PrintTables(const std::string& name_suffix, const std::wstring& sep)
    {
        std::wofstream of(outputFileName + "_" + name_suffix + ".csv");
        of << std::setprecision(1) << std::fixed;

        for(const HistogramDescriptor& hist : histograms) {
            if (hist.name != "m_sv") continue;
            PrintHistogram(of, sep, hist, false,false);
            PrintHistogram(of, sep, hist, true, false);
            PrintHistogram(of, sep, hist, false,true);
            PrintHistogram(of, sep, hist, true, true);
        }
        of.flush();
        of.close();
    }

    void PrintHistogram(std::wostream& of, const std::wstring& sep, const HistogramDescriptor& hist, bool includeOverflow,
                        bool includeError)
    {
        of << std::wstring(hist.title.begin(), hist.title.end());
        if(includeOverflow){
            if (!includeError)
                of << " with overflow" << sep;
            else
                of << "with overflow and error" << sep;
        }
        else if (includeError)
            of << "without overflow and error" << sep;
        else
            of << sep;

        for (const auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory& eventCategory = fullAnaDataEntry.first;
            of << eventCategory << sep;
        }
        of << std::endl;
        for (DataCategory& dataCategory : categories) {
            if (dataCategory.IsReference()) continue;
            of << std::wstring(dataCategory.title.begin(), dataCategory.title.end()) << sep;
            for (auto& fullAnaDataEntry : fullAnaData) {
                AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
                if( TH1D* histogram = anaData[dataCategory.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name) ) {
                    if (!includeError){
                        const double integral = includeOverflow ? histogram->Integral(0, histogram->GetNbinsX() + 1)
                                                                : histogram->Integral();
                        of << integral << sep;
                    }
                    else if (includeError){
                        double error = 0.;
                        const double integral = includeOverflow ?
                                    histogram->IntegralAndError(0, histogram->GetNbinsX() + 1,error)
                                  : histogram->IntegralAndError(0, histogram->GetNbinsX(),error);
                        of << integral << L" \u00B1 " << error << sep;
                    }
                }
                else
                    of << "NaN" << sep;
            }
            of << std::endl;
        }
        of << std::endl << std::endl;

    }

    void SetLegendStyle(TLegend* leg)
    {
      leg->SetFillStyle (0);
      leg->SetFillColor (0);
      leg->SetBorderSize(0);
    }

    void EstimateSumBkg(analysis::EventCategory /*eventCategory*/, AnaDataForDataCategory& anaData,
                             const analysis::HistogramDescriptor& hist)
    {
        const analysis::DataCategory& sumBkg = FindCategory("SUM");

        std::cout << "created sumBkg" << std::endl;

        const analysis::DataCategory& ewk = FindCategory("EWK");
        root_ext::SmartHistogram<TH1D>* hist_ewk;
        if(!(hist_ewk = anaData[ewk.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name))) return;
        TH1D& histogram = anaData[sumBkg.name].QCD[EventType_QCD::OS_Isolated].Clone(*hist_ewk);
        for (const analysis::DataCategory& category : categories){
            if (category.IsReference() || category.IsData() || category.IsSignal() || category.name == ewk.name || category.IsSumBkg()
                    || category.IsForLimitsOnly()) continue;

            TH1D* otherBkg_hist;
            if(!(otherBkg_hist = anaData[category.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name))) continue;
            histogram.Add(otherBkg_hist);
        }
    }


protected:
    std::string inputPath;
    std::string signalName, dataName;
    std::string outputFileName;
    DataCategoryCollection categories;
    std::vector<HistogramDescriptor> histograms;
    root_ext::SingleSidedPage page;
    FullAnaData fullAnaData;
    bool WjetsData;
    bool isBlind;
};

} // namespace analysis
