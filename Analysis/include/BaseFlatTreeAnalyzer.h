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
#include <list>
#include <TColor.h>
#include <TLorentzVector.h>

#include "AnalyzerData.h"
#include "FlatTree.h"
#include "AnalysisMath.h"
#include "PrintTools/include/RootPrintToPdf.h"
#include "Htautau_Summer13.h"
#include "exception.h"

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
enum class EventCategory { Inclusive, TwoJets_ZeroBtag, TwoJets_OneBtag, TwoJets_TwoBtag };

static const std::map<EventCategory, std::string> eventCategoryMapName =
          { { EventCategory::Inclusive, "Inclusive" }, { EventCategory::TwoJets_ZeroBtag, "TwoJets_ZeroBtag" },
          { EventCategory::TwoJets_OneBtag, "TwoJets_OneBtag"}, { EventCategory::TwoJets_TwoBtag, "TwoJets_TwoBtag" } };
typedef std::vector<EventCategory> EventCategoryVector;

std::ostream& operator<<(std::ostream& s, const EventCategory& eventCategory) {
    s << eventCategoryMapName.at(eventCategory);
    return s;
}

static const std::vector<double> mass_bins = { 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140., 150,
                                               160, 170, 180, 190, 200, 225, 250, 275, 300 };

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:

    TH1D_ENTRY_CUSTOM(m_sv, mass_bins)
    TH1D_ENTRY(m_bb, 15, 0, 600)
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
        std::cout << "Saving tables and printing stacked plots... " << std::endl;
        PrintTables("comma", ",");
        PrintTables("semicolon", ";");
        PrintStackedPlots();
    }

protected:
    void ProcessDataSource(const DataCategory& dataCategory, const DataSource& dataSource)
    {
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
            }
        }
    }

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
        if(event.mass_Bjets.size() >= 2) {
            std::vector<TLorentzVector> b_momentums(2);
            for(size_t n = 0; n < b_momentums.size(); ++n)
                b_momentums.at(n).SetPtEtaPhiM(event.pt_Bjets.at(n), event.eta_Bjets.at(n), event.phi_Bjets.at(n),
                                               event.mass_Bjets.at(n));
            const double mass_bb = (b_momentums.at(0) + b_momentums.at(1)).M();
            anaData.m_bb().Fill(mass_bb, weight);
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
                page.side.axis_titleX = hist.Xaxis_title;
                page.side.axis_titleY = hist.Yaxis_title;
                page.side.layout.has_stat_pad = true;
                page.side.layout.main_pad.right_top.x=0.75;
                page.side.layout.main_pad.right_top.y=0.85;
                page.side.layout.stat_pad.left_bottom.x=0.755;
                page.side.layout.stat_pad.left_bottom.y=0.3;

                std::shared_ptr<THStack> stack = std::shared_ptr<THStack>(new THStack(hist.name.c_str(),hist.title.c_str()));

                TLegend* leg = new TLegend ( 0, 0.6, 1, 1.0);
                leg->SetFillColor(0);
                leg->SetTextSize(0.055);
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
                    if(!anaData[category.name].QCD[EventType_QCD::OS_Isolated].Contains(hist.name)) continue;
                    if(category.IsReference()) continue;
                    TH1D& histogram = anaData[category.name].QCD[EventType_QCD::OS_Isolated].Get<TH1D>(hist.name);
                    histogram.SetLineColor(colorMapName.at("black"));
                    histogram.SetLineWidth(1.);
                    //InitHist(&histogram,"","",category.color,1001);
                    ReweightWithBinWidth(histogram);

                    std::string legend_option = "f";
                    if (!category.IsData()){
                        histogram.SetFillColor(category.color);
                        //histogram.SetBinErrorOption(0);
                        if(category.IsSignal()){
//                            histogram.Scale(10);
                        }
                        stack->Add(&histogram);
                    }
                    else {
                        legend_option = "LP";

                        if (isBlind){
                            TH1D* refilledData_histogram = refillData(&histogram);
                            data_histogram = refilledData_histogram;
                        }
                        else
                            data_histogram = &histogram;

                    }
                    if (category.IsSignal()){
                        //InitSignal(&histogram);
                        leg->AddEntry(&histogram , category.title.c_str() , "F" );
                        //leg->SetLineStyle(2);
                        leg->SetLineColor(category.color);
                    }
                    else
                        leg->AddEntry(&histogram, category.title.c_str(), legend_option.c_str());
                }
                //InitData(data_histogram);
                printer.PrintStack(page, *stack, *data_histogram, *leg, *ll);
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
            if(category.IsSignal()) {
                const size_t sub_name_pos = category.name.find(' ');
                const std::string sub_name = category.name.substr(sub_name_pos + 1);
                if(sub_name != signalName)
                    continue;
            }
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

    float blinding_SM(float mass){ return (100<mass && mass<150); }

    TH1D* refillData(TH1D* histogram)
    {
        TH1D* histData = (TH1D*)histogram->Clone();
        histData->Clear();
        for(int i=1; i<=histData->GetNbinsX(); ++i){
            histData->SetBinContent(i, blinding_SM (histData->GetBinCenter(i)) ? 0. : histData->GetBinContent(i));
        }
        return histData;

    }

    void PrintTables(const std::string& name_suffix, const std::string& sep)
    {
        std::ofstream of(outputFileName + "_" + name_suffix + ".csv");
        of << std::setprecision(1) << std::fixed;

        for(const auto& hist : histograms) {
            of << hist.title << sep;
            for (const auto& fullAnaDataEntry : fullAnaData) {
                const EventCategory& eventCategory = fullAnaDataEntry.first;
                of << eventCategory << sep;
            }
            of << std::endl;
            for (auto& dataCategory : categories) {
                of << dataCategory.title << sep;
                for (auto& fullAnaDataEntry : fullAnaData) {
                    AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
                    if( TH1D* histogram = anaData[dataCategory.name].QCD[EventType_QCD::OS_Isolated].GetPtr<TH1D>(hist.name) )
                        of << histogram->Integral() << sep;
                    else
                        of << "NaN" << sep;
                }
                of << std::endl;
            }
            of << std::endl << std::endl;
        }
    }

    void SetLegendStyle(TLegend* leg)
    {
      leg->SetFillStyle (0);
      leg->SetFillColor (0);
      leg->SetBorderSize(0);
    }

    void InitData(TH1* hist)
    {
      hist->SetMarkerStyle(20.);
      hist->SetMarkerSize (1.3);
      hist->SetLineWidth  ( 1.);
    }

    void InitSignal(TH1 *hist)
    {
      hist->SetLineStyle(2);
      hist->SetLineWidth(2.);
      hist->SetLineColor(kBlue);
    }

    void CMSPrelim(const char* dataset, const char* channel, double lowX, double lowY)
    {
      TPaveText* lumi  = new TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC");
      lumi->SetBorderSize(   0 );
      lumi->SetFillStyle(    0 );
      lumi->SetTextAlign(   12 );
      lumi->SetTextSize ( 0.035 );
      lumi->SetTextColor(    1 );
      lumi->SetTextFont (   62 );
      lumi->AddText(dataset);
      lumi->Draw();

      TPaveText* chan     = new TPaveText(lowX+0.68, lowY+0.061, lowX+0.80, lowY+0.161, "NDC");
      chan->SetBorderSize(   0 );
      chan->SetFillStyle(    0 );
      chan->SetTextAlign(   12 );
      chan->SetTextSize ( 0.04 );
      chan->SetTextColor(    1 );
      chan->SetTextFont (   62 );
      chan->AddText(channel);
      chan->Draw();
    }

    void InitHist(TH1 *hist, const char *xtit, const char *ytit, int color, int style)
    {
      hist->SetXTitle(xtit);
      hist->SetYTitle(ytit);
      hist->SetLineColor(kBlack);
      hist->SetLineWidth(    1.);
      hist->SetFillColor(color );
      hist->SetFillStyle(style );
      hist->SetTitleSize  (0.055,"Y");
      hist->SetTitleOffset(1.600,"Y");
      hist->SetLabelOffset(0.014,"Y");
      hist->SetLabelSize  (0.040,"Y");
      hist->SetLabelFont  (42   ,"Y");
      hist->SetTitleSize  (0.055,"X");
      hist->SetTitleOffset(1.300,"X");
      hist->SetLabelOffset(0.014,"X");
      hist->SetLabelSize  (0.040,"X");
      hist->SetLabelFont  (42   ,"X");
      hist->SetMarkerStyle(20);
      hist->SetMarkerColor(color);
      hist->SetMarkerSize (0.6);
      hist->GetYaxis()->SetTitleFont(42);
      hist->GetXaxis()->SetTitleFont(42);
      hist->SetTitle("");
      return;
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
