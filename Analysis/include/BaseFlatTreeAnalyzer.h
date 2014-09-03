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

#include "AnalyzerData.h"
#include "FlatTree.h"
#include "PrintTools/include/RootPrintToPdf.h"
#include "../include/exception.h"

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
    root_ext::Range xRange;
    unsigned rebinFactor;
    unsigned nBins;
    bool useLogY;
};

typedef std::vector<DataSource> DataSourceVector;

struct DataCategory {
    std::string name;
    std::string title;
    EColor color;

    DataSourceVector sources;
};

typedef std::vector<DataCategory> DataCategoryCollection;

static const std::map<std::string, EColor> colorMapName = {{"white",kWhite}, {"black",kBlack}, {"gray",kGray},
                                                           {"red",kRed}, {"green",kGreen}, {"blue",kBlue},
                                                           {"yellow",kYellow}, {"magenta",kMagenta}, {"cyan",kCyan},
                                                           {"orange",kOrange}, {"spring",kSpring}, {"teal",kTeal},
                                                           {"azure",kAzure}, {"violet",kViolet},{"pink",kPink}};

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
    s << "Name: " << hist.name << ", Title: " << hist.title << ", xmin: " << hist.xRange.min << ", xmax: " <<
         hist.xRange.max << ", rebinFactor: " << hist.rebinFactor << ", useLog: " << hist.useLogY  ;
    return s;
}

enum class EventType_QCD { Unknown, OS_Isolated, OS_NotIsolated, SS_Isolated, SS_NotIsolated };
enum class EventType_Wjets { Unknown, Signal, HighMissingEt };
enum class EventCategory { Unknown, TwoJets_ZeroBtag, TwoJets_OneBtag, TwoJets_TwoBtag };

static const std::map<EventCategory, std::string> eventCategoryMapName =
          { { EventCategory::Unknown, "Unknown" }, { EventCategory::TwoJets_ZeroBtag, "TwoJets_ZeroBtag" },
          { EventCategory::TwoJets_OneBtag, "TwoJets_OneBtag"}, { EventCategory::TwoJets_TwoBtag, "TwoJets_TwoBtag" } };

std::ostream& operator<<(std::ostream& s, const EventCategory& eventCategory) {
    s << eventCategoryMapName.at(eventCategory);
    return s;
}

class FlatAnalyzerData : public root_ext::AnalyzerData {
public:

    TH1D_ENTRY(m_sv, 30, 0.0, 300.0)
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
                         const std::string& outputFileName, const std::string& _signalName,
                         const std::string& _dataName)
        : inputPath(_inputPath), signalName(_signalName), dataName(_dataName), printer(outputFileName)
    {
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
                ProcessDataSource(category.name, source);
            }
        }

        std::cout << "Estimating QCD... " << std::endl;
        EstimateQCD();
        std::cout << "Estimating Wjets... " << std::endl;
        EstimateWjets();
        std::cout << "Printing stacked plots... " << std::endl;
        PrintStackedPlots();
    }

protected:
    void ProcessDataSource(const std::string& dataCategoryName, const DataSource& dataSource)
    {
        for(size_t current_entry = 0; current_entry < dataSource.tree->GetEntries(); ++current_entry) {
            dataSource.tree->GetEntry(current_entry);
            const ntuple::Flat& event = dataSource.tree->data;
            const EventType_QCD eventTypeQCD = DetermineEventTypeForQCD(event);
            const EventType_Wjets eventTypeWjets = DetermineEventTypeForWjets(event);
            const EventCategory eventCategory = DetermineEventCategory(event);
            if(eventCategory == EventCategory::Unknown) continue;
            const bool isData = dataCategoryName.find("DATA") != std::string::npos;
            const double weight = isData ? 1 : event.weight * dataSource.scale_factor;
            FillHistograms(dataCategoryName, event, eventTypeQCD, eventTypeWjets, eventCategory, weight);
        }
    }

    virtual EventType_QCD DetermineEventTypeForQCD(const ntuple::Flat& event) = 0;
    virtual EventType_Wjets DetermineEventTypeForWjets(const ntuple::Flat& event) = 0;
    virtual EventCategory DetermineEventCategory(const ntuple::Flat& event) = 0;

    virtual void EstimateQCD()
    {
        using analysis::EventType_QCD;
        using analysis::EventCategory;

        analysis::DataCategory qcd;
        qcd.name = "QCD";
        qcd.title = "QCD";
        qcd.color = kGray;
        for (auto& fullAnaDataEntry : fullAnaData){
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const analysis::HistogramDescriptor& hist : histograms){
                const analysis::DataCategory& data = FindCategory("DATA");
                if(!anaData[data.name].QCD[EventType_QCD::SS_Isolated].Contains(hist.name)) continue;
                TH1D& histogram = anaData[qcd.name].QCD[EventType_QCD::OS_Isolated].Clone(
                            anaData[data.name].QCD[EventType_QCD::SS_Isolated].Get<TH1D>(hist.name));
                for (const analysis::DataCategory& category : categories){
                    const bool isData = category.name.find("DATA") != std::string::npos;
                    const bool isSignal = category.name.find("SIGNAL") != std::string::npos;
                    if (isData || isSignal) continue;
                    if(!anaData[category.name].QCD[EventType_QCD::SS_Isolated].Contains(hist.name)) continue;
                    TH1D& nonQCD_hist = anaData[category.name].QCD[EventType_QCD::SS_Isolated].Get<TH1D>(hist.name);
                    histogram.Add(&nonQCD_hist,-1);
                }
                histogram.Scale(1.06);
            }
        }
        categories.push_back(qcd);
    }

    virtual void EstimateWjets() = 0;

    void FillHistograms(const std::string& dataCategoryName, const ntuple::Flat& event, EventType_QCD eventTypeQCD,
                        EventType_Wjets eventTypeWjets, EventCategory eventCategory, double weight)
    {
        fullAnaData[eventCategory][dataCategoryName].QCD[eventTypeQCD].m_sv().Fill(event.m_sv, weight);
        fullAnaData[eventCategory][dataCategoryName].Wjets[eventTypeWjets].m_sv().Fill(event.m_sv, weight);
    }

    void PrintStackedPlots()
    {
        for(auto& fullAnaDataEntry : fullAnaData) {
            const EventCategory eventCategory = fullAnaDataEntry.first;
            AnaDataForDataCategory& anaData = fullAnaDataEntry.second;
            for (const HistogramDescriptor& hist : histograms) {
                std::ostringstream ss_title;
                ss_title << eventCategory << ": " << hist.title;

                page.side.use_log_scaleY = hist.useLogY;
                page.side.xRange = hist.xRange;
                page.side.fit_range_x = false;
                page.title = ss_title.str();
                page.side.axis_titleX = hist.Xaxis_title;
                page.side.axis_titleY = hist.Yaxis_title;
                page.side.layout.has_stat_pad = true;

                std::shared_ptr<THStack> stack = std::shared_ptr<THStack>(new THStack(hist.name.c_str(),hist.title.c_str()));

                TLegend* leg = new TLegend ( 0, 0.6, 1, 1.0);
                leg->SetFillColor(0);
                leg->SetTextSize(0.1);
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
                ll->AddText(0.65, 0.6, text);

                TH1D* data_histogram = nullptr;
                for(const DataCategory& category : categories) {
                    if(!anaData[category.name].QCD[EventType_QCD::OS_Isolated].Contains(hist.name)) continue;
                    const bool isData = category.name.find("DATA") != std::string::npos;
                    TH1D& histogram = anaData[category.name].QCD[EventType_QCD::OS_Isolated].Get<TH1D>(hist.name);
                    histogram.SetLineColor(category.color);
                    std::string legend_option = "f";
                    if (!isData){
                        histogram.SetFillColor(category.color);
                        stack->Add(&histogram);
                    }
                    else {
                        legend_option = "p";
                        data_histogram = &histogram;
                    }
                    leg->AddEntry(&histogram, category.title.c_str(), legend_option.c_str());
                }
                printer.PrintStack(page, *stack, *data_histogram, *leg, *ll);
            }
        }
    }

    const analysis::DataCategory& FindCategory(const std::string& prefix) const
    {
        for (const analysis::DataCategory& category : categories){
            const bool foundPrefix = category.name.find(prefix) != std::string::npos;
            if (foundPrefix)
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
            if(category.name.find("SIGNAL") != std::string::npos) {
                const size_t sub_name_pos = category.name.find(' ');
                const std::string sub_name = category.name.substr(sub_name_pos + 1);
                if(sub_name != signalName)
                    continue;
            }
            if(category.name.find("DATA") != std::string::npos) {
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
            ss >> hist.xRange.min;
            ss >> hist.xRange.max;
            ss >> hist.rebinFactor;
            ss >> hist.nBins;
            ss >> hist.useLogY;
            histograms.push_back(hist);
          }
    }

protected:
    std::string inputPath;
    std::string signalName, dataName;
    DataCategoryCollection categories;
    std::vector<HistogramDescriptor> histograms;
    Printer printer;
    root_ext::SingleSidedPage page;
    FullAnaData fullAnaData;
};

} // namespace analysis