/*!
 * \file Print_Stack.C
 * \brief Print stack with specified name superimposing several files.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-04-03 created
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

#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "../include/RootPrintToPdf.h"


struct DataSource {
    std::string file_name;
    double scale_factor;
    TFile* file;
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
    s >> name ;
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

class Print_Stack {
public:
    typedef root_ext::PdfPrinter Printer;

    static TH1D* CreateHistogram(const DataCategory& category, const HistogramDescriptor& hist, bool isData)
    {
        std::ostringstream name;
        name << category.name << "_" << hist.name << "_hist";
        TH1D* histogram = new TH1D(name.str().c_str(), name.str().c_str(), hist.nBins, hist.xRange.min, hist.xRange.max);

        for(const DataSource& source : category.sources) {
            TH1D* source_histogram = Convert(source, hist, "syncTree", isData ? "" : "");
            histogram->Add(source_histogram, source.scale_factor);
        }
        return histogram;
    }

    static TH1D* Convert(const DataSource& source, const HistogramDescriptor& hist, const std::string& treeName,
                         const std::string& weightBranchName = "")
    {
        TTree* tree = static_cast<TTree*>(source.file->Get(treeName.c_str()));;
        if(!tree)
            throw std::runtime_error("tree not found.");
        std::ostringstream name;
        name << source.file_name << "_" << hist.name << "_hist";
        std::ostringstream command;
        command << "(" << hist.name;
        if(weightBranchName.size())
            command << "*" << weightBranchName;
        command << ")>>+" << name.str();

        TH1D* histogram = new TH1D(name.str().c_str(), name.str().c_str(), hist.nBins, hist.xRange.min, hist.xRange.max);
        tree->Draw(command.str().c_str(), " njetspt20>=2 && (nbtag >=2 && bcsv_1 > 0.679 && bcsv_2 > 0.679) && (q_1*q_2==-1)");
        return histogram;
    }


    Print_Stack(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                const std::string& outputFileName, const std::string& _signalName, const std::string& _dataName)
        : inputPath(_inputPath), signalName(_signalName), dataName(_dataName), printer(outputFileName)
    {
        ReadSourceCfg(source_cfg);
        ReadHistCfg(hist_cfg);
    }

    void Run()
    {
        std::cout << "Opening Sources... " << std::endl;
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
            }
        }

        std::cout << "Printing Histograms ... " << std::endl;
        for (const HistogramDescriptor& hist : histograms){
            //std::cout << hist << std::endl;
            page.side.use_log_scaleY = hist.useLogY;
            page.side.xRange = hist.xRange;
            page.side.fit_range_x = false;
            page.title = hist.title;
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

            //for (const DataSource& source : sources){
            TH1D* data_histogram = nullptr;
            for(const DataCategory& category : categories) {
                const bool isData = category.name.find("DATA") != std::string::npos;
                TH1D* histogram = CreateHistogram(category, hist,isData);
                histogram->SetLineColor(category.color);
                std::string legend_option = "f";
                if (!isData){
                    histogram->SetFillColor(category.color);
                    stack->Add(histogram);
                }
                else {
                    legend_option = "p";
                    data_histogram = histogram;
                }
                leg->AddEntry(histogram, category.title.c_str(), legend_option.c_str());
            }
            printer.PrintStack(page, *stack, *data_histogram, *leg, *ll);
        }
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

private:
    std::string inputPath;
    std::string signalName, dataName;
    DataCategoryCollection categories;
    std::vector<HistogramDescriptor> histograms;
    Printer printer;
    root_ext::SingleSidedPage page;
};
