/*!
 * \file Print_Stack.C
 * \brief Print stack with specified name superimposing several files.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-04-03 created
 */

#include <TTree.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "../include/RootPrintToPdf.h"

struct DataSource {
    std::string name;
    std::string file_name;
    double scale_factor;
    EColor color;
    TFile* file;
};

struct HistogramDescriptor {
    std::string name;
    std::string title;
    std::string Xaxis_title;
    std::string Yaxis_title;
    root_ext::Range xRange;
    unsigned rebinFactor;
    bool useLogY;
};


std::ostream& operator<<(std::ostream& s, const DataSource& source){
    s << "Name: " << source.name << ", File: " << source.file_name << ", SC: " << source.scale_factor <<
         ", color: " << source.color ;
    return s;
}

std::ostream& operator<<(std::ostream& s, const HistogramDescriptor& hist){
    s << "Name: " << hist.name << ", Title: " << hist.title << ", xmin: " << hist.xRange.min << ", xmax: " <<
         hist.xRange.max << ", rebinFactor: " << hist.rebinFactor << ", useLog: " << hist.useLogY  ;
    return s;
}

std::istream&operator>>(std::istream& s, EColor& color){
    const std::map<std::string, EColor> colorMapName = {{"white",kWhite}, {"black",kBlack}, {"gray",kGray},
                                                       {"red",kRed}, {"green",kGreen}, {"blue",kBlue},
                                                       {"yellow",kYellow}, {"magenta",kMagenta}, {"cyan",kCyan},
                                                       {"orange",kOrange}, {"spring",kSpring}, {"teal",kTeal},
                                                       {"azure",kAzure}, {"violet",kViolet},{"pink",kPink}};
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

class Print_Stack {
public:
    typedef root_ext::PdfPrinter Printer;

    Print_Stack(const std::string& source_cfg, const std::string& hist_cfg, const std::string& _inputPath,
                const std::string& outputFileName)
        : inputPath(_inputPath), printer(outputFileName)
    {
        ReadSourceCfg(source_cfg);
        ReadHistCfg(hist_cfg);
    }

    void Run()
    {
        std::cout << "Opening Sources... " << std::endl;
        for (DataSource& source : sources){
            //std::cout << source << std::endl;
            const std::string fullFileName = inputPath + "/" + source.file_name;
            source.file = new TFile(fullFileName.c_str(), "READ");
            if(source.file->IsZombie()) {
                std::ostringstream ss;
                ss << "Input file '" << source.file_name << "' not found.";
                throw std::runtime_error(ss.str());
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
            std::shared_ptr<THStack> stack = std::shared_ptr<THStack>(new THStack(hist.name.c_str(),hist.title.c_str()));

            TLegend* leg = new TLegend( 0.57, 0.65, 0.77, 0.92);
            leg->SetFillColor(0);
            leg->SetTextSize(0.035);
            TString lumist="19.6 fb^{-1}";
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
            for (unsigned n = 0; n<sources.size(); ++n){
                const DataSource& source = sources.at(n);
                TH1D* histogram = static_cast<TH1D*>(source.file->Get(hist.name.c_str()));
                if (!histogram){
                    std::ostringstream ss;
                    ss << "histogram '" << hist.name << "' not found in " << source.file_name;
                    throw std::runtime_error(ss.str());
                }
                histogram->Rebin(hist.rebinFactor);
                histogram->Scale(source.scale_factor);
                histogram->SetLineColor(source.color);
                std::string legend_option = "f";
                if (n < sources.size() - 1){
                    histogram->SetFillColor(source.color);
                    stack->Add(histogram);
                }
                else {
                    legend_option = "p";
                    data_histogram = histogram;
                }
                leg->AddEntry(histogram,source.name.c_str(),legend_option.c_str());
            }
            printer.PrintStack(page, *stack, *data_histogram, *leg, *ll);
        }
    }

private:
    void ReadSourceCfg(const std::string& cfg_name)
    {
        std::ifstream cfg(cfg_name);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            std::istringstream ss(cfgLine);
            DataSource source;
            ss >> source.name;
            ss >> source.file_name;
            ss >> source.scale_factor;
            ss >> source.color;
            sources.push_back(source);
          }
        if (sources.size() < 2)
            throw std::runtime_error("not enough sources");
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
            ss >> hist.useLogY;
            histograms.push_back(hist);
          }
    }

private:
    std::string inputPath;
    std::vector<DataSource> sources;
    std::vector<HistogramDescriptor> histograms;
    Printer printer;
    root_ext::SingleSidedPage page;
};
