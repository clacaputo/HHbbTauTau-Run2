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
    Color_t color;
    TFile* file;
};

struct HistogramDescriptor {
    std::string name;
    std::string title;
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
    s << "Name: " << hist.name << ", Title: " << hist.title << ", xmin: " << hist.xmin << ", xmax: " << hist.xmax <<
         ", rebinFactor: " << hist.rebinFactor << ", useLog: " << hist.useLogY  ;
    return s;
}

class Print_Stack {
public:

    typedef root_ext::PdfPrinter Printer;

    Print_Stack(const std::string& source_cfg, const std::string& hist_cfg, const std::string& outputFileName) :
        printer(outputFileName)
    {
        ReadSourceCfg(source_cfg);
        ReadHistCfg(hist_cfg);
    }

    void Run()
    {

        std::cout << "Opening Sources... " << std::endl;
        for (const DataSource& source : sources){
            //std::cout << source << std::endl;
            source.file = new TFile(source.file_name.c_str(),"READ");
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
            THStack* stack = new THStack(hist.name.c_str(),hist.title.c_str());

            for (const DataSource& source : sources){
                TH1D* histogram = static_cast<TH1D*>(source.file->Get(hist.name.c_str()));
                if (!histogram){
                    std::ostringstream ss;
                    ss << "histogram '" << hist.name << "' not found in " << source.file_name;
                    throw std::runtime_error(ss.str());
                }
                histogram->Rebin(hist.rebinFactor);

            }
        }
    }

private:

    void ReadSourceCfg(const std::string& cfg_name)
    {
        std::ifstream cfg(cfg_name);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size()) continue;
            std::istringstream ss(cfgLine);
            DataSource source;
            ss >> source.name;
            ss >> source.file_name;
            ss >> source.scale_factor;
            ss >> source.color;
            sources.push_back(source);
          }
    }

    void ReadHistCfg(const std::string& cfg_name)
    {
        std::ifstream cfg(cfg_name);
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size()) continue;
            std::istringstream ss(cfgLine);
            HistogramDescriptor hist;
            ss >> hist.name;
            ss >> hist.title;
            ss >> hist.xRange.min;
            ss >> hist.xRange.max;
            ss >> hist.rebinFactor;
            ss >> hist.useLogY;
            histograms.push_back(hist);
          }
    }


private:
    std::vector<DataSource> sources;
    std::vector<HistogramDescriptor> histograms;
    Printer printer;
    root_ext::SingleSidedPage page;

};

