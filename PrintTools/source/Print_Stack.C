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
};

struct HistogramDescriptor {
    std::string name;
    std::string title;
    double xmin, xmax;
    unsigned nbins;
    bool useLogY;
};


std::ostream& operator<<(std::ostream& s, const DataSource& source){
    s << "Name: " << source.name << ", File: " << source.file_name << ", SC: " << source.scale_factor <<
         ", color: " << source.color ;
    return s;
}

std::ostream& operator<<(std::ostream& s, const HistogramDescriptor& hist){
    s << "Name: " << hist.name << ", Title: " << hist.title << ", xmin: " << hist.xmin << ", xmax: " << hist.xmax <<
         ", nbins: " << hist.nbins << ", useLog: " << hist.useLogY  ;
    return s;
}

class Print_Stack {
public:
    typedef std::pair< std::string, std::string > FileTagPair;
    typedef root_ext::PdfPrinter Printer;
//    typedef root_ext::SimpleHistogramSource<TH1D, Double_t> MyHistogramSource;


    Print_Stack(const std::string& source_cfg, const std::string& hist_cfg)
    {
        ReadSourceCfg(source_cfg);
        ReadHistCfg(hist_cfg);
    }

    void Run()
    {
        std::cout << "Sources: " << std::endl;
        for (const auto& source : sources){
            std::cout << source << std::endl;
        }
        std::cout << "Histograms: " << std::endl;
        for (const auto& hist : histograms){
            std::cout << hist << std::endl;
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
            ss >> hist.xmin;
            ss >> hist.xmax;
            ss >> hist.nbins;
            ss >> hist.useLogY;
            histograms.push_back(hist);
          }
    }


private:
    std::vector<DataSource> sources;
    std::vector<HistogramDescriptor> histograms;

};

