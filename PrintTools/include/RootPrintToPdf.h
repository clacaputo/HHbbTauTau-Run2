/*!
 * \file RootPrintToPdf.h
 * \brief Print ROOT histograms to PDF.
 * \author Konstantin Androsov <konstantin.androsov@gmail.com>
 */

#pragma once

#include <sstream>
#include <stdexcept>

#include <TROOT.h>
#include <TStyle.h>
#include <TFile.h>
#include <Rtypes.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

#include "RootPrintSource.h"

namespace root_ext {

class PdfPrinter {
public:
    PdfPrinter(const std::string& _output_file_name)
        : output_file_name(_output_file_name)
    {
        canvas = new TCanvas();
        canvas->Print((output_file_name + "[").c_str());
    }

    template<typename Source>
    void Print(const Page& page, const Source& source)
    {
        gROOT->SetStyle(page.layout.global_style.c_str());
        gStyle->SetOptStat(page.layout.stat_options);
        gStyle->SetOptFit(page.layout.fit_options);
        canvas->cd();

        canvas->SetTitle(page.title.c_str());
        if(page.layout.has_title) {
            TPaveLabel *title = Adapter::NewPaveLabel(page.layout.title_box, page.title);
            title->SetTextFont(page.layout.title_font);
            title->Draw();
        }

        Page::RegionCollection page_regions = page.Regions();
        for(Page::RegionCollection::const_iterator iter = page_regions.begin(); iter != page_regions.end(); ++iter)
        {
            canvas->cd();
            DrawHistograms(*(*iter), source);
        }

        canvas->Draw();
        std::ostringstream print_options;
        print_options << "Title:" << page.title;
        canvas->Print(output_file_name.c_str(), print_options.str().c_str());
    }

    ~PdfPrinter()
    {
        canvas->Print((output_file_name+"]").c_str());
    }

private:
    template<typename Source>
    void DrawHistograms(const PageSide& page_side, const Source& source)
    {
        typedef root_ext::HistogramPlotter<typename Source::Histogram, typename Source::ValueType> Plotter;
        typedef root_ext::HistogramFitter<typename Source::Histogram, typename Source::ValueType> Fitter;

        TPad* stat_pad = 0;
        if(page_side.layout.has_stat_pad) {
            stat_pad = Adapter::NewPad(page_side.layout.stat_pad);
            stat_pad->Draw();
        }

        TPad *pad = Adapter::NewPad(page_side.layout.main_pad);
        if(page_side.use_log_scaleX)
            pad->SetLogx();
        if(page_side.use_log_scaleY)
            pad->SetLogy();
        pad->Draw();
        pad->cd();

        Plotter plotter(page_side.histogram_title, page_side.axis_titleX, page_side.axis_titleY);
        for(unsigned n = 0; n < source.Size(); ++n)
        {
            const typename Plotter::Entry entry = source.Get(n, page_side.histogram_name);
            plotter.Add(entry);
        }

        Fitter::SetRanges(plotter.Histograms(), page_side.fit_range_x, page_side.fit_range_y, page_side.xRange,
                          page_side.yRange, page_side.use_log_scaleY);
        plotter.Superpose(pad, stat_pad, page_side.layout.has_legend, page_side.layout.legend_pad,
                          page_side.draw_options);
    }

private:
    TCanvas* canvas;
    std::string output_file_name;
};

} // namespace root_ext