/*!
 * \file Print_Selection.C
 * \brief Print event selection histograms.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-24 created
 */

#include "../include/RootPrintToPdf.h"

class MyHistogramSource : public root_ext::SimpleHistogramSource<TH1D, Double_t> {
public:

protected:
    virtual void Prepare(TH1D* histogram, const std::string& display_name,
                         const PlotOptions& plot_options) const
    {
        SimpleHistogramSource::Prepare(histogram, display_name, plot_options);
    }

}; //1D plot


class Print_Selection {
public:
    typedef std::pair< std::string, std::string > FileTagPair;

    template<typename ...Args>
    Print_Selection(const std::string& _outputFileName, const Args& ...args):
       outputFileName(_outputFileName)
    {
        Initialize(args...);
    }

    void Run()
    {
        typedef root_ext::PdfPrinter Printer;
        MyHistogramSource source;

        for(const FileTagPair& fileTag : inputs) {
            TFile* file = new TFile(fileTag.first.c_str());
            source.Add(fileTag.second, file);
        }

        root_ext::SingleSidedPage page;
        page.side.histogram_name = "EventSelection_abs";
        page.title = page.side.histogram_title = "Event selection";
        page.side.axis_titleY = "N entries";
        page.side.use_log_scaleY = true;
        page.side.layout.legend_pad = root_ext::Box(0.6, 0.67, 0.87, 0.8);

        Printer printer(outputFileName);
        printer.Print(page, source);
    }

private:
    template<typename ...Args>
    void Initialize(const std::string& inputName, const Args& ...args)
    {
        const size_t split_index = inputName.find_first_of(':');
        const std::string fileName = inputName.substr(0, split_index);
        const std::string tagName = inputName.substr(split_index + 1);
        inputs.push_back(FileTagPair(fileName, tagName));
        Initialize(args...);
    }

    void Initialize() {}

private:

    std::string outputFileName;
    std::vector<FileTagPair> inputs;
};
