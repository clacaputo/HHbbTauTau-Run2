/*!
 * \file Print_Selection.C
 * \brief Print event selection histograms.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2014-02-24 created
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

#include <iostream>

#include "AnalysisBase/include/RootExt.h"
#include "../include/RootPrintToPdf.h"

class Print_Selection {
public:
    typedef std::pair< std::string, std::string > FileTagPair;
    typedef root_ext::PdfPrinter Printer;
    typedef root_ext::SimpleHistogramSource<TH1D, Double_t> MyHistogramSource;

    template<typename ...Args>
    Print_Selection(const std::string& outputFileName, const Args& ...args):
       printer(outputFileName)
    {
        Initialize(args...);
        for(const FileTagPair& fileTag : inputs) {
            auto file = root_ext::OpenRootFile(fileTag.first);
            source.Add(fileTag.second, file);
        }
    }

    void Run()
    {
        page.side.use_log_scaleY = true;
        page.side.layout.has_legend = false;

        PrintAll("Selection_event", "Event selection");
        PrintAll("Selection_muons", "Muon selection");
        PrintAll("Selection_taus", "Tau selection");
        PrintAll("Selection_electrons_bkg", "Bkg Electron selection");
        PrintAll("Selection_bjets", "b-jet selection");


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

    void PrintAll(const std::string& name, const std::string& title, std::string second_suffix = "")
    {
        try {
            Print(name, title, AddSuffix("", second_suffix), "Number of entries");
            Print(name, title, AddSuffix("effRel", second_suffix), "Relative efficiency");
            Print(name, title, AddSuffix("effAbs", second_suffix), "Absolute efficiency");
        } catch(std::runtime_error& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
        }
    }

    void Print(const std::string& name, const std::string& title,
               std::string name_suffix = "", std::string title_suffix = "")
    {
        page.side.histogram_name = AddSuffix(name, name_suffix);
        page.title = page.side.histogram_title = AddSuffix(title, title_suffix, ". ");
        printer.Print(page, source);
    }

    static std::string AddSuffix(const std::string& name, const std::string& suffix, std::string separator = "_")
    {
        std::ostringstream full_name;
        full_name << name;
        if(suffix.size())
            full_name << separator << suffix;
        return full_name.str();
    }

private:
    std::vector<FileTagPair> inputs;
    MyHistogramSource source;
    root_ext::SingleSidedPage page;
    Printer printer;
};
