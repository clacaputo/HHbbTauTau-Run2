/*!
 * \file H_BaselineSync.C
 * \brief Analyzer which recalls three ananlysis for Htautau using baseline selection.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-06-30 created
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

#include "../HHbbTauTau/Analysis/include/TreeExtractor.h"
#include "../HHbbTauTau/Analysis/include/Config.h"
#include "../HHbbTauTau/Analysis/source/HmutauBaseline_sync.C"
#include "../HHbbTauTau/Analysis/source/HetauBaseline_sync.C"
#include "../HHbbTauTau/Analysis/source/HtautauBaseline_sync.C"



class H_BaselineSync {
public:
    H_BaselineSync(const std::string& inputFileName, const std::string& outputMuTauFile, const std::string& outputETauFile,
                   const std::string& outputTauTauFile, const std::string& configFileName,
                   const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : config(configFileName), timer(config.ReportInterval()), maxNumberOfEvents(_maxNumberOfEvents),
          treeExtractor(_prefix == "none" ? "" : _prefix, inputFileName, config.extractMCtruth()),
          HmutauAnalyzer(inputFileName, outputMuTauFile, configFileName, "external", _maxNumberOfEvents),
          HetauAnalyzer(inputFileName, outputETauFile, configFileName, "external", _maxNumberOfEvents),
          HtautauAnalyzer(inputFileName, outputTauTauFile, configFileName, "external", _maxNumberOfEvents)
    { }

    virtual void Run()
    {
        size_t n = 0;
        auto _event = std::shared_ptr<analysis::EventDescriptor>(new analysis::EventDescriptor());
        for(; ( !maxNumberOfEvents || n < maxNumberOfEvents ) && treeExtractor.ExtractNext(*_event); ++n) {
            timer.Report(n);
            if(config.RunSingleEvent() && _event->eventId().eventId != config.SingleEventId()) continue;
            try {
                HmutauAnalyzer.ProcessEvent(_event);
                HetauAnalyzer.ProcessEvent(_event);
                HtautauAnalyzer.ProcessEvent(_event);
            } catch(cuts::cut_failed&){}
            if(config.RunSingleEvent()) break;
        }
        timer.Report(n, true);
    }


private:
    analysis::Config config;
    tools::Timer timer;
    size_t maxNumberOfEvents;
    std::shared_ptr<const analysis::EventDescriptor> event;
    analysis::TreeExtractor treeExtractor;
    HmutauBaseline_sync HmutauAnalyzer;
    HetauBaseline_sync HetauAnalyzer;
    HtautauBaseline_sync HtautauAnalyzer;
};

namespace make_tools {
template<>
struct Factory<H_BaselineSync> {
    static H_BaselineSync* Make(int argc, char *argv[])
    {
        std::cout << "Command line: ";
        for(int n = 0; n < argc; ++n)
            std::cout << argv[n] << " ";
        std::cout << std::endl;
        if(argc < 6 || argc > 8)
            throw std::runtime_error("Invalid number of command line arguments.");

        int n = 0;
        const std::string inputFileName = argv[++n];
        const std::string outputMuTauFile = argv[++n];
        const std::string outputETauFile = argv[++n];
        const std::string outputTauTauFile = argv[++n];
        const std::string configFileName = argv[++n];
        if(argc <= n)
            return new H_BaselineSync(inputFileName, outputMuTauFile, outputETauFile, outputTauTauFile, configFileName);

        const std::string prefix = argv[++n];
        if(argc <= n)
            return new H_BaselineSync(inputFileName, outputMuTauFile, outputETauFile, outputTauTauFile,
                                      configFileName, prefix);

        char c;
        size_t maxNumberOfEvents;
        std::istringstream ss_nEvents(argv[++n]);
        ss_nEvents >> c >> maxNumberOfEvents;
        if(c != '@')
            throw std::runtime_error("Bad command line format.");

        return new H_BaselineSync(inputFileName, outputMuTauFile, outputETauFile, outputTauTauFile,
                                  configFileName, prefix, maxNumberOfEvents);
    }
};
} // make_tools
