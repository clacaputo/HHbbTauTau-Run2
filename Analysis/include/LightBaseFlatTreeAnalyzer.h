/*!
 * \file LightBaseFlatTreeAnalyzer.h
 * \brief Definition of LightBaseFlatTreeAnalyzer class, the base class for separate studies on flat trees.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-08 created
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

#ifndef __APPLE__
#define override
#endif

#include <iostream>
#include <cmath>
#include <set>
#include <list>
#include <locale>
#include <numeric>
#include <algorithm>


#include <TColor.h>
#include <TLorentzVector.h>


#include "AnalysisBase/include/AnalyzerData.h"
#include "AnalysisBase/include/FlatEventInfo.h"
#include "AnalysisBase/include/AnalysisMath.h"
#include "AnalysisBase/include/exception.h"
#include "AnalysisBase/include/Particles.h"
#include "PrintTools/include/RootPrintToPdf.h"
//#include "KinFit.h"

#include "MVASelections/include/MvaReader.h"

#include "Htautau_Summer13.h"
#include "AnalysisCategories.h"

namespace analysis {

class LightFlatAnalyzerData : public root_ext::AnalyzerData {
public:
    LightFlatAnalyzerData(TFile& outputFile) : AnalyzerData(outputFile) {}

};

class LightBaseFlatTreeAnalyzer {
public:
    LightBaseFlatTreeAnalyzer(const std::string& inputFileName, const std::string& outputFileName)
        : inputFile(new TFile(inputFileName.c_str(),"READ")),
          outputFile(new TFile(outputFileName.c_str(), "RECREATE")),
          flatTree(new ntuple::FlatTree(*inputFile, "flatTree"))
    {
        TH1::SetDefaultSumw2();
    }

    virtual ~LightBaseFlatTreeAnalyzer() {}

    void Run()
    {
        using namespace cuts::Htautau_Summer13::btag;
        for(Long64_t current_entry = 0; current_entry < flatTree->GetEntries(); ++current_entry) {
            flatTree->GetEntry(current_entry);
            const ntuple::Flat& event = flatTree->data;
            const EventCategoryVector eventCategories = DetermineEventCategories(event.csv_Bjets, CSVM, CSVT);
            FlatEventInfo eventInfo(event, FlatEventInfo::BjetPair(0, 1));
            for (EventCategory eventCategory : eventCategories)
                AnalyzeEvent(eventInfo, eventCategory);
        }
    }



protected:
    virtual LightFlatAnalyzerData& GetAnaData() = 0;
    virtual void AnalyzeEvent(const FlatEventInfo& eventInfo, EventCategory eventCategory) = 0;

    TFile& GetOutputFile() { return *outputFile; }

    bool PassSelection(const FlatEventInfo& eventInfo, EventCategory /*category*/) const
    {
        using namespace cuts::Htautau_Summer13;
        const ntuple::Flat& event = *eventInfo.event;
        if(eventInfo.channel == Channel::MuTau)
            return event.byCombinedIsolationDeltaBetaCorrRaw3Hits_2 < MuTau::tauID::byCombinedIsolationDeltaBetaCorrRaw3Hits
                    && event.againstMuonTight_2
                    && event.mt_1 < MuTau::muonID::mt
                    && event.pfRelIso_1 < MuTau::muonID::pFRelIso
                    && event.q_1 * event.q_2 == -1;

        throw exception("Channel '") << eventInfo.channel << "' is not supported.";
    }

private:
    std::shared_ptr<TFile> inputFile, outputFile;
    std::shared_ptr<ntuple::FlatTree> flatTree;
};

} // namespace analysis
