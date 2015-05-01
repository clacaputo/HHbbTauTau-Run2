/*!
 * \class PileUp PileUp.C
 * \brief Class to calculate and apply PU reweighting.
 *
 * Based on PhysicsTools/Utilities/src/LumiReWeighting.cc
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \author Maria Teresa Grippo (University of Siena, INFN Pisa)
 * \date 2013-04-08 created
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

#include <TH1.h>

#include "AnalysisBase/include/TreeExtractor.h"
#include "AnalysisBase/include/RootExt.h"

class PileUp {
public:
    PileUp(const std::string& MC_File_name, const std::string& Data_File_name, const std::string& reweight_fileName,
                 const std::string& histName, const std::string& _mode, const std::string& _prefix = "none",
                        size_t _maxNumberOfEvents = 0)
        : treeExtractor(_prefix == "none" ? "" : _prefix, MC_File_name, false),
          outputFile(root_ext::CreateRootFile(reweight_fileName)),
          maxNumberOfEvents(_maxNumberOfEvents), mode(_mode)
    {

        auto Data_File = root_ext::OpenRootFile(Data_File_name);
        Data_distr = root_ext::ReadCloneObject(*Data_File, histName);
        root_ext::WriteObject(*Data_distr, outputFile.get());
        Data_distr->Scale( 1.0/ Data_distr->Integral() );
        Data_File->Close();
        nPU_MCdistr = new TH1D("MC_pileup", "MC nPU distribution", Data_distr->GetNbinsX(),
                               Data_distr->GetBinLowEdge(1),
                               Data_distr->GetBinLowEdge(Data_distr->GetNbinsX() + 1));
        nPU_MCdistr->SetDirectory(outputFile.get());
    }


    void Run()
    {
        size_t n = 0;
        analysis::EventDescriptor event;
        for(; !maxNumberOfEvents || n < maxNumberOfEvents; ++n) {
            if(!treeExtractor.ExtractNext(event))
                break;
            const ntuple::Event& eventInfo = event->eventInfo();
            for (unsigned n = 0; n < eventInfo.bunchCrossing.size(); ++n){
                if (eventInfo.bunchCrossing.at(n) == 0){
                    if(mode == "true")
                        nPU_MCdistr->Fill(eventInfo.trueNInt.at(n));
                    else if(mode == "observed")
                        nPU_MCdistr->Fill(eventInfo.nPU.at(n));
                    else
                        throw std::runtime_error("Unknown mode.");
                    break;
                }
            }
        }

        root_ext::WriteObject(*nPU_MCdistr);
        nPU_MCdistr->Scale( 1.0 / nPU_MCdistr->Integral() );

        TH1D* weights = root_ext::CloneObject(*Data_distr, "weights");
        weights->Divide(nPU_MCdistr);
        root_ext::WriteObject(*weights, outputFile.get());
    }

private:
    analysis::TreeExtractor treeExtractor;
    std::shared_ptr<TFile> outputFile;
    size_t maxNumberOfEvents;
    TH1D* Data_distr;
    TH1D* nPU_MCdistr;
    std::string mode;
};
