/*!
 * \class PileUp PileUp.h
 * \brief Class to calculate and apply PU reweighting.
 *
 * Based on PhysicsTools/Utilities/src/LumiReWeighting.cc
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2013-04-08 created
 */

#pragma once

#include "TH1D.h"
#include "../include/TreeExtractor.h"

class PileUp{
public:
    PileUp(const std::string& MC_File_name, const std::string& Data_File_name, const std::string& reweight_fileName,
                 const std::string& histName, const std::string& _mode, const std::string& _prefix = "none",
                        size_t _maxNumberOfEvents = 0)
        : treeExtractor(_prefix == "none" ? "" : _prefix, MC_File_name, false),
          outputFile(new TFile(reweight_fileName.c_str(),"RECREATE")),
          maxNumberOfEvents(_maxNumberOfEvents), mode(_mode)
    {

        TFile* Data_File = new TFile(Data_File_name.c_str(), "READ");
        outputFile->cd();
        Data_distr = static_cast<TH1D*>(Data_File->Get(histName.c_str())->Clone());
        outputFile->cd();
        Data_distr->Write();
        Data_distr->Scale( 1.0/ Data_distr->Integral() );
        Data_File->Close();
        outputFile->cd();
        nPU_MCdistr = new TH1D("MC_pileup", "MC nPU distribution", Data_distr->GetNbinsX(), 0, Data_distr->GetNbinsX());
    }


    void Run()
    {
        size_t n = 0;
        analysis::EventDescriptor event;
        for(; !maxNumberOfEvents || n < maxNumberOfEvents; ++n) {
            if(!treeExtractor.ExtractNext(event))
                break;
            const ntuple::Event& eventInfo = event.eventInfo();
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

        outputFile->cd();

        nPU_MCdistr->Write();
        nPU_MCdistr->Scale( 1.0 / nPU_MCdistr->Integral() );


        TH1D* weights = static_cast<TH1D*>(Data_distr->Clone("weights"));
        weights->Divide(nPU_MCdistr);
        weights->Write();

        outputFile->Close();
    }

private:
    analysis::TreeExtractor treeExtractor;
    TFile* outputFile;
    size_t maxNumberOfEvents;
    TH1D* Data_distr;
    TH1D* nPU_MCdistr;
    std::string mode;
};
