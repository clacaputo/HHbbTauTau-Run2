/*
 * \file HHbbTauTau_analyzer.C
 * \brief Analysis HHbbTauTau.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-12 created
 */

#include <vector>
using namespace std;

#define ExtractTree_cxx
#include "../include/ExtractTree.h"
#include "../include/AnalyzerData.h"
#include "../include/Particles.h"
#include "TH2.h"
#include <TStyle.h>
#include <TCanvas.h>
#include <boost/shared_ptr.hpp>

void ExtractTree::Loop(){} //just declared


class SignalAnalyzerData : public root_ext::AnalyzerData {
public:
    SignalAnalyzerData(const std::string& outputFileName) : AnalyzerData(outputFileName) {}
    TH1D_ENTRY(Pt_tau, 100, 0, 100)
    TH1D_ENTRY(N_realtau, 100, 0, 100)

};

class HHbbTauTau_analyzer
{
public:
    HHbbTauTau_analyzer(const std::string& inputFileName, const std::string& outputFileName,
                        Long64_t _maxNumberOfEvents = 0)
        : anaData(outputFileName), maxNumberOfEvents(_maxNumberOfEvents)
    {
        TFile* inputFile = new TFile(inputFileName.c_str(),"READ");
        TTree* inputTree = dynamic_cast<TTree*> (inputFile->Get("treeCreator/vhtree"));
        eventTree = boost::shared_ptr<ExtractTree>(new ExtractTree(inputTree));
        anaData.getOutputFile().cd();
        std::cout << "starting analyzer" << std::endl;
    }

    void Run()
    {
        if (eventTree->fChain == nullptr) return;

        const Long64_t nentries = maxNumberOfEvents ? std::min(eventTree->fChain->GetEntriesFast(),maxNumberOfEvents)
                                                    : eventTree->fChain->GetEntriesFast();
        for (Long64_t jentry=0; jentry<nentries;jentry++) {
           const Long64_t ientry = eventTree->LoadTree(jentry);
           if (ientry < 0)
               throw std::runtime_error("cannot read entry");
           eventTree->fChain->GetEntry(jentry);
           ProcessEvent();
        }
    }

private:

    void ProcessEvent()
    {
        for (Int_t n = 0; n < eventTree->nTau; ++n){
            anaData.Pt_tau().Fill(eventTree->Tau_pt[n]);
        }

        unsigned tau_counter = 0;
        for (Int_t n = 0; n < eventTree->nGenParticle; ++n){
            if (std::abs(eventTree->GenParticle_pdgId[n]) == particles::tau.RawCode()
                    && eventTree->GenParticle_status[n] == particles::Decayed_or_fragmented)
                ++tau_counter;
        }
        anaData.N_realtau().Fill(tau_counter);
    }


private:
    boost::shared_ptr<ExtractTree> eventTree;
    SignalAnalyzerData anaData;
    Long64_t maxNumberOfEvents;
};

