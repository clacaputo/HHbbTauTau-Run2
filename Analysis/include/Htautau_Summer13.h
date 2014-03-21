/*!
 * \file Htautau_Summer13.h
 * \brief Higgs in tautau recommended baseline selection cuts for Summer13
 *        see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-13 created
 */

#pragma once

namespace cuts {
namespace Htautau_Summer13 {

const double DeltaR_signalLeptons = 0.5;

namespace tauID {
    namespace ETau {
        const double pt = 20; // twiki TauIDRecommendation
        const double eta = 2.3; // twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double decayModeFinding = 0.5; // AN-10-82
        const double againstMuonLoose = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronMediumMVA5 = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                      // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                      // MVA3 is recommended, but it does not exists any more
        const double LooseCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                      // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
    namespace MuTau {
        const double pt = 20; // twiki TauIDRecommendation
        const double eta = 2.3; // twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double decayModeFinding = 0.5; // AN-10-82
        const double againstMuonTight = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double LooseCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                      // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
    namespace TauTau {
        const double pt = 45; // AN-13-189
        const double eta = 2.1; // AN-13-189
        const double decayModeFinding = 0.5; // AN-10-82
        const double againstMuonLoose = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLooseMVA5 = 0.5; // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                     // recommended only for sub-leading tau
                                                     // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                     // MVA3 is recommended, but it does not exists any more
        const double MediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
}

namespace muonID {
    namespace MuTau {
        const double pt = 20; // twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double eta = 2.1; // twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const bool isTrackerMuon = true; //
        const bool isGlobalMuonPromptTight = true;
        const bool isPFMuon = true;
        const int nChambers = 1;
        const int nMatched_Stations = 1;
        const int trackerLayersWithMeasurement = 5;
        const int pixHits = 0;
        const double globalChiSquare = 10;
        const double dB = 0.045;
        const double pFRelIso = 0.1;
    }
}

namespace electronID{
    const double pt = 10;
    const double eta_high = 2.5;
    const double eta_CrackVeto_low = 1.4442;
    const double eta_CrackVeto_high = 1.566;
    const double dz_pv = 0.2;
    const int missingHits = 1;
    const float hasMatchedConv = 0.5;
    const double dB = 0.045;
    const double ref_pt = 20;
    const double scEta_min[2] = {0.8, 1.479};
    const double MVApogNonTrig[2][3] = {{0.925, 0.915, 0.965},{0.905,0.955, 0.975}};
}

namespace btag{
    const double pt = 30;
    const double eta = 2.4;
    const double CSVL = 0.244; // loose
    const double CSVM = 0.679; //medium
    const double CSV = CSVM; // recommended
}

} // Htautau_Summer13
} // cuts
