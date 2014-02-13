/*!
 * \file Htautau_Summer13 Htautau_Summer13.h
 * \brief Higgs in tautau recommended baseline selection cuts for Summer13
 *        see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-13 created
 */

#pragma once

namespace cuts
{
    const double DeltaR_signalLeptons = 0.5;

    namespace tauID{
        const double pt = 20;
        const double eta = 2.3;
        const int decayModeFinding = 1;
        const int byLooseIsolationDeltaBetaCorr = 1;
        const int againstMuonTight = 1;
        const int againstElectronLoose = 1;
    }

    namespace muonID{
        const double pt = 20;
        const double eta = 2.1;
        const bool isTrackerMuon = true;
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

    namespace electronID{
        const double pt = 10;
        const double eta_high = 2.5;
        const double eta_CrackVeto_low = 1.4442;
        const double eta_CrackVeto_high = 1.566;
        const double dz_pv = 0.2;
        const int missingHits = 1;
        const int hasMatchedConv = 1;
        const double dB = 0.045;
        const double ref_pt = 20;
        const unsigned nCategories = 3;
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
}
