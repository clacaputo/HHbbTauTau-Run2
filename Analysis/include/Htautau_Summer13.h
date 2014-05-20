/*!
 * \file Htautau_Summer13.h
 * \brief Higgs in tautau recommended baseline selection cuts for Summer13
 *        see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorkingSummer2013
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-13 created
 */

#pragma once
#include <string>
#include <vector>
#include <map>
#include <TMath.h>

#include "Tools.h"
#include "TreeProduction/interface/Tau.h"

namespace cuts {
namespace Htautau_Summer13 {
// AN-2013/178 H->etau,mutau
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/analyzers/DiLeptonAnalyzer.py
const double DeltaR_betweenSignalObjects = 0.5; // >

namespace MuTau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const std::vector<std::string> hltPaths =
            { "HLT_IsoMu17_eta2p1_LooseIsoPFTau20", "HLT_IsoMu18_eta2p1_LooseIsoPFTau20" };
    }

    namespace muonID {
        const double pt = 20; // > twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double eta = 2.1; // < twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state

//        const bool isTightMuon = true; // = HiggsToTauTauWorkingSummer2013#Muon_ID
        //def of isTightMuon: twiki SWGuideMuonId#Tight_Muon
        const bool isGlobalMuonPromptTight = true;
              // = https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/MuonReco/src/MuonSelectors.cc#567 and 590 definition
              // = isGlobalMuon && normalizedChi2<10 && numberOfValidMuonHits > 0
        const bool isPFMuon = true; // = def of isTightMuon
        const int nMatched_Stations = 1; // > def of isTightMuon
        const int pixHits = 0; // > def of isTightMuon
        const int trackerLayersWithMeasurement = 5; // > def of isTightMuon
        const double dz = 0.2; // < twiki HiggsToTauTauWorkingSummer2013#Muon_ID
        const double dB = 0.045; // < twiki HiggsToTauTauWorkingSummer2013#Muon_ID

        const double pFRelIso = 0.1; // < twiki SWGuideMuonId#Muon_Isolation_AN1
    }

    namespace tauID {
        const double pt = 20; // > twiki TauIDRecommendation
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonTight = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
   }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
    namespace ZmumuVeto {
        const double pt = 15; // >
        const double eta = 2.4; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // < same definition as db
        const bool isGlobalMuon = true; // = already applied at Tree selection level as skim - MuonBlock.cc
        const bool isTrackerMuon = true; // =
        const bool isPFMuon = true; // =
        const double pfRelIso = 0.3; // <
        const double deltaR = 0.15; // >
        const bool haveOppositeCharge = true; // =
    }
}

namespace ETau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const std::vector<std::string> hltPaths =
            { "HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20",
              "HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20" };
    }

    namespace electronID{
        const double pt = 24; // >  HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double eta_high = 2.1; // <  HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double dz = 0.2; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const int missingHits = 1; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const bool hasMatchedConversion = false; // =  HiggsToTauTauWorkingSummer2013#Electron_ID
        const double dB = 0.045; // <  HiggsToTauTauWorkingSummer2013#Electron_ID
        const double pFRelIso = 0.1; // < twiki HiggsToTauTauWorkingSummer2013#Electron_Muon_Isolation
        const double scEta_min[2] = { 0.8, 1.479 }; // tight HiggsToTauTauWorkingSummer2013#Electron_ID
        const double MVApogNonTrig[3] = { 0.925, 0.975, 0.985 }; // tight HiggsToTauTauWorkingSummer2013#Electron_ID
    }

    namespace tauID {
        const double pt = 20; // > twiki TauIDRecommendation
        const double eta = 2.3; // < twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronMediumMVA3 = 0.5; //  > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                      // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                      // MVA3 is recommended, but it does not exists any more
        const double byCombinedIsolationDeltaBetaCorrRaw3Hits = 1.5;
                                                      // GeV < twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation

        // > https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/physicsobjects/Tau.py#L62
        const std::vector<double> againstElectronMediumMVA3_customValues = {
            0.933,0.921,0.944,0.945,0.918,0.941,0.981,0.943,0.956,0.947,0.951,0.95,0.897,0.958,0.955,0.942
        };
    }

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
    namespace ZeeVeto {
        const double pt = 15; // >
        const double eta = 2.5; // <
        const double dz = 0.2; // <
        const double d0 = 0.045; // <
        const double pfRelIso = 0.3; // <
        const double deltaR = 0.15; // >
        const bool haveOppositeCharge = true; // =

        // WP95 definition for Veto twiki EgammaCutBasedIdentification
        const double barrel_eta_high = 1.479; // <=
        const double endcap_eta_low = 1.479; // >
        const double endcap_eta_high = 2.5; // <
        const size_t barrel_index = 0, endcap_index = 1;
        const double sigma_ieta_ieta[] = { 0.01, 0.03 }; // <
        const double delta_eta[] = { 0.007, 0.01 }; // <
        const double delta_phi[] = { 0.8, 0.7 }; // <
        const double HoverE[] = { 0.15, std::numeric_limits<double>::max() }; // <
        const double dZ_vtx[] = { 0.2, 0.2 }; // <
    }
}

namespace TauTau {
    namespace trigger {
        // twiki HiggsToTauTauWorkingSummer2013#Tau_Tau_Final_state
        const std::map< std::string, bool > hltPathsMap =
        {{"HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30",true},
         {"HLT_DoubleMediumIsoPFTau35_Trk5_eta2p1", false},
         {"HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1",false} };

        const std::vector<std::string> hltPaths = tools::collect_map_keys(hltPathsMap);
    }

    namespace tauID {
        const double pt = 45; // > AN-2013/189 H->tautau full hadronic
        const double eta = 2.1; // < AN-2013/189 H->tautau full hadronic
        const double decayModeFinding = 0.5; // > AN-2010/082 Z->tautau
        const double againstMuonLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLoose = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
        const double againstElectronLooseMVA3 = 0.5; // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
                                                     // recommended only for sub-leading pt tau
                                                     // twiki SWGuidePFTauID#Tau_ID_2014_preparation_for_AN1
                                                     // MVA3 is recommended, but it does not exists any more for new tauID
        const double byMediumCombinedIsolationDeltaBetaCorr3Hits = 0.5;
                                                     // > twiki HiggsToTauTauWorkingSummer2013#Tau_ID_Isolation
    }
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/RootTools/python/physicsobjects/HTauTauElectron.py
namespace electronVeto {
    const double pt = 10; // >
    const double eta_high = 2.5; // <
    const double dz = 0.2; // <
    const double d0 = 0.045; // <
    const double pFRelIso = 0.3; // <
    const double ref_pt = 20; // twiki HiggsToTauTauWorkingSummer2013#Electron_ID
    const double scEta_min[2] = {0.8, 1.479}; // loose HiggsToTauTauWorkingSummer2013#Electron_ID
    const double MVApogNonTrig[2][3] = {{0.925, 0.915, 0.965},{0.905,0.955, 0.975}};
                                              // loose HiggsToTauTauWorkingSummer2013#Electron_ID
}

// AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Electron_Tau_Final_state
// twiki HiggsToTauTauWorkingSummer2013#Muon_Tau_Final_state
// https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/TauTauAnalyzer.py
namespace muonVeto {
    const double pt = 10; // >
    const double eta = 2.4; // <
    const double dz = 0.2; // <
    const double d0 = 0.045; // <

    //def of isTightMuon: twiki SWGuideMuonId#Tight_Muon
    const bool isGlobalMuonPromptTight = true;
          // = https://cmssdt.cern.ch/SDT/lxr/source/DataFormats/MuonReco/src/MuonSelectors.cc#567 and 590 definition
          // = isGlobalMuon && normalizedChi2<10 && numberOfValidMuonHits > 0
    const bool isPFMuon = true; // = def of isTightMuon
    const int nMatched_Stations = 1; // > def of isTightMuon
    const int pixHits = 0; // > def of isTightMuon
    const int trackerLayersWithMeasurement = 5; // > def of isTightMuon

    const double pfRelIso = 0.3; // <
}

namespace jetID {
    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 30; // >
    const double eta = 4.7; // <
    const double puLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >

    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const double pfLooseID = true; // =
}

namespace btag {
    // twiki BTagPerformanceOP#B_tagging_Operating_Points_for_5
    const double CSVL = 0.244; // > loose
    const double CSVM = 0.679; // > medium
    const double CSVT = 0.898; // > tight

    // AN-2013/188 H->tautau physics objects && twiki HiggsToTauTauWorkingSummer2013#Jets
    const double pt = 20; // >
    const double eta = 2.4; // <
    const double CSV = CSVM; // >
    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/python/proto/analyzers/VBFAnalyzer.py
    const double puLooseID = true; // =
    const double pfLooseID = true; // =
    const double deltaR_signalObjects = 0.5; // >
}

// AN-2013/188 H->tautau physics objects
namespace vertex {
    const double ndf = 4; // >
    const double z = 24.0; // < cm
    const double r = 2.0; // < cm
    const bool chooseHighestSumPt2 = true; // =
}

namespace tauCorrections {

    const double DecayModeWeight = 0.88; // = HiggsToTauTauWorkingSummer2013#TauTau_scale_factors
                                         // for 1-prong no pi 0 taus

    const double deltaR = 0.3; // < Updated to be compatible with H->tautau code

    // For taus that matched MC truth.
    // Original corrections from HiggsToTauTauWorkingSummer2013. Updated to be compatible with H->tautau code.
    inline double MomentumScaleFactor(bool hasMCmatch, double pt, ntuple::tau_id::hadronicDecayMode decayMode)
    {
        if(!hasMCmatch) return 1.0;
        if(decayMode == ntuple::tau_id::kOneProng1PiZero)
            //return 1.025 + 0.001 * std::min(std::max(pt - 45.0, 0.0), 10.0);
            return 1.012;
        if(decayMode == ntuple::tau_id::kThreeProng0PiZero)
            //scaleFactor = 1.012 + 0.001 * std::min(std::max(pt - 32.0, 0.0), 18.0);
            return 1.012;
        return 1.0;
    }
}

//namespace triggerEfficiencyCorrections {
//    // https://github.com/rmanzoni/HTT/blob/master/CMGTools/H2TauTau/interface/TriggerEfficiency.h
//    inline double eff2012IsoTau19fb(double pt, double eta)
//    {
//        return ( 808.411 * ( 0.764166 * 0.5 * (TMath::Erf((pt-33.2236)/2./0.97289 /sqrt(pt))+1.))
//        + 4428.0 * ( 0.802387 * 0.5 * (TMath::Erf((pt-38.0971)/2./0.82842 /sqrt(pt))+1.))
//        + 1783.003 * ( 0.818051 * 0.5 * (TMath::Erf((pt-37.3669)/2./0.74847 /sqrt(pt))+1.))
//        + 5109.155 * ( 0.796086 * 0.5 * (TMath::Erf((pt-37.3302)/2./0.757558/sqrt(pt))+1.))
//        + 4131. * ( 0.828182 * 0.5 * (TMath::Erf((pt-37.6596)/2./0.830682/sqrt(pt))+1.))
//        + 3143. * ( 0.833004 * 0.5 * (TMath::Erf((pt-37.634 )/2./0.777843/sqrt(pt))+1.)) )
//        /(808.411+4428.0+1783.003+5109.155+4131+3143);
//    }

//    double efficiency(double m, double m0, double sigma, double alpha, double n, double norm)
//     {
//        static const double sqrtPiOver2 = 1.2533141373;
//        static const double sqrt2 = 1.4142135624;
//        double sig = fabs(sigma);
//        double t = (m - m0)/sig;
//        if(alpha < 0)
//            t = -t;
//        double absAlpha = fabs(alpha/sig);
//        double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
//        double b = absAlpha - n/absAlpha;
//        double ApproxErf;
//        double arg = absAlpha / sqrt2;
//        if (arg > 5.) ApproxErf = 1;
//        else if (arg < -5.) ApproxErf = -1;
//        else ApproxErf = TMath::Erf(arg);

//        double leftArea = (1 + ApproxErf) * sqrtPiOver2;
//        double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
//        double area = leftArea + rightArea;
//        if( t <= absAlpha ) {
//            arg = t / sqrt2;
//            if(arg > 5.) ApproxErf = 1;
//            else if (arg < -5.) ApproxErf = -1;
//            else ApproxErf = TMath::Erf(arg);
//            return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
//        }
//        return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) - 1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
//    }

//    ///////////////////////////////trigger eff_Data//////////////////////////////

//    Double_t TauLegB18[5]= {18.604910,    0.276042,    0.137039,    2.698437,    0.940721};
//    Double_t TauLegE18[5]= {18.701715,    0.216523,    0.148111,    2.245081,    0.895320};

//      double MuLegEtaM3[] = {15.9977,    7.64004e-05,     6.4951e-08,     1.57403,     0.865325};
//      double MuLegEtaM2[] = {17.3974,     0.804001,    1.47145,     1.24295,     0.928198};
//      double MuLegEtaM1[] = {16.4307,     0.226312,     0.265553,     1.55756,     0.974462};

//      double MuLegEtaP1[] = {17.313,     0.662731,     1.3412,     1.05778,     1.26624};
//      double MuLegEtaP2[] = {16.9966,     0.550532,     0.807863,     1.55402,     0.885134};
//      double MuLegEtaP3[] = {15.9962,     0.000106195,     4.95058e-08,     1.9991,     0.851294} ;

//      ///////////////////////////////trigger eff_mc//////////////////////////////

//      Double_t TauLegB18_MC[5]= {18.537441,    1.385790,    3.102076,    1.002486,    6.275127};
//      Double_t TauLegE18_MC[5]= {18.393366,    1.526254,    2.021678,    124.741631,    0.894280};


//      double MuLegEtaM3_mc[] = {16.0051, 2.45144e-05, 4.3335e-09, 1.66134, 0.87045};
//      double MuLegEtaM2_mc[] = {17.3135, 0.747636, 1.21803, 1.40611, 0.934983};
//      double MuLegEtaM1_mc[] = {15.9556, 0.0236127, 0.00589832, 1.75409, 0.981338};

//      double MuLegEtaP1_mc[] = {15.9289, 0.0271317, 0.00448573, 1.92101, 0.978625};
//      double MuLegEtaP2_mc[] = {16.5678, 0.328333, 0.354533, 1.67085, 0.916992};
//      double MuLegEtaP3_mc[] = {15.997, 7.90069e-05, 4.40036e-08, 1.66272, 0.884502};

//      ///////////////////////////////Mu ID & ISO//////////////////////////////

//      double MuID_eta1[]={0.9818,0.9852}; //1st bin,2nd bin
//      double MuID_eta2[]={0.9829,0.9852}; //1st bin,2nd bin
//      double MuID_eta3[]={0.9869,0.9884}; //1st bin,2nd bin

//      double MuIso_eta1[]={0.9494,0.9883}; //1st bin,2nd bin
//      double MuIso_eta2[]={0.9835 ,0.9937}; //1st bin,2nd bin

//      double MuIso_eta3[]={0.9923 ,0.9996}; //1st bin,2nd bin

//      ///////////////////////////////Mu ID & ISO//////////////////////////////
//    double weight_HLT;
//      int runNumber = iEvent.id().run();
//      double  weightTotal;
//      std::vector<double> weight_MuLeg18, weight_MuLeg18DATA, weight_MuLeg18MC;
//      std::vector<double> weight_MuLeg17;
//      double weight_MuID;
//      double weight_MuIso;

//             if(etaMu <-1.2){weight_MuLeg18DATA.push_back(efficiency(ptMu, MuLegEtaM3[0], MuLegEtaM3[1], MuLegEtaM3[2], MuLegEtaM3[3], MuLegEtaM3[4])); }
//             if(-0.8<etaMu <-1.2){weight_MuLeg18DATA.push_back(efficiency(ptMu,MuLegEtaM2[0],MuLegEtaM2[1],MuLegEtaM2[2],MuLegEtaM2[3],MuLegEtaM2[4])); }
//             if(-0.8<etaMu <0){ weight_MuLeg18DATA.push_back(efficiency(ptMu,MuLegEtaM1[0],MuLegEtaM1[1],MuLegEtaM1[2],MuLegEtaM1[3],MuLegEtaM1[4]));     }
//             if(0<etaMu <0.8){weight_MuLeg18DATA.push_back(efficiency(ptMu,MuLegEtaP1[0],MuLegEtaP1[1],MuLegEtaP1[2],MuLegEtaP1[3],MuLegEtaP1[4]));     }
//             if(0.8<etaMu <1.2){weight_MuLeg18DATA.push_back(efficiency(ptMu,MuLegEtaP2[0],MuLegEtaP2[1],MuLegEtaP2[2],MuLegEtaP2[3],MuLegEtaP2[4]));     }
//             if(etaMu > 1.2){weight_MuLeg18DATA.push_back(efficiency(ptMu,MuLegEtaP3[0],MuLegEtaP3[1],MuLegEtaP3[2],MuLegEtaP3[3],MuLegEtaP3[4]));     }

//             if(etaMu <-1.2){weight_MuLeg18MC.push_back(efficiency(ptMu, MuLegEtaM3_mc[0], MuLegEtaM3_mc[1], MuLegEtaM3_mc[2], MuLegEtaM3_mc[3], MuLegEtaM3_mc[4])); }
//             if(-0.8<etaMu <-1.2){weight_MuLeg18MC.push_back(efficiency(ptMu,MuLegEtaM2_mc[0],MuLegEtaM2_mc[1],MuLegEtaM2_mc[2],MuLegEtaM2_mc[3],MuLegEtaM2_mc[4])); }
//             if(-0.8<etaMu <0){ weight_MuLeg18MC.push_back(efficiency(ptMu,MuLegEtaM1_mc[0],MuLegEtaM1_mc[1],MuLegEtaM1_mc[2],MuLegEtaM1_mc[3],MuLegEtaM1_mc[4]));     }
//             if(0<etaMu <0.8){weight_MuLeg18MC.push_back(efficiency(ptMu,MuLegEtaP1_mc[0],MuLegEtaP1_mc[1],MuLegEtaP1_mc[2],MuLegEtaP1_mc[3],MuLegEtaP1_mc[4]));     }
//             if(0.8<etaMu <1.2){weight_MuLeg18MC.push_back(efficiency(ptMu,MuLegEtaP2_mc[0],MuLegEtaP2_mc[1],MuLegEtaP2_mc[2],MuLegEtaP2_mc[3],MuLegEtaP2_mc[4]));     }
//             if(etaMu > 1.2){weight_MuLeg18MC.push_back(efficiency(ptMu,MuLegEtaP3_mc[0],MuLegEtaP3_mc[1],MuLegEtaP3_mc[2],MuLegEtaP3_mc[3],MuLegEtaP3_mc[4]));     }


//             weight_MuLeg18.push_back((weight_MuLeg18DATA.at(0)/weight_MuLeg18MC.at(0)));


//             ////////muon id and iso corrections//////

//             if(TMath::Abs(etaMu)<0.8 && ptMu<30) {weight_MuID= MuID_eta1[0];  weight_MuIso= MuIso_eta1[0];}
//             if((0.8<TMath::Abs(etaMu)<1.2) && ptMu<30) {weight_MuID= MuID_eta2[0]; weight_MuIso= MuIso_eta2[0];}
//             if((1.2<TMath::Abs(etaMu)<2.1) && ptMu<30) {weight_MuID= MuID_eta3[0]; weight_MuIso= MuIso_eta3[0];}

//             if(TMath::Abs(etaMu)<0.8 && ptMu>30) {weight_MuID= MuID_eta1[1];  weight_MuIso= MuIso_eta1[1];}
//             if((0.8<TMath::Abs(etaMu)<1.2) && ptMu>30) {weight_MuID= MuID_eta2[1]; weight_MuIso= MuIso_eta2[1];}
//             if((1.2<TMath::Abs(etaMu)<2.1) && ptMu>30) {weight_MuID= MuID_eta3[1]; weight_MuIso= MuIso_eta3[1];}


//    ////////////////////////////////////////////////// END Compute Trigger correction factors for Muon leg//////////////////////////////////////////////////

//     for(edm::View<pat::Tau>::const_iterator patTau=taus->begin(); patTau!=taus->end(); ++patTau){


//       double etaTau = patTau->eta();
//       double phiTau = patTau->phi();
//       double ptTau = patTau->pt();
//       std::cout<<" tau pt="<<patTau->pt()<<" eta="<<patTau->eta()<<" phi="<<patTau->phi()<<std::endl;

//       //////////////////////////////////////////////////Compute Trigger correction factors for tau leg //////////////////////////////////////////////////
//       if(TMath::Abs(etaTau)<1.5){
//         weight_TauLeg18DATA.push_back(efficiency(ptTau,TauLegB18[0],TauLegB18[1],TauLegB18[2],TauLegB18[3],TauLegB18[4]));
//         weight_TauLeg18MC.push_back(efficiency(ptTau,TauLegB18_MC[0],TauLegB18_MC[1],TauLegB18_MC[2],TauLegB18_MC[3],TauLegB18_MC[4]));
//         weight_TauLeg18.push_back( weight_TauLeg18DATA.at(0)/weight_TauLeg18MC.at(0));
//       }

//       else{
//         weight_TauLeg18DATA.push_back(efficiency(ptTau,TauLegE18[0],TauLegE18[1],TauLegE18[2],TauLegE18[3],TauLegE18[4]));
//         weight_TauLeg18MC.push_back(efficiency(ptTau,TauLegE18_MC[0],TauLegE18_MC[1],TauLegE18_MC[2],TauLegE18_MC[3],TauLegE18_MC[4]));
//         weight_TauLeg18.push_back(( weight_TauLeg18DATA.at(0)/weight_TauLeg18MC.at(0) ));
//       }



//       double weight_isomu18=(weight_TauLeg18.at(0))*(weight_MuLeg18.at(0));
//       //    double weight_isomu17=(weight_TauLeg17.at(0)*weight_MuLeg17.at(0));
//       double p1=890.183/12944.6;
//       double p2=12054/12944.6;



//       if(runNumber==1){
//         weight_HLT=weight_isomu18;
//         weightTotal= weight_HLT*weight_MuIso*weight_MuID;
//         //      std::cout<<"--------------------------nel loop----------------------- weight hlt="<<weight_HLT<<" runNumb"<<runNumber<<" -----------------------------------------------"<<std::endl;
//       }
//       else{weight_HLT=1;weight_MuIso=1;weight_MuID=1;weightTotal=1;}

//       //////////////////////////////////////////////////Compute Trigger correction factors for tau leg//////////////////////////////////////////////////

//       histContainer_["hTauPt"]->Fill(patTau->pt(), weight*weightTotal);
//       histContainer_["hTauEta"]->Fill(patTau->eta(), weight*weightTotal);
//       histContainer_["hTauPhi"]->Fill(patTau->phi(), weight*weightTotal);
//       histContainer_["hTauCharge"]->Fill(patTau->charge(), weight*weightTotal);
//       histContainer_["hTauVisMass"]->Fill(patTau->mass(), weight*weightTotal);

//    }
//}

} // Htautau_Summer13
} // cuts
