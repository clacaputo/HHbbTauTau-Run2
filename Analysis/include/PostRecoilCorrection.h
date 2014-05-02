/*!
 * \file PostRecoilCorrection.h
 * \brief Definition of wrapper for RecoilCorrector code to apply post-recoil correction.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-05-02 created
 */

#pragma once

#include <TLorentzVector.h>

#include "HHbbTauTau/RecoilCorrector_v7/RecoilCorrector.hh"
#include "HHbbTauTau/TreeProduction/interface/MET.h"

namespace analysis {
ntuple::MET ApplyPostRecoilCorrection(const ntuple::MET& originalMET, const TLorentzVector& resonantMomentum,
                                      const TLorentzVector& resonantMomentumMC)
{
    double iU1, iU2;
    double met = originalMET.pt_uncorrected;
    double metphi = originalMET.phi_uncorrected;
    RecoilCorrector corrector("RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_20pv_njet.root");
    corrector.addMCFile("RecoilCorrector_v7/recoilfits/recoilfit_higgs53X_20pv_njet.root");
    corrector.addDataFile("RecoilCorrector_v7/recoilfits/recoilfit_datamm53X_20pv_njet.root");
    corrector.CorrectType1(met, metphi, resonantMomentumMC.Pt(), resonantMomentumMC.Phi(),
                           resonantMomentum.Pt(), resonantMomentum.Phi(), iU1, iU2, 0);
//    std::cerr << "I'm here" << std::endl;
    ntuple::MET correctedMET(originalMET);
    correctedMET.pt = met;
    correctedMET.phi = metphi;
    return correctedMET;
}
} // analysis
