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
                                      const TLorentzVector& resonantMomentumMC, const size_t njets)
{
    //from Riccardo
    static const std::string fileCorrectTo = "RecoilCorrector_v7/recoilfits/recoilfit_ztt53X_20pv_njet.root";
    static const std::string fileZmmData = "RecoilCorrector_v7/recoilfits/recoilfit_datamm53XRR_2012_njet.root";
    static const std::string fileZmmMC = "RecoilCorrector_v7/recoilfits/recoilfit_zmm53XRR_2012_njet.root";
    double iU1, iU2;
    double met = originalMET.pt;
    double metphi = originalMET.phi;
    RecoilCorrector corrector(fileCorrectTo);
    corrector.addDataFile(fileZmmData);
    corrector.addMCFile(fileZmmMC);

    corrector.CorrectType2(met, metphi, resonantMomentumMC.Pt(), resonantMomentumMC.Phi(),
                           resonantMomentum.Pt(), resonantMomentum.Phi(), iU1, iU2, 0, 0, njets);

    ntuple::MET correctedMET(originalMET);
    correctedMET.pt = met;
    correctedMET.phi = metphi;
    return correctedMET;
}
} // analysis
