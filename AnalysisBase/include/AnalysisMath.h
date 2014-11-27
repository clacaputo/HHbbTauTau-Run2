/*!
 * \file AnalysisMath.h
 * \brief Common math functions and definitions suitable for analysis purposes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-08 created
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

#include <iostream>
#include <cmath>

#include <TLorentzVector.h>
#include <TH1.h>

#include "AnalysisTypes.h"

namespace analysis {

//see AN-13-178
inline double Calculate_MT(const TLorentzVector& lepton_momentum, double met_pt, double met_phi)
{
    const double delta_phi = TVector2::Phi_mpi_pi(lepton_momentum.Phi() - met_phi);
    return std::sqrt( 2.0 * lepton_momentum.Pt() * met_pt * ( 1.0 - std::cos(delta_phi) ) );
}

// from DataFormats/TrackReco/interface/TrackBase.h
inline double Calculate_dxy(const TVector3& legV, const TVector3& PV, const TLorentzVector& leg)
{
    return ( - (legV.x()-PV.x()) * leg.Py() + (legV.y()-PV.y()) * leg.Px() ) / leg.Pt();
}

// from DataFormats/TrackReco/interface/TrackBase.h
inline double Calculate_dz(const TVector3& trkV, const TVector3& PV, const TVector3& trkP)
{
  return (trkV.z() - PV.z()) - ( (trkV.x() - PV.x()) * trkP.x() + (trkV.y() - PV.y()) * trkP.y() ) / trkP.Pt()
                               * trkP.z() / trkP.Pt();
}

inline TLorentzVector MakeLorentzVectorPtEtaPhiM(Double_t pt, Double_t eta, Double_t phi, Double_t m)
{
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, m);
    return v;
}

inline TLorentzVector MakeLorentzVectorPtEtaPhiE(Double_t pt, Double_t eta, Double_t phi, Double_t e)
{
    TLorentzVector v;
    v.SetPtEtaPhiE(pt, eta, phi, e);
    return v;
}

inline PhysicalValue Integral(const TH1D& histogram, bool include_overflows = true)
{
    typedef std::pair<Int_t, Int_t> limit_pair;
    const limit_pair limits = include_overflows ? limit_pair(0, histogram.GetNbinsX() + 1)
                                                : limit_pair(1, histogram.GetNbinsX());

    double error = 0;
    const double integral = histogram.IntegralAndError(limits.first, limits.second, error);
    return PhysicalValue(integral, error);
}

inline void RenormalizeHistogram(TH1D& histogram, const PhysicalValue norm, bool include_overflows = true)
{
    histogram.Scale(norm.value/Integral(histogram,include_overflows).value);
}

} // namespace analysis
