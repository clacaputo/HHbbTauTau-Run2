/*!
 * \file KinFitTree.h
 * \brief Definiton of ntuple::KinFitTree and ntuple::KinFit classes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-24 created
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
#include <limits>
#include "TreeProduction/interface/SmartTree.h"

#define KINFIT_DATA() \
    /* Event Variables */ \
    SIMPLE_VAR(Int_t, evt) /* Event */ \
    SIMPLE_VAR(Int_t, channel) /* Analysis channel as defined in analysis::Channel */ \
    SIMPLE_VAR(Int_t, eventRegion) /*Sideband region for estimation of bkg*/ \
    SIMPLE_VAR(Int_t, eventCategory) /*event category*/ \
    \
    \
    /* First signal lepton :  muon for MuTau, electron for ETau, leading (in pT) tau for TauTau */ \
    SIMPLE_VAR(Float_t, pt_1) /* pT */ \
    SIMPLE_VAR(Float_t, phi_1) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_1) /* Eta */ \
    SIMPLE_VAR(Float_t, energy_1) /* Energy */ \
    \
    \
    /* Second lepton : hadronic tau for MuTau & ETau, trailing (in pT) tau for TauTau */ \
    SIMPLE_VAR(Float_t, pt_2) /* pT */ \
    SIMPLE_VAR(Float_t, phi_2) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_2) /* Eta */ \
    SIMPLE_VAR(Float_t, energy_2) /* Energy */ \
    \
    \
    /* First bjet : leading in CSV */ \
    SIMPLE_VAR(Float_t, pt_b1) /* pT */ \
    SIMPLE_VAR(Float_t, phi_b1) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_b1) /* Eta */ \
    SIMPLE_VAR(Float_t, energy_b1) /* Energy */ \
    SIMPLE_VAR(Float_t, csv_b1) /* csv */ \
    \
    \
    /* Second bjet : subleading in CSV */ \
    SIMPLE_VAR(Float_t, pt_b2) /* pT */ \
    SIMPLE_VAR(Float_t, phi_b2) /* Phi */ \
    SIMPLE_VAR(Float_t, eta_b2) /* Eta */ \
    SIMPLE_VAR(Float_t, energy_b2) /* Energy */ \
    SIMPLE_VAR(Float_t, csv_b2) /* csv */ \
    \
    \
    /* Met related variables */ \
    SIMPLE_VAR(Float_t, mvamet) /* mvamet */ \
    SIMPLE_VAR(Float_t, mvametphi) /* mvamet Phi */ \
    /* MVAMet covariance matrices */ \
    SIMPLE_VAR(Float_t, mvacov00) /* mva met covariance matrix 00 */ \
    SIMPLE_VAR(Float_t, mvacov01) /* mva met covariance matrix 01 */ \
    SIMPLE_VAR(Float_t, mvacov10) /* mva met covariance matrix 10 */ \
    SIMPLE_VAR(Float_t, mvacov11) /* mva met covariance matrix 11 */ \
    \
    \
    /* H_tautau variables from SVfit*/ \
    SIMPLE_VAR(Float_t, m_Htt) /* Mass of H_tautau corrected by svFit using integration method MC */ \
    SIMPLE_VAR(Float_t, pt_Htt) /* Pt of H_tautau corrected by svFit using integration method MC */ \
    SIMPLE_VAR(Float_t, eta_Htt) /* eta of original H_tautau without MVAMET */ \
    SIMPLE_VAR(Float_t, phi_Htt) /* phi of original H_tautau without MVAMET */ \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, KinFit, KINFIT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, KinFitTree, KINFIT_DATA, KinFit, "KinFit", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, KinFitTree, KINFIT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef KINFIT_DATA

