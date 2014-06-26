/*!
 * \file MET.h
 * \brief Definiton of ntuple::METTree and ntuple::MET classes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-03-25 created
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

#include "SmartTree.h"
#include "TMatrixD.h"

#define MET_DATA() \
    SIMPLE_VAR(Float_t, pt) \
    SIMPLE_VAR(Float_t, phi) \
    SIMPLE_VAR(Float_t, sumEt) \
    SIMPLE_VAR(Float_t, pt_uncorrected) \
    SIMPLE_VAR(Float_t, phi_uncorrected) \
    SIMPLE_VAR(Float_t, sumEt_uncorrected) \
    VECTOR_VAR(Float_t, significanceMatrix) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, MET, MET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, METTree, MET_DATA, MET, "METs", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, METTree, MET_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef MET_DATA


namespace ntuple {
    std::vector<Float_t> SignificanceMatrixToVector(const TMatrixD& m)
    {
        if(m.GetNrows() != 2 || m.GetNcols() != 2)
            throw std::runtime_error("invalid matrix dimensions");
        std::vector<Float_t> v(4);
        v[0] = m[0][0];
        v[1] = m[0][1];
        v[2] = m[1][0];
        v[3] = m[1][1];
        return v;
    }

    TMatrixD VectorToSignificanceMatrix(const std::vector<Float_t>& v)
    {
        if(v.size() != 4)
            throw std::runtime_error("invalid input vector for significance matrix");
        TMatrixD m(2, 2);
        m[0][0] = v[0];
        m[0][1] = v[1];
        m[1][0] = v[2];
        m[1][1] = v[3];
        return m;
    }
}
