/*!
 * \file Vertex.h
 * \brief Definiton of ntuple::VertexTree and ntuple::Vertex classes.
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

#define VERTEX_DATA() \
    SIMPLE_VAR(Float_t, x) \
    SIMPLE_VAR(Float_t, y) \
    SIMPLE_VAR(Float_t, z) \
    SIMPLE_VAR(Float_t, xErr) \
    SIMPLE_VAR(Float_t, yErr) \
    SIMPLE_VAR(Float_t, zErr) \
    SIMPLE_VAR(Float_t, rho) \
    SIMPLE_VAR(Float_t, chi2) \
    SIMPLE_VAR(Float_t, ndf) \
    SIMPLE_VAR(UInt_t, ntracks) \
    SIMPLE_VAR(UInt_t, ntracksw05) \
    SIMPLE_VAR(Bool_t, isfake) \
    SIMPLE_VAR(Bool_t, isvalid) \
    SIMPLE_VAR(Float_t, sumPt) \
    SIMPLE_VAR(Float_t, sumPtSquared) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Vertex, VERTEX_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS_WITH_EVENT_ID(ntuple, VertexTree, VERTEX_DATA, Vertex, "vertices", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_WITH_EVENT_ID_INITIALIZE(ntuple, VertexTree, VERTEX_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef VERTEX_DATA
