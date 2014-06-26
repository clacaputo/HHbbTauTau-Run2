/*!
 * \file Event.h
 * \brief Definiton of ntuple::EventTree and ntuple::Event classes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-03-24 created
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

#define EVENT_DATA() \
    SIMPLE_VAR(UInt_t, run) \
    SIMPLE_VAR(UInt_t, EventId) \
    SIMPLE_VAR(UInt_t, lumis) \
    SIMPLE_VAR(Int_t, bunch) \
    SIMPLE_VAR(Int_t, orbit) \
    SIMPLE_VAR(UInt_t, unixTime) \
    SIMPLE_VAR(UInt_t, microsecondOffset) \
    SIMPLE_VAR(Bool_t, isdata) \
    SIMPLE_VAR(Bool_t, isPhysDeclared)  \
    SIMPLE_VAR(Bool_t, isBPTX0) \
    SIMPLE_VAR(Bool_t, isBSCMinBias) \
    SIMPLE_VAR(Bool_t, isBSCBeamHalo) \
    SIMPLE_VAR(Bool_t, isPrimaryVertex) \
    VECTOR_VAR(Int_t, nPU) \
    VECTOR_VAR(Int_t, bunchCrossing) \
    VECTOR_VAR(Float_t, trueNInt) \
    /**/

#define SIMPLE_VAR(type, name) DECLARE_SIMPLE_BRANCH_VARIABLE(type, name)
#define VECTOR_VAR(type, name) DECLARE_VECTOR_BRANCH_VARIABLE(type, name)
DATA_CLASS(ntuple, Event, EVENT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) SIMPLE_DATA_TREE_BRANCH(type, name)
#define VECTOR_VAR(type, name) VECTOR_DATA_TREE_BRANCH(type, name)
TREE_CLASS(ntuple, EventTree, EVENT_DATA, Event, "events", false)
#undef SIMPLE_VAR
#undef VECTOR_VAR

#define SIMPLE_VAR(type, name) ADD_SIMPLE_DATA_TREE_BRANCH(name)
#define VECTOR_VAR(type, name) ADD_VECTOR_DATA_TREE_BRANCH(name)
TREE_CLASS_INITIALIZE(ntuple, EventTree, EVENT_DATA)
#undef SIMPLE_VAR
#undef VECTOR_VAR
#undef EVENT_DATA
