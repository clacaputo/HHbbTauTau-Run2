/*!
 * \file DebugTools.h
 * \brief Common tools and definitions suitable for debug purposes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-03-29 created
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

#include <vector>

namespace debug {

#define PRINT_SIZEOF(name) \
    std::cout << "Sizeof " #name " = " << sizeof(name) << std::endl

inline void PrintCommonTypeSizes()
{
    PRINT_SIZEOF(Float_t);
    PRINT_SIZEOF(Double_t);
    PRINT_SIZEOF(Int_t);
    PRINT_SIZEOF(UInt_t);
    PRINT_SIZEOF(Bool_t);
    PRINT_SIZEOF(Long64_t);
    PRINT_SIZEOF(int);
    PRINT_SIZEOF(unsigned);
    PRINT_SIZEOF(float);
    PRINT_SIZEOF(double);
    PRINT_SIZEOF(bool);
    PRINT_SIZEOF(size_t);

    typedef std::vector<double>::size_type std_collection_size_type;
    PRINT_SIZEOF(std_collection_size_type);
}
#undef PRINT_SIZEOF

} // debug
