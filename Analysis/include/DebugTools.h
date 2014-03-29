/*!
 * \file DebugTools.h
 * \brief Common tools and definitions suitable for debug purposes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-03-29 created
 */

#pragma once

#include<vector>

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
