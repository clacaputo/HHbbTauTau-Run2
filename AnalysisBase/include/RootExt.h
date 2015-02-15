/*!
 * \file RootExt.h
 * \brief Common CERN ROOT extensions.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \date 2015-02-15 combined from different files
 *
 * Copyright 2013-2015 Konstantin Androsov <konstantin.androsov@gmail.com>
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

#include <memory>

#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TFile.h>
#include <Compression.h>

#include "exception.h"

namespace root_ext {

std::shared_ptr<TFile> CreateRootFile(const std::string& file_name)
{
    std::shared_ptr<TFile> file(TFile::Open(file_name.c_str(), "RECREATE", "", ROOT::kZLIB * 100 + 9));
    if(file->IsZombie())
        throw analysis::exception("File '") << file_name << "' not created.";
    return file;
}

std::shared_ptr<TFile> OpenRootFile(const std::string& file_name)
{
    std::shared_ptr<TFile> file(TFile::Open(file_name.c_str(), "READ"));
    if(file->IsZombie())
        throw analysis::exception("File '") << file_name << "' not opened.";
    return file;
}

} // namespace root_ext


std::ostream& operator<<(std::ostream& s, const TVector3& v) {
    s << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")";
    return s;
}

std::ostream& operator<<(std::ostream& s, const TLorentzVector& v) {
    s << "(pt=" << v.Pt() << ", eta=" << v.Eta() << ", phi=" << v.Phi() << ", E=" << v.E() << ", m=" << v.M() << ")";
    return s;
}

// Based on TMatrixD::Print code.
std::ostream& operator<<(std::ostream& s, const TMatrixD& matrix)
{
    if (!matrix.IsValid()) {
        s << "Matrix is invalid";
        return s;
    }

    //build format
    const char *format = "%11.4g ";
    char topbar[100];
    snprintf(topbar,100,format,123.456789);
    Int_t nch = strlen(topbar)+1;
    if (nch > 18) nch = 18;
    char ftopbar[20];
    for (Int_t i = 0; i < nch; i++) ftopbar[i] = ' ';
    Int_t nk = 1 + Int_t(std::log10(matrix.GetNcols()));
    snprintf(ftopbar+nch/2,20-nch/2,"%s%dd","%",nk);
    Int_t nch2 = strlen(ftopbar);
    for (Int_t i = nch2; i < nch; i++) ftopbar[i] = ' ';
    ftopbar[nch] = '|';
    ftopbar[nch+1] = 0;

    s << matrix.GetNrows() << "x" << matrix.GetNcols() << " matrix";

    Int_t cols_per_sheet = 5;
    if (nch <= 8) cols_per_sheet =10;
    const Int_t ncols  = matrix.GetNcols();
    const Int_t nrows  = matrix.GetNrows();
    const Int_t collwb = matrix.GetColLwb();
    const Int_t rowlwb = matrix.GetRowLwb();
    nk = 5+nch*std::min(cols_per_sheet, matrix.GetNcols());
    for (Int_t i = 0; i < nk; i++)
        topbar[i] = '-';
    topbar[nk] = 0;
    for (Int_t sheet_counter = 1; sheet_counter <= ncols; sheet_counter += cols_per_sheet) {
        s << "\n     |";
        for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++) {
            char ftopbar_out[100];
            snprintf(ftopbar_out, 100, ftopbar, j+collwb-1);
            s << ftopbar_out;
        }
        s << "\n" << topbar << "\n";
        if (matrix.GetNoElements() <= 0) continue;
        for (Int_t i = 1; i <= nrows; i++) {
            char row_out[100];
            snprintf(row_out, 100, "%4d |",i+rowlwb-1);
            s << row_out;
            for (Int_t j = sheet_counter; j < sheet_counter+cols_per_sheet && j <= ncols; j++) {
                snprintf(row_out, 100, format, matrix(i+rowlwb-1,j+collwb-1));
                s << row_out;
            }
            s << "\n";
        }
    }
    return s;
}
