/*!
 * \file Scan.C
 * \brief ROOT macro to scan analysis tree file for the specified event.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \date 2014-11-25 created
 *
 * Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>
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

void Scan(UInt_t event_id, const char* file_name)
{
    TFile f(file_name, "READ");
    std::ostringstream ss;
    ss << "EventId==" << event_id;
    const TSQLResult* result = events->Query("EventId", ss.str().c_str());
    if(result->GetRowCount() > 0)
        std::cerr << "found in " << file_name << std::endl;
}

