/*!
 * \file exception.h
 * \brief Definition of the base exception class for the analysis namespace.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
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

#pragma once

#include <exception>
#include <string>
#include <sstream>

namespace analysis {
class exception : public std::exception {
public:
    explicit exception() noexcept {}
    explicit exception(const std::string& message) noexcept { s_msg << message; }
    exception(const exception& other) noexcept { s_msg << other.message(); }
    virtual ~exception() noexcept {}
    virtual const char* what() const noexcept { msg = s_msg.str(); return msg.c_str(); }
    const std::string& message() const noexcept { msg = s_msg.str(); return msg; }

    template<typename T>
    exception& operator << (const T& t) { s_msg << t; return *this; }

private:
    mutable std::string msg;
    std::ostringstream s_msg;
};

} // analysis
