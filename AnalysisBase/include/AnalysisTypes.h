/*!
 * \file AnalysisMath.h
 * \brief Common simple types for analysis purposes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-09 created
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
#include <iomanip>
#include <map>
#include <cmath>

#include "exception.h"

namespace analysis {

enum class Channel { ETau = 0, MuTau = 1, TauTau = 2 };

namespace detail {
std::map<Channel, std::string> ChannelNameMap = {
    { Channel::ETau, "eTau" }, { Channel::MuTau, "muTau" }, { Channel::TauTau, "tauTau" }
};
} // namespace detail

std::ostream& operator<< (std::ostream& s, const Channel& c)
{
    s << detail::ChannelNameMap.at(c);
    return s;
}

std::istream& operator>> (std::istream& s, Channel& c)
{
    std::string name;
    s >> name;
    for(const auto& map_entry : detail::ChannelNameMap) {
        if(map_entry.second == name) {
            c = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown channel name '") << name << "'.";
}

template<typename T>
T sqr(const T& x) { return x * x; }

namespace detail {
template<typename char_type>
struct PhysicalValueErrorSeparator;

template<>
struct PhysicalValueErrorSeparator<char> {
    static std::string Get() { return " +/- "; }
};

template<>
struct PhysicalValueErrorSeparator<wchar_t> {
    static std::wstring Get() { return L" \u00B1 "; }
};

}

struct PhysicalValue {

    double value;
    double error;

    PhysicalValue() : value(0), error(0) {}
    PhysicalValue(double _value, double _error) : value(_value), error(_error)
    {
        if(error < 0)
            throw exception("Negative error = ") << error << ".";
    }

    PhysicalValue operator+(const PhysicalValue& other) const
    {
        const double new_error = std::sqrt(sqr(error) + sqr(other.error));
        return PhysicalValue(value + other.value, new_error);
    }

    PhysicalValue operator-(const PhysicalValue& other) const
    {
        const double new_error = std::sqrt(sqr(error) + sqr(other.error));
        return PhysicalValue(value - other.value, new_error);
    }

    PhysicalValue operator*(const PhysicalValue& other) const
    {
        const double new_error = std::sqrt(sqr(other.value * error) + sqr(value * other.error));
        return PhysicalValue(value * other.value, new_error);
    }

    PhysicalValue operator/(const PhysicalValue& other) const
    {
        const double new_error = std::sqrt(sqr(error) + sqr(value * other.error / sqr(other.value)))
                / std::abs(other.value);
        return PhysicalValue(value / other.value, new_error);
    }

    bool operator<(const PhysicalValue& other) const { return value < other.value; }

    bool IsCompatible(const PhysicalValue& other) const
    {
        const double delta = std::abs(value - other.value);
        const double max_error = std::max(error, other.error);
        return delta < max_error;
    }

    template<typename char_type>
    std::basic_string<char_type> ToString(bool print_error) const
    {
        static const int number_of_significant_digits_in_error = 2;
        const int precision = error ? std::floor(std::log10(error)) - number_of_significant_digits_in_error + 1
                                    : -15;
        const double ten_pow_p = std::pow(10.0, precision);
        const double error_rounded = std::ceil(error / ten_pow_p) * ten_pow_p;
        const double value_rounded = std::round(value / ten_pow_p) * ten_pow_p;
        const int decimals_to_print = std::max(0, -precision);
        std::basic_ostringstream<char_type> ss;
        ss << std::setprecision(decimals_to_print) << std::fixed << value_rounded;
        if(print_error)
            ss << detail::PhysicalValueErrorSeparator<char_type>::Get() << error_rounded;
        return ss.str();
    }
};

typedef std::pair<PhysicalValue, PhysicalValue> PhysicalValuePair;

std::ostream& operator<<(std::ostream& s, const PhysicalValue& v)
{
    s << v.ToString<char>(true);
    return s;
}

std::wostream& operator<<(std::wostream& s, const PhysicalValue& v)
{
    s << v.ToString<wchar_t>(true);
    return s;
}

} // namespace analysis
