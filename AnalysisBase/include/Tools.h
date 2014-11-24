/*!
 * \file Tools.h
 * \brief Common tools and definitions suitable for general purposes.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-05-05 created
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
#include <set>
#include <chrono>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace tools {

class Timer {
public:
    typedef std::chrono::high_resolution_clock clock;
    Timer(unsigned _report_interval)
        : start(clock::now()), block_start(start), report_interval(_report_interval)
    {
        std::cout << "Starting analyzer..." << std::endl;
    }

    void Report(size_t event_id, bool final_report = false)
    {
        using namespace std::chrono;
        const auto now = clock::now();
        const auto since_last_report = duration_cast<seconds>(now - block_start).count();
        if(!final_report && since_last_report < report_interval) return;

        const auto since_start = duration_cast<seconds>(now - start).count();
        const double speed = ((double) event_id) / since_start;
        if(final_report)
            std::cout << "Total: ";
        std::cout << "time = " << since_start << " seconds, events processed = " << event_id
                  << ", average speed = " << std::setprecision(1) << std::fixed << speed << " events/s" << std::endl;
        block_start = now;
    }

private:
    clock::time_point start, block_start;
    unsigned report_interval;
};

template<typename Type>
std::vector<Type> join_vectors(const std::vector< const std::vector<Type>* >& inputVectors)
{
    size_t totalSize = 0;
    for(auto inputVector : inputVectors) {
        if(!inputVector)
            throw std::runtime_error("input vector is nullptr");
        totalSize += inputVector->size();
    }

    std::vector<Type> result;
    result.reserve(totalSize);
    for(auto inputVector : inputVectors)
        result.insert(result.end(), inputVector->begin(), inputVector->end());

    return result;
}

template<typename Container, typename T>
size_t find_index(const Container& container, const T& value)
{
    const auto iter = std::find(container.begin(), container.end(), value);
    return std::distance(container.begin(), iter);
}

template<typename Map>
std::set< typename Map::key_type > collect_map_keys(const Map& map)
{
    std::set< typename Map::key_type > result;
    std::transform(map.begin(), map.end(), std::inserter(result, result.end()),
                   [](const typename Map::value_type& pair) { return pair.first; } );
    return result;
}

template<typename Map>
std::set< typename Map::mapped_type > collect_map_values(const Map& map)
{
    std::set< typename Map::mapped_type > result;
    std::transform(map.begin(), map.end(), std::inserter(result, result.end()),
                   [](const typename Map::value_type& pair) { return pair.second; } );
    return result;
}

} // tools
