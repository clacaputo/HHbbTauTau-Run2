/*!
 * \file ProgressReporter.h
 * \brief Definition of the class to report current execution progress.
 * \author Konstantin Androsov (University of Siena, INFN Pisa)
 * \date 2014-05-05 created
 * \date 2015-02-14 separated from Tools.h
 *
 * Copyright 2014, 2015 Konstantin Androsov <konstantin.androsov@gmail.com>
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

#include <sstream>
#include <iomanip>
#include <chrono>

namespace analysis {
namespace tools {
class ProgressReporter {
private:
    typedef std::chrono::system_clock clock;

public:
    static std::string TimeStamp(const clock::time_point& time_point)
    {
        std::ostringstream ss;
        const std::time_t time_t_point = clock::to_time_t(time_point);
		//ss << "[" << std::put_time(std::localtime(&time_t_point), "%F %T") << "] ";
		char mbstr[100];
		if (std::strftime(mbstr,sizeof(mbstr),"%F %T",std::localtime(&time_t_point)))
			ss << "[" << mbstr << "] ";
        return ss.str();
    }

    ProgressReporter(unsigned _report_interval, std::ostream& _output)
        : start(clock::now()), block_start(start), report_interval(_report_interval), output(&_output)
    {
        *output << TimeStamp(start) << "Starting analyzer..." << std::endl;
    }

    void Report(size_t event_id, bool final_report = false)
    {
        using namespace std::chrono;

        const auto now = clock::now();
        const unsigned since_last_report = duration_cast<seconds>(now - block_start).count();
        if(!final_report && since_last_report < report_interval) return;

        const unsigned since_start = duration_cast<seconds>(now - start).count();
        const double speed = double(event_id) / since_start;
        *output << TimeStamp(now);
        if(final_report)
            *output << "Total: ";
        *output << "time = " << since_start << " seconds, events processed = " << event_id
                  << ", average speed = " << std::setprecision(1) << std::fixed << speed << " events/s" << std::endl;
        const unsigned since_start_residual = since_start % report_interval;
        block_start = now - seconds(since_start_residual);
    }

private:
    clock::time_point start, block_start;
    unsigned report_interval;
    std::ostream* output;
};
} // namespace tools
} // namespace analysis
