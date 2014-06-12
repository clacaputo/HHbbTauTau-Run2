#pragma once


#include <string>
#include <vector>

namespace NSVfitStandalone {

std::string replace_string(const std::string&, const std::string&, const std::string&, unsigned, unsigned, int&);

std::string format_vstring(const std::vector<std::string>&);
std::string format_vdouble(const std::vector<double>&);
std::string format_vfloat(const std::vector<float>&);
std::string format_vint(const std::vector<int>&);
std::string format_vunsigned(const std::vector<unsigned>&);
std::string format_vbool(const std::vector<bool>&);

} // namespace NSVfitStandalone
