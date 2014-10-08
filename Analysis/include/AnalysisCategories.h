/*!
 * \file AnalysisCategories.h
 * \brief Definition of data and event categories used in HH->bbTauTau analysis.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-09-16 created
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
#include <list>
#include <map>
#include <set>
#include <cmath>

#include <TFile.h>
#include <Rtypes.h>

#include "AnalysisBase/include/exception.h"

namespace analysis {

typedef std::map<std::string, double> DataSourceScaleFactorMap;

enum class DataCategoryType { Signal, Background, Data, DYJets, ZL, ZJ, ZTT, ZTT_MC, Embedded, Limits, Sum, QCD, WJets,
                            EWK};
static const std::map<DataCategoryType, std::string> dataCategoryTypeNameMap = {
    { DataCategoryType::Signal, "SIGNAL" }, { DataCategoryType::Background, "BACKGROUND" },
    { DataCategoryType::Data, "DATA" }, { DataCategoryType::DYJets, "DY_JETS" }, { DataCategoryType::ZL, "ZL" },
    { DataCategoryType::ZJ, "ZJ" }, { DataCategoryType::ZTT, "ZTT" }, { DataCategoryType::ZTT_MC, "ZTT_MC" },
    { DataCategoryType::Embedded, "EMBEDDED" }, { DataCategoryType::Limits, "LIMITS" },
    { DataCategoryType::Sum, "SUM" }, { DataCategoryType::QCD, "QCD" }, { DataCategoryType::WJets, "W_JETS" },
    { DataCategoryType::EWK, "EWK" }
};

std::ostream& operator<< (std::ostream& s, const DataCategoryType& dataCategoryType) {
    s << dataCategoryTypeNameMap.at(dataCategoryType);
    return s;
}
std::istream& operator>> (std::istream& s, DataCategoryType& dataCategoryType) {
    std::string name;
    s >> name;
    for(const auto& map_entry : dataCategoryTypeNameMap) {
        if(map_entry.second == name) {
            dataCategoryType = map_entry.first;
            return s;
        }
    }
    throw exception("Unknown data category type '") << name << "'.";
}

struct DataCategory {
    std::string name;
    std::string title;
    std::string datacard;
    EColor color;
    double limits_sf;
    bool draw;
    unsigned draw_sf;

    std::set<DataCategoryType> types;
    std::set<std::string> channels;
    DataSourceScaleFactorMap sources_sf;

    DataCategory()
        : color(kBlack), limits_sf(std::nan("")), draw(false), draw_sf(1) {}

    bool IsSignal() const { return types.count(DataCategoryType::Signal); }
    bool IsBackground() const { return types.count(DataCategoryType::Background); }
    bool IsData() const { return types.count(DataCategoryType::Data); }
};

typedef std::map<std::string, DataCategory> DataCategoryMap;
typedef std::set<const DataCategory*> DataCategoryPtrSet;
typedef std::map<DataCategoryType, DataCategoryPtrSet> DataCategoryTypeMap;

struct DataSource {
    std::string file_name;
    DataCategoryPtrSet data_categories;

    double scale_factor(const DataCategory* category) const
    {
        return category->sources_sf.at(file_name);
    }
};

typedef std::set<const DataSource*> DataSourcePtrSet;
typedef std::map<std::string, DataSource> DataSourceMap;

class DataCategoryCollection {
public:
    DataCategoryCollection(const std::string& sources_cfg_name, const std::string& signal_list,
                           const std::string& channel_name)
    {
        DataCategory category;
        std::ifstream cfg(sources_cfg_name);
        size_t line_number = 0;
        while(ReadNextCategory(cfg, line_number, category)) {
            if(categories.count(category.name))
                throw exception("Category with name '") << category.name << "' is already defined.";
            if(category.channels.size() && !category.channels.count(channel_name)) continue;
            categories[category.name] = category;
            all_categories.insert(&categories[category.name]);
            for(DataCategoryType type : category.types)
                categories_by_type[type].insert(&categories[category.name]);
        }
        for(const auto& category_entry : categories) {
            for(const auto& sources_sf_entry : category_entry.second.sources_sf) {
                sources[sources_sf_entry.first].file_name = sources_sf_entry.first;
                sources[sources_sf_entry.first].data_categories.insert(&category_entry.second);
                all_sources.insert(&sources[sources_sf_entry.first]);
            }
        }
        const auto& signal_names = ParseSignalList(signal_list);
        for(const auto& signal_name : signal_names)
            categories[signal_name].draw = true;
    }

    const DataSourcePtrSet& GetSources() const { return all_sources; }
    const DataCategoryPtrSet& GetAllCategories() const { return all_categories; }
    const DataCategoryPtrSet& GetCategories(DataCategoryType dataCategoryType) const
    {
        return categories_by_type.at(dataCategoryType);
    }

    const DataCategory& GetUniqueCategory(DataCategoryType dataCategoryType) const
    {
        if(!categories_by_type.at(dataCategoryType).size())
            throw exception("Unique category for data category type '") << dataCategoryType << "' not found.";
        if(categories_by_type.at(dataCategoryType).size() != 1)
            throw exception("More than one category for data category type '") << dataCategoryType << "'.";
        return *(*categories_by_type.at(dataCategoryType).begin());
    }


    const DataCategory& FindCategory(const std::string& name) const
    {
        if(!categories.count(name))
            throw exception("Data category '") << name << "' not found.";
        return categories.at(name);
    }

private:
    static bool ReadNextCategory(std::istream& cfg, size_t line_number, DataCategory& category)
    {
        category = DataCategory();
        bool category_started = false;
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            ++line_number;
            if (cfgLine.at(0) == '#' || (!cfgLine.size() && !category_started)) continue;
            if(!cfgLine.size())
                return true;
            if(!category_started && cfgLine.at(0) == '[') {
                const size_t pos = cfgLine.find(']');
                if(pos == std::string::npos)
                    throw exception("bad source config syntax in line ") << line_number;
                category.name = cfgLine.substr(1, pos - 1);
                category_started = true;
            } else if(category_started) {
                ReadParameterLine(cfgLine, line_number, category);
            } else
                throw exception("bad source config syntax in line ") << line_number;
        }
        return category_started;
    }

    static void ReadParameterLine(const std::string& cfgLine, size_t line_number, DataCategory& category)
    {
        static const char separator = ',';

        const size_t pos = cfgLine.find(separator);
        const std::string param_name = cfgLine.substr(0, pos - 1);
        const std::string param_value = cfgLine.substr(pos + 2);
        std::istringstream ss(param_value);
        if(param_name == "type") {
            DataCategoryType type;
            ss >> type;
            category.types.insert(type);
            category.draw |= type == DataCategoryType::Background || type == DataCategoryType::Data;
        } else if(param_name == "title") {
            category.title = param_value;
        } else if(param_name == "color") {
            ss >> category.color;
        } else if(param_name == "file") {
            std::string file_name;
            double scale_factor;
            ss >> file_name;
            ss >> scale_factor;
            category.sources_sf[file_name] = scale_factor;
        } else if(param_name == "limits_sf") {
            ss >> category.limits_sf;
        } else if(param_name == "draw_sf") {
            ss >> category.draw_sf;
        } else if(param_name == "channel") {
            std::string channel_name;
            ss >> channel_name;
            category.channels.insert(channel_name);
        } else if(param_name == "datacard") {
            ss >> category.datacard;
        } else
            throw exception("Unsupported parameter '") << param_name << "' in configuration line " << line_number;
    }

    static std::set<std::string> ParseSignalList(const std::string& signal_list)
    {
        static const char separator = ',';

        std::set<std::string> result;
        size_t prev_pos = 0;
        for(size_t pos = signal_list.find(separator); pos != std::string::npos;
                                                      pos = signal_list.find(separator, prev_pos)) {
            const std::string signal_name = signal_list.substr(prev_pos, pos - 1);
            result.insert(signal_name);
            prev_pos = pos + 1;
        }
        return result;
    }

private:
    DataSourceMap sources;
    DataSourcePtrSet all_sources;
    DataCategoryMap categories;
    DataCategoryPtrSet all_categories;
    DataCategoryTypeMap categories_by_type;
};


std::ostream& operator<<(std::ostream& s, const DataSource& source){
    s << "File: " << source.file_name << " for data categories:\n";
    for(const DataCategory* category : source.data_categories)
        s << category->name << ", SF: " << source.scale_factor(category) << "\n";
    return s;
}

std::ostream& operator<<(std::ostream& s, const DataCategory& category){
    s << "Name: " << category.name << ", Title: '" << category.title << "', Color: " << category.color << std::endl;
    for(const auto& source : category.sources_sf)
        s << "File: " << source.first << ", SF: " << source.second << "\n";
    return s;
}


enum class EventType_QCD { Unknown, OS_Isolated, OS_NotIsolated, SS_Isolated, SS_NotIsolated };
enum class EventType_Wjets { Unknown, Signal, HighMt };
enum class EventCategory { Inclusive, OneJet_ZeroBtag, OneJet_OneBtag, TwoJets_ZeroBtag, TwoJets_OneBtag, TwoJets_TwoBtag };

static const std::map<EventCategory, std::string> eventCategoryMapName =
          { { EventCategory::Inclusive, "Inclusive" }, { EventCategory::OneJet_ZeroBtag, "1jet0btag" },
            { EventCategory::OneJet_OneBtag, "1jet1btag" }, { EventCategory::TwoJets_ZeroBtag, "2jets0btag" },
          { EventCategory::TwoJets_OneBtag, "2jets1btag"}, { EventCategory::TwoJets_TwoBtag, "2jets2btag" } };
typedef std::vector<EventCategory> EventCategoryVector;

std::ostream& operator<<(std::ostream& s, const EventCategory& eventCategory) {
    s << eventCategoryMapName.at(eventCategory);
    return s;
}

std::wostream& operator<<(std::wostream& s, const EventCategory& eventCategory) {
    const std::string str = eventCategoryMapName.at(eventCategory);
    s << std::wstring(str.begin(), str.end());
    return s;
}

} // namespace analysis
