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

#include "AnalysisBase/include/AnalysisTypes.h"
#include "AnalysisBase/include/Tools.h"

namespace analysis {

typedef std::map<std::string, double> DataSourceScaleFactorMap;

enum class DataCategoryType { Signal, Background, Data, DYJets, DYJets_incl, DYJets_excl, ZL, ZJ, ZTT, ZTT_MC, ZTT_L,
                              Embedded, Limits, Composit, QCD, WJets, WJets_MC, WJets_MC_incl, WJets_MC_excl};
static const std::map<DataCategoryType, std::string> dataCategoryTypeNameMap = {
    { DataCategoryType::Signal, "SIGNAL" }, { DataCategoryType::Background, "BACKGROUND" },
    { DataCategoryType::Data, "DATA" }, { DataCategoryType::DYJets, "DY_JETS" },
    { DataCategoryType::DYJets_incl, "DY_JETS_incl" },{ DataCategoryType::DYJets_excl, "DY_JETS_excl" },
    { DataCategoryType::ZL, "ZL" }, { DataCategoryType::ZJ, "ZJ" }, { DataCategoryType::ZTT, "ZTT" },
    { DataCategoryType::ZTT_MC, "ZTT_MC" }, { DataCategoryType::ZTT_L, "ZTT_L" }, { DataCategoryType::Embedded, "EMBEDDED" },
    { DataCategoryType::Limits, "LIMITS" }, { DataCategoryType::Composit, "COMPOSIT" }, { DataCategoryType::QCD, "QCD" },
    { DataCategoryType::WJets, "W_JETS" }, { DataCategoryType::WJets_MC, "W_JETS_MC" },
    { DataCategoryType::WJets_MC_incl, "W_JETS_MC_incl" }, { DataCategoryType::WJets_MC_excl, "W_JETS_MC_excl" }
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
    std::set<Channel> channels;
    std::set<std::string> sub_categories;
    DataSourceScaleFactorMap sources_sf;

    DataCategory()
        : color(kBlack), limits_sf(1.0), draw(false), draw_sf(1) {}

    bool IsSignal() const { return types.count(DataCategoryType::Signal); }
    bool IsBackground() const { return types.count(DataCategoryType::Background); }
    bool IsData() const { return types.count(DataCategoryType::Data); }
    bool IsComposit() const { return types.count(DataCategoryType::Composit); }
};

typedef std::map<std::string, DataCategory> DataCategoryMap;
typedef std::set<const DataCategory*> DataCategoryPtrSet;
typedef std::vector<const DataCategory*> DataCategoryPtrVector;
typedef std::map<DataCategoryType, DataCategoryPtrSet> DataCategoryTypeMap;

class DataCategoryCollection {
public:
    DataCategoryCollection(const std::string& sources_cfg_name, const std::string& signal_list, Channel channel_id)
    {
        DataCategory category;
        std::ifstream cfg(sources_cfg_name);
        size_t line_number = 0;
        while(ReadNextCategory(cfg, line_number, category)) {
            CheckCategoryValidity(category);
            if(category.channels.size() && !category.channels.count(channel_id)) continue;
            categories[category.name] = category;
            all_categories.push_back(&categories[category.name]);
            for(DataCategoryType type : category.types)
                categories_by_type[type].insert(&categories[category.name]);
            for(const auto& source_entry : category.sources_sf)
                all_sources.insert(source_entry.first);
        }
        const auto& signal_names = ParseSignalList(signal_list);
        for(const auto& signal_name : signal_names) {
            if(!categories.count(signal_name))
                throw exception("Undefined signal '") << signal_name << "'.";
            categories[signal_name].draw = true;
        }
    }

    const DataCategoryPtrVector& GetAllCategories() const { return all_categories; }
    const DataCategoryPtrSet& GetCategories(DataCategoryType dataCategoryType) const
    {
        if(!categories_by_type.count(dataCategoryType))
            return empty_category_set;
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
    void CheckCategoryValidity(const DataCategory& category) const
    {
        if(categories.count(category.name))
            throw exception("Category with name '") << category.name << "' is already defined.";
        if(category.sub_categories.size() && !category.types.count(DataCategoryType::Composit))
            throw exception("Not composit category '") << category.name << "' may not contain sub-categories.";
        if(category.types.count(DataCategoryType::Composit) && category.sources_sf.size())
            throw exception("Composit category '") << category.name << "' may not contain direct file definitions.";
        for(const auto& sub_category : category.sub_categories) {
            if(!categories.count(sub_category))
                throw exception("Sub-category '") << sub_category << "' for category '"
                                                  << category.name << "' is not defined.";
            if(categories.at(sub_category).types.count(DataCategoryType::Composit))
                throw exception("Invalid sub-category '") << sub_category << "' for category '" << category.name
                                                      << "'. Composit category hierarchy is not supported.";
        }
        for(const auto& source_entry : category.sources_sf) {
            if(all_sources.count(source_entry.first))
                throw exception("Source '") << source_entry.first << "' is already part of the other data category.";
        }
    }

    static bool ReadNextCategory(std::istream& cfg, size_t& line_number, DataCategory& category)
    {
        category = DataCategory();
        bool category_started = false;
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            ++line_number;
            if ((cfgLine.size() && cfgLine.at(0) == '#') || (!cfgLine.size() && !category_started)) continue;
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
        static const char separator = ':';

        const size_t pos = cfgLine.find(separator);
        if(pos == std::string::npos)
            throw exception("bad source config syntax for a parameter in line ") << line_number;
        const std::string param_name = cfgLine.substr(0, pos);
        if(pos + 2 >= cfgLine.size())
            throw exception("empty parameter value in source config in line ") << line_number;
        const std::string param_value = cfgLine.substr(pos + 2);
        std::istringstream ss(param_value);
        ss >> std::boolalpha;
        if(param_name == "type") {
            DataCategoryType type;
            ss >> type;
            category.types.insert(type);
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
        } else if(param_name == "draw") {
            ss >> category.draw;
        } else if(param_name == "channel") {
            Channel channel_id;
            ss >> channel_id;
            category.channels.insert(channel_id);
        } else if(param_name == "datacard") {
            ss >> category.datacard;
        } else if(param_name == "subcategory") {
            category.sub_categories.insert(param_value);
        } else
            throw exception("Unsupported parameter '") << param_name << "' in configuration line " << line_number;
    }

    static std::set<std::string> ParseSignalList(const std::string& signal_list)
    {
        static const char separator = ',';

        std::set<std::string> result;
        size_t prev_pos = 0;
        for(bool next = true; next;) {
            const size_t pos = signal_list.find(separator, prev_pos);
            next = pos != std::string::npos;
            const size_t last_pos = next ? pos - 1 : std::string::npos;
            const std::string signal_name = signal_list.substr(prev_pos, last_pos);
            result.insert(signal_name);
            prev_pos = pos + 1;
        }
        return result;
    }

private:
    std::set<std::string> all_sources;
    DataCategoryMap categories;
    DataCategoryPtrVector all_categories;
    DataCategoryTypeMap categories_by_type;
    DataCategoryPtrSet empty_category_set;
};

std::ostream& operator<<(std::ostream& s, const DataCategory& category){
    s << "Name: " << category.name << ", Title: '" << category.title << "', Color: " << category.color << std::endl;
    for(const auto& source : category.sources_sf)
        s << "File: " << source.first << ", SF: " << source.second << "\n";
    return s;
}

enum class EventRegion { Unknown = 0, OS_Isolated = 1, OS_NotIsolated = 2, SS_Isolated = 3, SS_NotIsolated = 4,
                         OS_Iso_HighMt = 5, SS_Iso_HighMt = 6, OS_NotIso_HighMt = 7, SS_NotIso_HighMt = 8 };
enum class EventCategory { Inclusive = 0, OneJet_ZeroBtag = 1, OneJet_OneBtag = 2, TwoJets_ZeroBtag = 3,
                           TwoJets_OneBtag = 4, TwoJets_TwoBtag = 5 };

namespace detail {
static const std::map<EventCategory, std::string> eventCategoryNamesMap =
          { { EventCategory::Inclusive, "Inclusive" }, { EventCategory::OneJet_ZeroBtag, "1jet0btag" },
            { EventCategory::OneJet_OneBtag, "1jet1btag" }, { EventCategory::TwoJets_ZeroBtag, "2jets0btag" },
          { EventCategory::TwoJets_OneBtag, "2jets1btag"}, { EventCategory::TwoJets_TwoBtag, "2jets2btag" } };

static const std::map<EventRegion, std::string> eventRegionNamesMap =
          { { EventRegion::Unknown, "Unknown"}, { EventRegion::OS_Isolated, "OS_Isolated"},
            { EventRegion::OS_NotIsolated, "OS_NotIsolated"}, { EventRegion::SS_Isolated, "SS_Isolated"},
            { EventRegion::SS_NotIsolated, "SS_NotIsolated"}, { EventRegion::OS_Iso_HighMt, "OS_Iso_HighMt"},
            { EventRegion::SS_Iso_HighMt, "SS_Iso_HighMt"} , { EventRegion::OS_NotIso_HighMt, "OS_NotIso_HighMt"},
            { EventRegion::SS_NotIso_HighMt, "SS_NotIso_HighMt"} };
} // namespace detail

typedef std::vector<EventCategory> EventCategoryVector;
typedef std::set<EventCategory> EventCategorySet;

static const EventCategorySet AllEventCategories = tools::collect_map_keys(detail::eventCategoryNamesMap);
static const EventCategorySet OneJetEventCategories = { EventCategory::OneJet_ZeroBtag, EventCategory::OneJet_OneBtag };
static const EventCategorySet TwoJetsEventCategories =
        { EventCategory::TwoJets_ZeroBtag, EventCategory::TwoJets_OneBtag, EventCategory::TwoJets_TwoBtag };

typedef std::set<EventRegion> EventRegionSet;
typedef std::map<EventRegion, EventRegion> EventRegionMap;
static const EventRegionMap HighMt_LowMt_RegionMap = { {EventRegion::OS_Iso_HighMt, EventRegion::OS_Isolated},
                                                     {EventRegion::SS_Iso_HighMt, EventRegion::SS_Isolated},
                                                       {EventRegion::OS_NotIso_HighMt, EventRegion::OS_NotIsolated},
                                                       {EventRegion::SS_NotIso_HighMt, EventRegion::SS_NotIsolated} };
static const EventRegionSet HighMtRegions = {EventRegion::OS_Iso_HighMt, EventRegion::SS_Iso_HighMt,
                                             EventRegion::OS_NotIso_HighMt, EventRegion::SS_NotIso_HighMt};

static const EventRegionSet QcdRegions = {EventRegion::OS_Isolated, EventRegion::SS_Isolated,
                                             EventRegion::OS_NotIsolated, EventRegion::SS_NotIsolated};

static const EventRegionSet AllEventRegions = tools::collect_map_keys(detail::eventRegionNamesMap);

std::ostream& operator<<(std::ostream& s, const EventCategory& eventCategory) {
    s << detail::eventCategoryNamesMap.at(eventCategory);
    return s;
}

std::wostream& operator<<(std::wostream& s, const EventCategory& eventCategory) {
    const std::string str = detail::eventCategoryNamesMap.at(eventCategory);
    s << std::wstring(str.begin(), str.end());
    return s;
}

std::ostream& operator<<(std::ostream& s, const EventRegion& eventRegion) {
    s << detail::eventRegionNamesMap.at(eventRegion);
    return s;
}


EventCategoryVector DetermineEventCategories(const std::vector<float>& csv_Bjets, double CSVM, double CSVT)
{
    EventCategoryVector categories;
    categories.push_back(EventCategory::Inclusive);

    typedef std::pair<size_t, size_t> Key;

    static const std::map< Key, EventCategory> categories_map {
        { { 1, 0 }, EventCategory::OneJet_ZeroBtag }, { { 1, 1 }, EventCategory::OneJet_OneBtag },
        { { 2, 0 }, EventCategory::TwoJets_ZeroBtag }, { { 2, 1 }, EventCategory::TwoJets_OneBtag },
        { { 2, 2 }, EventCategory::TwoJets_TwoBtag }
    };

    const size_t n_jets = std::min<size_t>(csv_Bjets.size(), 2);
    size_t n_btag = 0;
    if(n_jets == 1 && csv_Bjets.at(0) > CSVT) ++n_btag;
    if(n_jets == 2 && csv_Bjets.at(0) > CSVM) ++n_btag;
    if(n_jets == 2 && csv_Bjets.at(1) > CSVM) ++n_btag;

    const Key key(n_jets, n_btag);
    if(categories_map.count(key))
        categories.push_back(categories_map.at(key));

    return categories;
}

} // namespace analysis
