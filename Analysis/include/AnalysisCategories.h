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

#include <TFile.h>
#include <Rtypes.h>

namespace analysis {

struct DataSource {
    std::string file_name;
    double scale_factor;
    TFile* file;
    std::shared_ptr<ntuple::FlatTree> tree;
};

typedef std::vector<DataSource> DataSourceVector;

struct DataCategory;
typedef std::list<DataCategory> DataCategoryCollection;

struct DataCategory {
    std::string name;
    std::string title;
    EColor color;

    DataSourceVector sources;

public:
    bool IsData() const { return NameContains("DATA"); }
    bool IsSignal() const { return NameContains("SIGNAL"); }
    bool IsReference() const { return NameContains("REFERENCE"); }
    bool IsVirtual() const { return NameContains("VIRTUAL"); }
    bool IsForLimitsOnly() const { return NameContains("LIMITS"); }
    bool IsSumBkg() const { return NameContains("SUM"); }

    bool NameContains(const std::string& substring) const { return name.find(substring) != std::string::npos; }

public:

    static DataCategoryCollection ReadFromFile(const std::string& cfg_name, const std::string& dataName)
    {
        DataCategoryCollection categories;
        std::ifstream cfg(cfg_name);
        std::shared_ptr<DataCategory> currentCategory;
        while (cfg.good()) {
            std::string cfgLine;
            std::getline(cfg,cfgLine);
            if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
            if (cfgLine.at(0) == '[') {
                if(currentCategory)
                    categories.push_back(*currentCategory);
                currentCategory = std::shared_ptr<DataCategory>(new DataCategory());
                const size_t pos = cfgLine.find(']');
                currentCategory->name = cfgLine.substr(1, pos - 1);
                std::getline(cfg, currentCategory->title);
                std::string colorLine;
                std::getline(cfg,colorLine);
                std::istringstream ss(colorLine);
                ss >> currentCategory->color;
            }
            else if (currentCategory) {
                std::istringstream ss(cfgLine);
                DataSource source;
                ss >> source.file_name;
                ss >> source.scale_factor;
                currentCategory->sources.push_back(source);
            }
            else
                throw std::runtime_error("bad source file format");
          }
        if(currentCategory)
            categories.push_back(*currentCategory);
        DataCategoryCollection filteredCategories;
        for(const DataCategory& category : categories) {
    //            if(category.IsSignal()) {
    //                const size_t sub_name_pos = category.name.find(' ');
    //                const std::string sub_name = category.name.substr(sub_name_pos + 1);
    //                if(sub_name != signalName)
    //                    continue;
    //            }
            if(category.IsData()) {
                const size_t sub_name_pos = category.name.find(' ');
                const std::string sub_name = category.name.substr(sub_name_pos + 1);
                if(sub_name != dataName)
                    continue;
            }
            filteredCategories.push_back(category);
        }
        return filteredCategories;
    }
};

std::ostream& operator<<(std::ostream& s, const DataSource& source){
    s << "File: " << source.file_name << ", SF: " << source.scale_factor;
    return s;
}

std::ostream& operator<<(std::ostream& s, const DataCategory& category){
    s << "Name: " << category.name << ", Title: '" << category.title << "', Color: " << category.color << std::endl;
    for(const DataSource& source : category.sources)
        s << source << std::endl;
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
