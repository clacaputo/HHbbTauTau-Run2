/*!
 * \file FlatAnalyzerDataCollection.h
 * \brief Collection of histogram containers for flat tree analyzers.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2015-02-07 created
 *
 * Copyright 2015 Konstantin Androsov <konstantin.androsov@gmail.com>,
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

#include "FlatAnalyzerData.h"

namespace analysis {

struct FlatAnalyzerDataId;
typedef std::set<FlatAnalyzerDataId> FlatAnalyzerDataIdSet;

struct FlatAnalyzerDataId {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;
    std::string dataCategoryName;

    FlatAnalyzerDataId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventRegion(EventRegion::OS_Isolated), eventEnergyScale(EventEnergyScale::Central)
    {}

    FlatAnalyzerDataId(EventCategory _eventCategory, EventSubCategory _eventSubCategory, EventRegion _eventRegion,
                       EventEnergyScale _eventEnergyScale, const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale), dataCategoryName(_dataCategoryName)
    {}

    bool operator< (const FlatAnalyzerDataId& other) const
    {
        if(eventCategory < other.eventCategory) return true;
        if(eventCategory > other.eventCategory) return false;
        if(eventSubCategory < other.eventSubCategory) return true;
        if(eventSubCategory > other.eventSubCategory) return false;
        if(eventRegion < other.eventRegion) return true;
        if(eventRegion > other.eventRegion) return false;
        if(eventEnergyScale < other.eventEnergyScale) return true;
        if(eventEnergyScale > other.eventEnergyScale) return false;
        return dataCategoryName < other.dataCategoryName;
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        std::ostringstream ss;
        ss << eventCategory << separator << eventSubCategory << separator << eventRegion << separator
           << eventEnergyScale << separator << dataCategoryName;
        return ss.str();
    }

};

template<typename ...Types>
struct FlatAnalyzerDataMetaId;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string> {
    EventCategory eventCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;
    std::string dataCategoryName;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventRegion(EventRegion::OS_Isolated),
          eventEnergyScale(EventEnergyScale::Central) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventRegion _eventRegion, EventEnergyScale _eventEnergyScale,
                           const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale), dataCategoryName(_dataCategoryName)  {}

    FlatAnalyzerDataId MakeId(EventSubCategory eventSubCategory) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        static const std::string wildcard = "*";
        std::ostringstream ss;
        ss << eventCategory << separator << wildcard << separator << eventRegion << separator
           << eventEnergyScale << separator << dataCategoryName;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventRegion, EventEnergyScale, std::string> FlatAnalyzerDataMetaId_noSub;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventRegion, std::string> {
    EventCategory eventCategory;
    EventRegion eventRegion;
    std::string dataCategoryName;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventRegion(EventRegion::OS_Isolated) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventRegion _eventRegion, const std::string& _dataCategoryName)
        : eventCategory(_eventCategory), eventRegion(_eventRegion), dataCategoryName(_dataCategoryName)  {}

    FlatAnalyzerDataId MakeId(EventSubCategory eventSubCategory, EventEnergyScale eventEnergyScale) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    FlatAnalyzerDataMetaId_noSub MakeMetaId(EventEnergyScale eventEnergyScale) const
    {
        return FlatAnalyzerDataMetaId_noSub(eventCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        static const std::string wildcard = "*";
        std::ostringstream ss;
        ss << eventCategory << separator << wildcard << separator << eventRegion << separator
           << wildcard << separator << dataCategoryName;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventRegion, std::string> FlatAnalyzerDataMetaId_noSub_noES;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion, EventEnergyScale> {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventRegion eventRegion;
    EventEnergyScale eventEnergyScale;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventRegion(EventRegion::OS_Isolated), eventEnergyScale(EventEnergyScale::Central) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventSubCategory _eventSubCategory, EventRegion _eventRegion,
                           EventEnergyScale _eventEnergyScale)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventRegion(_eventRegion),
          eventEnergyScale(_eventEnergyScale) {}

    FlatAnalyzerDataId MakeId(const std::string& dataCategoryName) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        std::ostringstream ss;
        ss << eventCategory << separator << eventSubCategory << separator << eventRegion << separator
           << eventEnergyScale;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventRegion, EventEnergyScale>
    FlatAnalyzerDataMetaId_noName;

template<>
struct FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale> {
    EventCategory eventCategory;
    EventSubCategory eventSubCategory;
    EventEnergyScale eventEnergyScale;

    FlatAnalyzerDataMetaId()
        : eventCategory(EventCategory::Inclusive), eventSubCategory(EventSubCategory::NoCuts),
          eventEnergyScale(EventEnergyScale::Central) {}

    FlatAnalyzerDataMetaId(EventCategory _eventCategory, EventSubCategory _eventSubCategory,
                           EventEnergyScale _eventEnergyScale)
        : eventCategory(_eventCategory), eventSubCategory(_eventSubCategory), eventEnergyScale(_eventEnergyScale) {}

    FlatAnalyzerDataId MakeId(EventRegion eventRegion, const std::string& dataCategoryName) const
    {
        return FlatAnalyzerDataId(eventCategory, eventSubCategory, eventRegion, eventEnergyScale, dataCategoryName);
    }

    FlatAnalyzerDataMetaId_noName MakeMetaId(EventRegion eventRegion) const
    {
        return FlatAnalyzerDataMetaId_noName(eventCategory, eventSubCategory, eventRegion, eventEnergyScale);
    }

    std::string GetName() const
    {
        static const std::string separator = "/";
        static const std::string wildcard = "*";
        std::ostringstream ss;
        ss << eventCategory << separator << eventSubCategory << separator << wildcard << separator
           << eventEnergyScale;
        return ss.str();
    }
};

typedef FlatAnalyzerDataMetaId<EventCategory, EventSubCategory, EventEnergyScale>
    FlatAnalyzerDataMetaId_noRegion_noName;

template<typename ...Types>
std::ostream& operator<< (std::ostream& s, const FlatAnalyzerDataMetaId<Types...>& meta_id)
{
    s << meta_id.GetName();
    return s;
}

typedef std::shared_ptr<FlatAnalyzerData> FlatAnalyzerDataPtr;
typedef std::map<FlatAnalyzerDataId, FlatAnalyzerDataPtr> FlatAnalyzerDataMap;

class FlatAnalyzerDataCollection {
public:
    FlatAnalyzerDataCollection(const std::string& outputFileName, bool store)
    {
        if(store)
            outputFile = std::shared_ptr<TFile>(root_ext::AnalyzerData::CreateFile(outputFileName));
    }

    FlatAnalyzerData& Get(const FlatAnalyzerDataId& id, Channel channel)
    {
        auto& anaData = anaDataMap[id];
        if(!anaData)
            anaData = MakeAnaData(id, channel);
        return *anaData;
    }

    void Fill(const FlatAnalyzerDataId& id, Channel channel, const FlatEventInfo& eventInfo, double weight)
    {
        auto& anaData = Get(id, channel);
        anaData.Fill(eventInfo, weight);
    }

    void FillSubCategories(const FlatAnalyzerDataMetaId_noSub& meta_id, Channel channel,
                           const FlatEventInfo& eventInfo, double weight)
    {
        using namespace cuts::massWindow;

        const double mass_tautau = eventInfo.event->m_sv_MC;
        Fill(meta_id.MakeId(EventSubCategory::NoCuts), channel, eventInfo, weight);
        const bool inside_mass_window = mass_tautau > m_tautau_low && mass_tautau < m_tautau_high
                && eventInfo.Hbb.M() > m_bb_low && eventInfo.Hbb.M() < m_bb_high;
        if(inside_mass_window)
            Fill(meta_id.MakeId(EventSubCategory::MassWindow), channel, eventInfo, weight);
        else
            Fill(meta_id.MakeId(EventSubCategory::OutsideMassWindow), channel, eventInfo, weight);
        if(eventInfo.fitResults.has_valid_mass)
            Fill(meta_id.MakeId(EventSubCategory::KinematicFitConverged), channel, eventInfo, weight);
        if(eventInfo.fitResults.has_valid_mass && inside_mass_window)
            Fill(meta_id.MakeId(EventSubCategory::KinematicFitConvergedWithMassWindow), channel, eventInfo, weight);
        if(eventInfo.fitResults.has_valid_mass && !inside_mass_window)
            Fill(meta_id.MakeId(EventSubCategory::KinematicFitConvergedOutsideMassWindow), channel, eventInfo, weight);
    }

    void FillEnergyScales(const FlatAnalyzerDataMetaId_noSub_noES& meta_id, Channel channel,
                          const FlatEventInfo& eventInfo, double weight, const std::set<EventEnergyScale>& energyScales)
    {
        for (EventEnergyScale energyScale : energyScales)
            FillSubCategories(meta_id.MakeMetaId(energyScale), channel, eventInfo, weight);
    }

    void FillCentralAndJetRelatedScales(const FlatAnalyzerDataMetaId_noSub_noES& meta_id, Channel channel,
                                        const FlatEventInfo& eventInfo, double weight)
    {
        static const std::set<EventEnergyScale> energyScales =
            { EventEnergyScale::Central, EventEnergyScale::JetUp, EventEnergyScale::JetDown,
              EventEnergyScale::BtagEfficiencyUp, EventEnergyScale::BtagEfficiencyDown,
              EventEnergyScale::BtagFakeUp, EventEnergyScale::BtagFakeDown };
        FillEnergyScales(meta_id, channel, eventInfo, weight, energyScales);
    }

    void FillAllEnergyScales(const FlatAnalyzerDataMetaId_noSub_noES& meta_id, Channel channel,
                             const FlatEventInfo& eventInfo, double weight)
    {
        FillEnergyScales(meta_id, channel, eventInfo, weight, AllEventEnergyScales);
    }

private:
    FlatAnalyzerDataPtr MakeAnaData(const FlatAnalyzerDataId& id, Channel channel) const
    {
        if(channel == Channel::ETau || channel == Channel::MuTau) {
            if(id.eventCategory == EventCategory::TwoJets_TwoBtag
                    || id.eventCategory == EventCategory::TwoJets_TwoLooseBtag)
                return Make<FlatAnalyzerData_semileptonic_2tag>(id);
            return Make<FlatAnalyzerData_semileptonic>(id);
        }

        if(channel == Channel::TauTau) {
            if(id.eventCategory == analysis::EventCategory::TwoJets_TwoBtag
                    || id.eventCategory == analysis::EventCategory::TwoJets_TwoLooseBtag)
                return Make<FlatAnalyzerData_tautau_2tag>(id);
            return Make<FlatAnalyzerData_tautau>(id);
        }

        throw exception("Can't make analyzer data for ") << channel << " channel.";
    }

    template<typename AnaData>
    FlatAnalyzerDataPtr Make(const FlatAnalyzerDataId& id) const
    {
        const bool fill_all = id.eventEnergyScale == EventEnergyScale::Central;
        if(outputFile) {
            const std::string dir_name = id.GetName();
            return FlatAnalyzerDataPtr(new AnaData(*outputFile, dir_name, fill_all));
        }
        return FlatAnalyzerDataPtr(new AnaData(fill_all));
    }

private:
    std::shared_ptr<TFile> outputFile;
    FlatAnalyzerDataMap anaDataMap;
};

} // namespace analysis
