#pragma once

#include <vector>
#include <map>

namespace MVA_Selections {

enum EventCategory { UnknownCategory, OneJet_ZeroBtag, OneJet_OneBtag, TwoJets_ZeroBtag, TwoJets_OneBtag, TwoJets_TwoBtag };
inline EventCategory EventCategoryFromString(const std::string& category_name)
{
    static std::map<std::string, EventCategory> category_name_map;
    if(!category_name_map.size()) {
        category_name_map["1jet0btag"] = OneJet_ZeroBtag;
        category_name_map["1jet1btag"] = OneJet_OneBtag;
        category_name_map["2jets0btag"] = TwoJets_ZeroBtag;
        category_name_map["2jets1btag"] = TwoJets_OneBtag;
        category_name_map["2jets2btag"] = TwoJets_TwoBtag;
    }
    if(!category_name_map.count(category_name))
        throw std::runtime_error("Unknown category name");
    return category_name_map[category_name];
}

enum Channel { ETau, MuTau, TauTau };

typedef std::map<std::string, size_t> var_map;
typedef std::vector<std::string> str_vector;
typedef std::pair<Channel, EventCategory> ChannelCategoryPair;
typedef std::map<ChannelCategoryPair, str_vector> var_list;
typedef std::map<ChannelCategoryPair, var_map> var_map_list;

inline var_map MakeVarMap(const str_vector& vars)
{
    var_map result;
    for(size_t n = 0; n < vars.size(); ++n)
        result[vars[n]] = n;
    return result;
}

const str_vector& Input_Variables(Channel channel, EventCategory category)
{
    static var_list l;
    if(!l.size()) {
        {
            str_vector& v = l[ChannelCategoryPair(MuTau, TwoJets_ZeroBtag)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
            v.push_back("pt_b2");
            v.push_back("DR_bb");
            v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            v.push_back("DR_HBBHTT");
            v.push_back("Pt_Hbb");
            v.push_back("DeltaPhi_METTT");
            v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ChannelCategoryPair(MuTau, TwoJets_OneBtag)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
            v.push_back("pt_b2");
            v.push_back("DR_bb");
            v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            v.push_back("Pt_Htt");
            v.push_back("DR_HBBHTT");
            v.push_back("Pt_Hbb");
            v.push_back("DeltaPhi_METTT");
            v.push_back("PtH");
            v.push_back("mT2");
            v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
        {
            str_vector& v = l[ChannelCategoryPair(MuTau, TwoJets_TwoBtag)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            //v.push_back("pt_b1");
            v.push_back("pt_b2");
            v.push_back("DR_bb");
            //v.push_back("DPhi_BBMET");
            v.push_back("DR_ll");
            //v.push_back("Pt_Htt");
            //v.push_back("DR_HBBHTT");
            //v.push_back("Pt_Hbb");
            //v.push_back("DeltaPhi_METTT");
            //v.push_back("PtH");
            v.push_back("mT2");
            //v.push_back("mT1");
            v.push_back("Pt_Htt_MET");
        }
    }
    const ChannelCategoryPair key(channel, category);
    if(!l.count(key))
        throw std::runtime_error("Unknown combination channel-category");
    return l[key];
}

const var_map& Input_Variables_Map(Channel channel, EventCategory category)
{
    static var_map_list l;
    const ChannelCategoryPair key(channel, category);
    if(!l.count(key))
        l[key] = MakeVarMap(Input_Variables(channel, category));
    return l[key];
}

} // namespace MVA_Selections
