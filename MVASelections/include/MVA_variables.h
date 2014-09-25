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

enum MvaMethod { BDT, BDTMitFisher, BDTD };
inline MvaMethod MvaMethodFromString(const std::string& method_name)
{
    static std::map<std::string, MvaMethod> method_name_map;
    if(!method_name_map.size()) {
        method_name_map["BDT"] = BDT;
        method_name_map["BDTMitFisher"] = BDTMitFisher;
        method_name_map["BDTD"] = BDTD;
    }
    if(!method_name_map.count(method_name))
        throw std::runtime_error("Unknown MVA method");
    return method_name_map[method_name];
}

struct ParamId {
    Channel channel;
    EventCategory category;
    MvaMethod method;

    ParamId(Channel _channel, EventCategory _category, MvaMethod _method)
        : channel(_channel), category(_category), method(_method) {}

    bool operator< (const ParamId& other) const
    {
        if(channel < other.channel) return true;
        if(channel > other.channel) return false;
        if(category < other.category) return true;
        if(category > other.category) return false;
        return method < other.method;
    }
};


typedef std::map<std::string, size_t> var_map;
typedef std::vector<std::string> str_vector;
typedef std::map<ParamId, str_vector> var_list;
typedef std::map<ParamId, var_map> var_map_list;

inline var_map MakeVarMap(const str_vector& vars)
{
    var_map result;
    for(size_t n = 0; n < vars.size(); ++n)
        result[vars[n]] = n;
    return result;
}

const str_vector& Input_Variables(Channel channel, EventCategory category, MvaMethod method)
{
    static var_list l;
    if(!l.size()) {
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_ZeroBtag, BDT)];
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
            str_vector& v = l[ParamId(MuTau, TwoJets_ZeroBtag, BDTMitFisher)];
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
            str_vector& v = l[ParamId(MuTau, TwoJets_ZeroBtag, BDTD)];
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
            str_vector& v = l[ParamId(MuTau, TwoJets_OneBtag, BDT)];
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
            str_vector& v = l[ParamId(MuTau, TwoJets_OneBtag, BDTMitFisher)];
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
            str_vector& v = l[ParamId(MuTau, TwoJets_OneBtag, BDTD)];
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
            str_vector& v = l[ParamId(MuTau, TwoJets_TwoBtag, BDT)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
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
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_TwoBtag, BDTMitFisher)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
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
        {
            str_vector& v = l[ParamId(MuTau, TwoJets_TwoBtag, BDTD)];
            v.push_back("pt_mu");
            v.push_back("pt_tau");
            v.push_back("pt_b1");
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
            //v.push_back("Pt_Htt_MET");
        }
    }
    const ParamId key(channel, category, method);
    if(!l.count(key))
        throw std::runtime_error("Unknown combination channel-category");
    return l[key];
}

const var_map& Input_Variables_Map(Channel channel, EventCategory category, MvaMethod method)
{
    static var_map_list l;
    const ParamId key(channel, category, method);
    if(!l.count(key))
        l[key] = MakeVarMap(Input_Variables(channel, category, method));
    return l[key];
}

} // namespace MVA_Selections
