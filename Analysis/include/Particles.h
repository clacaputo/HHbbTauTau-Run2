/*!
 * \file Particles.h
 * \brief Enumaration of the PDG MC particle codes.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2013-04-05 created
 */

#pragma once

#include <cmath>
#include <string>
#include <map>
#include <stdexcept>
#include <sstream>

#include <iostream>

namespace particles {
class ParticleCode {
private:
    typedef std::map<int, std::string> CodeToNameMap;
    typedef std::map<std::string, int> NameToCodeMap;
    static CodeToNameMap& Names()
    {
        static CodeToNameMap names;
        return names;
    }

    static NameToCodeMap& Codes()
    {
        static NameToCodeMap codes;
        return codes;
    }

    static bool& Freeze()
    {
        static bool freeze = false;
        return freeze;
    }

public:
    ParticleCode() : raw_code(0) {}
    ParticleCode(int _raw_code)
        : raw_code(_raw_code) {
        if(Names().find(raw_code) == Names().end()) {
            std::ostringstream ss;
            ss << "Unknown particle id " << raw_code << ".";
            std::cerr << ss.str() << std::endl;
            //throw std::runtime_error(ss.str());
            raw_code = 0;
        }
    }
    ParticleCode(const std::string& name)
    {
        const NameToCodeMap::const_iterator iter = Codes().find(name);
        if(iter == Codes().end()) {
            std::ostringstream ss;
            ss << "Unknown particle name '" << name << "'.";
            throw std::runtime_error(ss.str());
        }
        raw_code = iter->second;
    }
    ParticleCode(int _raw_code, const std::string& name)
        : raw_code(_raw_code) {
        if(Freeze())
            throw std::runtime_error("New particle code can not be defined.");
        if(Names().find(raw_code) != Names().end()) {
            std::ostringstream ss;
            ss << "Particle id " << raw_code << " is already defined.";
            throw std::runtime_error(ss.str());
        }
        if(Codes().find(name) != Codes().end()) {
            std::ostringstream ss;
            ss << "Particle name '" << name << "' is already taken.";
            throw std::runtime_error(ss.str());
        }

        Names()[raw_code] = name;
        Codes()[name] = raw_code;
        if(!raw_code)
            Freeze() = true;
    }

    int RawCode() const { return raw_code; }
    const std::string& Name() const { return Names().find(raw_code)->second; }
    bool operator<(const ParticleCode& other) const { return raw_code < other.raw_code; }
    bool operator>(const ParticleCode& other) const { return raw_code > other.raw_code; }
    bool operator==(const ParticleCode& other) const { return raw_code == other.raw_code; }
    bool operator!=(const ParticleCode& other) const { return raw_code != other.raw_code; }
private:
    int raw_code;
};

#define PARTICLE(name, code) \
    static const ParticleCode name(code, #name)

PARTICLE(d, 1);
PARTICLE(u, 2);
PARTICLE(s, 3);
PARTICLE(c, 4);
PARTICLE(b, 5);
PARTICLE(t, 6);
PARTICLE(b_prime, 7);
PARTICLE(t_prime, 8);
PARTICLE(e, 11);
PARTICLE(nu_e, 12);
PARTICLE(mu, 13);
PARTICLE(nu_mu, 14);
PARTICLE(tau, 15);
PARTICLE(nu_tau, 16);
PARTICLE(tau_prime, 17);
PARTICLE(nu_tau_prime, 18);
PARTICLE(g, 21);
PARTICLE(gamma, 22);
PARTICLE(pi_zero, 111);
PARTICLE(rho_zero, 113);
PARTICLE(K_zero_L, 130);
PARTICLE(pi, 211);
PARTICLE(rho, 213);
PARTICLE(eta, 221);
PARTICLE(omega, 223);
PARTICLE(K_zero_S, 310);
PARTICLE(K_zero, 311);
PARTICLE(K, 321);
PARTICLE(K_star, 323);
PARTICLE(p, 2212);
PARTICLE(n, 2112);
PARTICLE(uu1, 2203);
PARTICLE(sigma_star_0, 3214);
PARTICLE(sigma_minus, 3112);
PARTICLE(MC_internal_92, 92);
PARTICLE(Delta_plusplus, 2224);
PARTICLE(Delta_plus, 2214);
PARTICLE(Delta_zero, 2114);
PARTICLE(Delta_minus, 1114);
PARTICLE(Higgs, 25);
PARTICLE(K_star_0, 313);
PARTICLE(B_0, 511);
PARTICLE(Lambda, 3122);
PARTICLE(eta_prime, 331);
PARTICLE(ud0, 2101);
PARTICLE(xi_minus, 3312);
PARTICLE(phi, 333);
PARTICLE(sigma_plus, 3222);
PARTICLE(sigma_0, 3212);
PARTICLE(ud1, 2103);
PARTICLE(D_star_0, 423);
PARTICLE(D_0, 421);
PARTICLE(sigma_star_minus, 3114);
PARTICLE(sigma_star_plus, 3224);
PARTICLE(D_plus, 411);
PARTICLE(a1_plus, 20213);
PARTICLE(D_star_plus, 413);
PARTICLE(B_star_0, 513);
PARTICLE(xi_0, 3322);
PARTICLE(B_star_plus, 523);
PARTICLE(Lambda_b_zero, 5122);
PARTICLE(Lambda_c_plus, 4122);
PARTICLE(D_s_plus, 431);
PARTICLE(D_s_star_plus, 433);
PARTICLE(xi_star_0, 3324);
PARTICLE(omega_minus, 3334);
PARTICLE(xi_star_minus, 3314);
PARTICLE(sigma_c_star_0, 4114);
PARTICLE(B_plus, 521);
PARTICLE(MC_internal_91, 91);
PARTICLE(sigma_c_star_plus, 4214);
PARTICLE(sigma_c_plus, 4212);
PARTICLE(xi_c_prime_0, 4312);
PARTICLE(D_s2_star_plus, 435);
PARTICLE(sigma_c_star_plusplus, 4224);
PARTICLE(K1_plus, 10323);
PARTICLE(K1_0, 20313);
PARTICLE(B_s_star_0, 533);
PARTICLE(B_s_0, 531);
PARTICLE(sigma_b_star_plus, 5224);
PARTICLE(xi_c_zero, 4132);
PARTICLE(xi_b_minus, 5132);
PARTICLE(sigma_b_0, 5212);
PARTICLE(sigma_b_star_minus, 5114);
PARTICLE(sigma_c_plusplus, 4222);
PARTICLE(xi_b_0, 5232);
PARTICLE(xi_b_star_0, 5324);
PARTICLE(D1_0_H, 20423);
PARTICLE(D2_star_plus, 415);
PARTICLE(xi_c_plus, 4232);
PARTICLE(D1_plus_H, 20413);
PARTICLE(sigma_b_star_0, 5214);
PARTICLE(xi_c_star_zero, 4314);
PARTICLE(D0_star_plus, 10411);
PARTICLE(D1_plus, 10413);
PARTICLE(D1_0, 10423);
PARTICLE(J_psi, 443);
PARTICLE(D_s0_star_plus, 10431);
PARTICLE(sigma_c_0, 4112);
PARTICLE(D2_star_0, 425);
PARTICLE(D_s1_plus_2536, 10433);
PARTICLE(f0, 10221);
PARTICLE(sigma_b_plus, 5222);
PARTICLE(xi_c_star_plus, 4324);
PARTICLE(p_diffr_plus, 9902210);
PARTICLE(deuteron, 1000010020);
PARTICLE(eta_c_1S, 441);
PARTICLE(dd1, 1103);
PARTICLE(su0, 3201);
PARTICLE(su1, 3203);
PARTICLE(D0_star_0, 10421);
PARTICLE(xi_b_prime_minus, 5312);
PARTICLE(chi_c1, 20443);
PARTICLE(sd0, 3101);
PARTICLE(D_s1_plus_2460,20433);
PARTICLE(sigma_b_minus, 5112);
PARTICLE(junction, 88);
PARTICLE(W_plus, 24);
PARTICLE(Z, 23);
PARTICLE(xi_c_plus_prime, 4322);
PARTICLE(upsilon_1s, 553);
PARTICLE(B_c_plus, 541);
PARTICLE(B_c_star_plus, 543);
PARTICLE(xi_b_prime_0, 5322);
PARTICLE(NONEXISTENT, 0);

#undef PARTICLE

enum Status {
    FinalStateParticle = 1, Decayed_or_fragmented = 2, HardInteractionProduct = 3
};

enum ParticleType {
    Particle = 1, AntiParticle = -1
};

namespace detail {
template<typename T>
struct NameCollection;
} // detail

template<typename T>
class NameProvider {
public:
    typedef std::map<int, std::string> NameMap;

    static const std::string& Name(const T& value)
    {
        const NameMap::const_iterator iter = detail::NameCollection<T>::Names().find(value);
        Check(iter, (int) value);
        return iter->second;
    }

    static T Convert(int id)
    {
        const NameMap::const_iterator iter = detail::NameCollection<T>::Names().find(id);
        Check(iter, id);
        return (T) id;
    }

    static T Convert(const std::string& name)
    {
        NameMap::const_iterator iter = detail::NameCollection<T>::Names().begin();
        for(; iter != detail::NameCollection<T>::Names().end(); ++iter) {
            if(iter->second == name)
                break;
        }
        Check(iter, name);
        return (T) iter->first;
    }

private:
    NameProvider() {}

    template<typename Value>
    static void Check(const NameMap::const_iterator& iter, const Value& value)
    {
        if(iter == detail::NameCollection<T>::Names().end())
        {
            std::ostringstream ss;
            ss << "Unknown " << detail::NameCollection<T>::TypeName() << " '" << value << "'.";
            throw std::runtime_error(ss.str());
        }
    }
};

namespace detail {
template<>
struct NameCollection<Status> {
    static const char* TypeName() { return "particle status"; }
    static const NameProvider<Status>::NameMap& Names()
    {
        static NameProvider<Status>::NameMap names;
        if(!names.size())
        {
            names[FinalStateParticle] = "FinalStateParticle";
            names[Decayed_or_fragmented] = "Decayed_or_fragmented";
            names[HardInteractionProduct] = "HardInteractionProduct";
        }
        return names;
    }
};

template<>
struct NameCollection<ParticleType> {
    static const char* TypeName() { return "particle type"; }
    static const NameProvider<ParticleType>::NameMap& Names()
    {
        static NameProvider<ParticleType>::NameMap names;
        if(!names.size())
        {
            names[Particle] = "";
            names[AntiParticle] = "anti-";
        }
        return names;
    }
};
} // detail

struct PdgParticle {
    ParticleCode Code;
    ParticleType Type;
    PdgParticle() : Code(NONEXISTENT), Type(Particle) {}
    PdgParticle(int id) : Code(std::abs(id)), Type(id >= 0 ? Particle : AntiParticle) { }
    std::string Name() const { return NameProvider<ParticleType>::Name(Type) + Code.Name(); }
    bool operator<(const PdgParticle& other) const
    {
        if(Code < other.Code) return true;
        if(Code > other.Code) return false;
        return Type < other.Type;
    }
};
} // particles
