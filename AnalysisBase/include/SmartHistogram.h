/*!
 * \file SmartHistogram.h
 * \brief Definition of class SmartHistogram that allows to create ROOT-compatible histograms.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-02-20 created
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

#include <deque>
#include <string>
#include <limits>
#include <stdexcept>
#include <memory>
#include <fstream>

#include <TObject.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

namespace root_ext {

class AbstractHistogram {
public:
    AbstractHistogram(const std::string& _name)
        : name(_name) {}

    virtual ~AbstractHistogram() {}

    virtual void WriteRootObject() = 0;
    virtual void DetachFromFile() {}

    const std::string& Name() const { return name; }

private:
    std::string name;
};

namespace detail {

template<typename ValueType>
class Base1DHistogram : public AbstractHistogram {
public:
    typedef typename std::deque<ValueType>::const_iterator const_iterator;

    Base1DHistogram(const std::string& name) : AbstractHistogram(name) {}

    const std::deque<ValueType>& Data() const { return data; }
    const size_t size() const { return data.size(); }
    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }

    void Fill(const ValueType& value)
    {
        data.push_back(value);
    }

    virtual void WriteRootObject()
    {
        std::unique_ptr<TTree> rootTree(new TTree(Name().c_str(), Name().c_str()));
        ValueType branch_value;
        rootTree->Branch("values", &branch_value);
        for(const ValueType& value : data) {
            branch_value = value;
            rootTree->Fill();
        }
        rootTree->Write("", TObject::kWriteDelete);
    }

private:
    std::deque<ValueType> data;
};

template<typename NumberType>
class Base2DHistogram : public AbstractHistogram {
public:
    struct Value {
        NumberType x, y;
        Value() {}
        Value(NumberType _x, NumberType _y) : x(_x), y(_y) {}
    };

    typedef typename std::deque<Value>::const_iterator const_iterator;

    Base2DHistogram(const std::string& name) : AbstractHistogram(name) {}

    const std::deque<Value>& Data() const { return data; }
    const size_t size() const { return data.size(); }
    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }

    void Fill(const NumberType& x, const NumberType& y)
    {
        data.push_back(Value(x, y));
    }

    virtual void WriteRootObject()
    {
        std::unique_ptr<TTree> rootTree(new TTree(Name().c_str(), Name().c_str()));
        NumberType branch_value_x, branch_value_y;
        rootTree->Branch("x", &branch_value_x);
        rootTree->Branch("y", &branch_value_y);
        for(const Value& value : data) {
            branch_value_x = value.x;
            branch_value_y = value.y;
            rootTree->Fill();
        }
        rootTree->Write("", TObject::kWriteDelete);
    }

private:
    std::deque<Value> data;
};

} // namespace detail

template<typename ValueType>
class SmartHistogram;

template<>
class SmartHistogram<double> : public detail::Base1DHistogram<double> {
public:
    SmartHistogram(const std::string& name) : Base1DHistogram<double>(name) {}
};

template<>
class SmartHistogram<float> : public detail::Base1DHistogram<float> {
public:
    SmartHistogram(const std::string& name) : Base1DHistogram<float>(name) {}
};

template<>
class SmartHistogram<int> : public detail::Base1DHistogram<int> {
public:
    SmartHistogram(const std::string& name) : Base1DHistogram<int>(name) {}
};

template<>
class SmartHistogram<unsigned> : public detail::Base1DHistogram<unsigned> {
public:
    SmartHistogram(const std::string& name) : Base1DHistogram<unsigned>(name) {}
};

template<>
class SmartHistogram<bool> : public detail::Base1DHistogram<bool> {
public:
    SmartHistogram(const std::string& name) : Base1DHistogram<bool>(name) {}
};

template<>
class SmartHistogram< detail::Base2DHistogram<double>::Value > : public detail::Base2DHistogram<double> {
public:
    SmartHistogram(const std::string& name) : Base2DHistogram<double>(name) {}
};

template<>
class SmartHistogram< detail::Base2DHistogram<float>::Value > : public detail::Base2DHistogram<float> {
public:
    SmartHistogram(const std::string& name) : Base2DHistogram<float>(name) {}
};

template<>
class SmartHistogram< detail::Base2DHistogram<int>::Value > : public detail::Base2DHistogram<int> {
public:
    SmartHistogram(const std::string& name) : Base2DHistogram<int>(name) {}
};

template<>
class SmartHistogram< detail::Base2DHistogram<bool>::Value > : public detail::Base2DHistogram<bool> {
public:
    SmartHistogram(const std::string& name) : Base2DHistogram<bool>(name) {}
};

template<>
class SmartHistogram<TH1D> : public TH1D, public AbstractHistogram {
public:
    SmartHistogram(const std::string& name, size_t nbins, double low, double high)
        : TH1D(name.c_str(), name.c_str(), nbins, low, high), AbstractHistogram(name), use_log_y(false),
          max_y_sf(1) {}

    SmartHistogram(const std::string& name, const std::vector<double>& bins)
        : TH1D(name.c_str(), name.c_str(), static_cast<int>(bins.size()) - 1, bins.data()), AbstractHistogram(name),
          use_log_y(false), max_y_sf(1) {}

    SmartHistogram(const std::string& name, size_t nbins, double low, double high, const std::string& x_axis_title,
                   const std::string& y_axis_title, bool _use_log_y, double _max_y_sf)
        : TH1D(name.c_str(), name.c_str(), nbins, low, high), AbstractHistogram(name), use_log_y(_use_log_y),
          max_y_sf(_max_y_sf)
    {
        SetXTitle(x_axis_title.c_str());
        SetYTitle(y_axis_title.c_str());
    }

    SmartHistogram(const std::string& name, const std::vector<double>& bins, const std::string& x_axis_title,
                   const std::string& y_axis_title, bool _use_log_y, double _max_y_sf)
        : TH1D(name.c_str(), name.c_str(), static_cast<int>(bins.size()) - 1, bins.data()), AbstractHistogram(name),
          use_log_y(_use_log_y), max_y_sf(_max_y_sf)
    {
        SetXTitle(x_axis_title.c_str());
        SetYTitle(y_axis_title.c_str());
    }

    virtual void WriteRootObject()
    {
        this->Write();
    }

    virtual void DetachFromFile()
    {
        this->SetDirectory(nullptr);
    }

    bool UseLogY() const { return use_log_y; }
    double MaxYDrawScaleFactor() const { return max_y_sf; }
    std::string GetXTitle() const { return GetXaxis()->GetTitle(); }
    std::string GetYTitle() const { return GetYaxis()->GetTitle(); }

private:
    bool use_log_y;
    double max_y_sf;
};

template<>
class SmartHistogram<TH2D> : public TH2D, public AbstractHistogram {
public:
    SmartHistogram(const std::string& name,
                   size_t nbinsx, double xlow, double xup,
                   size_t nbinsy, double ylow, double yup)
        : TH2D(name.c_str(), name.c_str(), nbinsx, xlow, xup, nbinsy, ylow, yup), AbstractHistogram(name) {}

    virtual void WriteRootObject()
    {
        this->Write();
    }

    virtual void DetachFromFile()
    {
        this->SetDirectory(nullptr);
    }
};

template<typename ValueType>
struct HistogramFactory {
    template<typename ...Args>
    static SmartHistogram<ValueType>* Make(const std::string& name, Args... args)
    {
        return new SmartHistogram<ValueType>(name, args...);
    }
};

template<>
struct HistogramFactory<TH1D> {
private:
    struct HistogramParameters {
        size_t nbins;
        double low;
        double high;
    };

    typedef std::map<std::string, HistogramParameters> HistogramParametersMap;

    static const HistogramParameters& GetParameters(const std::string& name)
    {
        static const std::string configName = "Analysis/config/histograms.cfg";
        static bool parametersLoaded = false;
        static HistogramParametersMap parameters;
        if(!parametersLoaded) {
            std::ifstream cfg(configName);
            while (cfg.good()) {
                std::string cfgLine;
                std::getline(cfg,cfgLine);
                if (!cfgLine.size() || cfgLine.at(0) == '#') continue;
                std::istringstream ss(cfgLine);
                std::string param_name;
                HistogramParameters param;
                ss >> param_name;
                ss >> param.nbins;
                ss >> param.low;
                ss >> param.high;
                if(parameters.count(param_name)) {
                    std::ostringstream ss_error;
                    ss_error << "Redefinition of default parameters for histogram '" << param_name << "'.";
                    throw std::runtime_error(ss_error.str());
                }
                parameters[param_name] = param;
              }
            parametersLoaded = true;
        }
        std::string best_name = name;
        for(size_t pos = name.find_last_of('_'); !parameters.count(best_name) && pos != 0 && pos != std::string::npos;
            pos = name.find_last_of('_', pos - 1))
            best_name = name.substr(0, pos);

        if(!parameters.count(best_name)) {
            std::ostringstream ss_error;
            ss_error << "Not found default parameters for histogram '" << name;
            if(best_name != name)
                ss_error << "' or '" << best_name;
            ss_error << "'. Please, define it in '" << configName << "'.";
            throw std::runtime_error(ss_error.str());
        }
        return parameters.at(best_name);
    }

public:
    template<typename ...Args>
    static SmartHistogram<TH1D>* Make(const std::string& name, Args... args)
    {
        return new SmartHistogram<TH1D>(name, args...);
    }

    static SmartHistogram<TH1D>* Make(const std::string& name)
    {
        const HistogramParameters& params = GetParameters(name);
        return new SmartHistogram<TH1D>(name, params.nbins, params.low, params.high);
    }

};


} // root_ext
