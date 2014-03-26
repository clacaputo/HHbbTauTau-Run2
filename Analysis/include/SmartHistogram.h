/*!
 * \file SmartHistogram.h
 * \brief Definition of class SmartHistogram that allows to create ROOT-compatible histograms.
 * \author Konstantin Androsov (INFN Pisa, Siena University)
 * \author Maria Teresa Grippo (INFN Pisa, Siena University)
 * \date 2014-02-20 created
 */

#pragma once

#include <deque>
#include <string>
#include <limits>
#include <stdexcept>

#include <TObject.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

namespace root_ext {

class AbstractHistogram {
public:
    AbstractHistogram(const std::string& _name)
        : name(_name) {}

    virtual ~AbstractHistogram() {}

    virtual void WriteRootObject() = 0;

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
        : TH1D(name.c_str(), name.c_str(), nbins, low, high), AbstractHistogram(name) {}

    virtual void WriteRootObject()
    {
        this->Write();
    }
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
};

template<typename ValueType>
struct HistogramFactory {
    template<typename ...Args>
    static SmartHistogram<ValueType>* Make(const std::string& name, Args... args)
    {
        return new SmartHistogram<ValueType>(name, args...);
    }
};

} // root_ext
