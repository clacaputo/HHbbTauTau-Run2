/*!
 * \file LightBaseBigTreeAnalyzer.h
 * \brief Definition of LightBaseBigTreeAnalyzer class, the base class for separate studies on big trees.
 * \author Konstantin Androsov (Siena University, INFN Pisa)
 * \author Maria Teresa Grippo (Siena University, INFN Pisa)
 * \date 2014-10-24 created
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

#include "Analysis/source/HHbbetau_FlatTreeProducer.C"
#include "Analysis/source/HHbbmutau_FlatTreeProducer.C"
#include "Analysis/source/HHbbtautau_FlatTreeProducer.C"

namespace analysis {

class LightBaseBigTreeAnalyzer : public HHbbetau_FlatTreeProducer, public HHbbmutau_FlatTreeProducer,
                                 public HHbbtautau_FlatTreeProducer {
public:
    LightBaseBigTreeAnalyzer(const std::string& channelName, const std::string& inputFileName,
                             const std::string& outputFileName, const std::string& configFileName,
                             const std::string& _prefix = "none", size_t _maxNumberOfEvents = 0)
        : BaseFlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          HHbbetau_FlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          HHbbmutau_FlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents),
          HHbbtautau_FlatTreeProducer(inputFileName, outputFileName, configFileName, _prefix, _maxNumberOfEvents)
    {
        std::istringstream ss_channel(channelName);
        ss_channel >> channel;
        if(channel != Channel::ETau && channel != Channel::MuTau && channel != Channel::TauTau)
            throw exception("Unsupported channel ") << channel;
        writeFlatTree = false;
    }

    virtual analysis::BaseAnalyzerData& GetAnaData() override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::GetAnaData();
        if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::GetAnaData();
        return HHbbtautau_FlatTreeProducer::GetAnaData();
    }

    virtual void ProcessEvent(std::shared_ptr<const analysis::EventDescriptor> _event) override
    {
        BaseAnalyzer::ProcessEvent(_event);

        const SelectionResults& selection = ApplyBaselineSelection();

        AnalyzeSelection(selection);
        if(channel == Channel::ETau)
            AnalyzeETauSelection(HHbbetau_FlatTreeProducer::selection);
        else if(channel == Channel::MuTau)
            AnalyzeMuTauSelection(HHbbmutau_FlatTreeProducer::selection);
        else
            AnalyzeTauTauSelection(HHbbtautau_FlatTreeProducer::selection);
    }

protected:
    virtual void AnalyzeSelection(const SelectionResults& selection) {}
    virtual void AnalyzeETauSelection(const SelectionResults_etau& selection) {}
    virtual void AnalyzeMuTauSelection(const SelectionResults_mutau& selection) {}
    virtual void AnalyzeTauTauSelection(const SelectionResults_tautau& selection) {}

    virtual analysis::SelectionResults& ApplyBaselineSelection() override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::ApplyBaselineSelection();
        else if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::ApplyBaselineSelection();
        return HHbbtautau_FlatTreeProducer::ApplyBaselineSelection();
    }

    virtual void CalculateTriggerWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateTriggerWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateTriggerWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateTriggerWeights(higgs);
    }

    virtual void CalculateIsoWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateIsoWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateIsoWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateIsoWeights(higgs);
    }

    virtual void CalculateIdWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateIdWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateIdWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateIdWeights(higgs);
    }

    virtual void CalculateDMWeights(const analysis::Candidate& higgs) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::CalculateDMWeights(higgs);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::CalculateDMWeights(higgs);
        else
            HHbbtautau_FlatTreeProducer::CalculateDMWeights(higgs);
    }

    virtual analysis::Candidate SelectTau(size_t id, cuts::ObjectSelector* objectSelector,
                                          root_ext::AnalyzerData& _anaData,
                                          const std::string& selection_label) override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::SelectTau(id, objectSelector, _anaData, selection_label);
        else if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::SelectTau(id, objectSelector, _anaData, selection_label);
        return HHbbtautau_FlatTreeProducer::SelectTau(id, objectSelector, _anaData, selection_label);
    }

    virtual analysis::Candidate SelectSignalTau(size_t id, cuts::ObjectSelector* objectSelector,
                                                root_ext::AnalyzerData& _anaData,
                                                const std::string& selection_label) override
    {
        if(channel == Channel::ETau)
            return HHbbetau_FlatTreeProducer::SelectSignalTau(id, objectSelector, _anaData, selection_label);
        else if(channel == Channel::MuTau)
            return HHbbmutau_FlatTreeProducer::SelectSignalTau(id, objectSelector, _anaData, selection_label);
        return HHbbtautau_FlatTreeProducer::SelectSignalTau(id, objectSelector, _anaData, selection_label);
    }

    virtual void FillFlatTree(const analysis::SelectionResults& selection) override
    {
        if(channel == Channel::ETau)
            HHbbetau_FlatTreeProducer::FillFlatTree(selection);
        else if(channel == Channel::MuTau)
            HHbbmutau_FlatTreeProducer::FillFlatTree(selection);
        else
            HHbbtautau_FlatTreeProducer::FillFlatTree(selection);
    }

private:
    Channel channel;
};

} // namespace analysis
