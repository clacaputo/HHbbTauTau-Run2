## patTrigger.py - configuration file that defines parameters related to PAT Trigger objects.

import FWCore.ParameterSet.Config as cms

# Trigger and Trigger matching
from PhysicsTools.PatAlgos.tools.trigTools import *

def applyTriggerParameters(process):
    switchOnTrigger(process)

    return
