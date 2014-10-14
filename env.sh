#!/bin/bash
#
#  \file env.sh
#  \brief Setup analysis environment.
#  \author Konstantin Androsov (Siena University, INFN Pisa)
#  \author Maria Teresa Grippo (Siena University, INFN Pisa)
#
#  Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>,
#                 Maria Teresa Grippo <grippomariateresa@gmail.com>
#
#  This file is part of X->HH->bbTauTau.
#
#  X->HH->bbTauTau is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  X->HH->bbTauTau is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with X->HH->bbTauTau.  If not, see <http://www.gnu.org/licenses/>.

#cmsenv
eval `scramv1 runtime -sh`
export CMSSW_ROOT_PATH=$( cd $CMSSW_DATA_PATH/.. ; pwd )

export GCC_PATH=$( scram tool info gcc-ccompiler | grep GCC_CCOMPILER_BASE | sed s/GCC_CCOMPILER_BASE=// )
export GCC_INCLUDE_PATH=$( find $GCC_PATH/include/c++/ -maxdepth 1 -type d | tail -n 1 )

export BOOST_BASE=$( scram tool info boost | grep BOOST_BASE | sed s/BOOST_BASE=// )
export BOOST_INCLUDE_PATH=$BOOST_BASE/include
export ROOT_INCLUDE_PATH=$ROOTSYS/include
export HEPMC_INCLUDE_PATH=$( scram tool info hepmc | grep HEPMC_BASE | sed s/HEPMC_BASE=// )/include
