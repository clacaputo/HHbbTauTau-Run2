#!/bin/bash
cmsenv
export CMSSW_ROOT_PATH=$( cd $CMSSW_DATA_PATH/.. ; pwd )

export GCC_PATH=$( find $CMSSW_ROOT_PATH/external/gcc/ -maxdepth 1 -type d | tail -n 1 )
export GCC_INCLUDE_PATH=$( find $GCC_PATH/include/c++/ -maxdepth 1 -type d | tail -n 1 )

export BOOST_INCLUDE_PATH=$( find $CMSSW_ROOT_PATH/external/boost/ -maxdepth 1 -type d | tail -n 1 )/include
export ROOT_INCLUDE_PATH=$( find $CMSSW_ROOT_PATH/lcg/root/ -maxdepth 1 -type d | tail -n 1 )/include
export HEPMC_INCLUDE_PATH=$( find $CMSSW_ROOT_PATH/external/hepmc/ -maxdepth 1 -type d | tail -n 1 )/include
