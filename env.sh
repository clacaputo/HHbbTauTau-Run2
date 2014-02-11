#!/bin/bash
cmsenv
export GCC_INCLUDE_PATH=$( echo $CMSSW_RELEASE_BASE/../../../external/gcc/*/include/c++/* )
export BOOST_INCLUDE_PATH=$( echo $CMSSW_RELEASE_BASE/../../../external/boost/*/include )
export ROOT_INCLUDE_PATH=$( echo $CMSSW_RELEASE_BASE/../../../lcg/root/*/include )
export HEPMC_INCLUDE_PATH=$( echo $CMSSW_RELEASE_BASE/../../../external/hepmc/*/include )
