#!/bin/bash
#
#  \file scan_trees.sh
#  \brief Scan root files in the given directory for the specified event.
#  \author Konstantin Androsov (Siena University, INFN Pisa)
#  \date 2014-11-25 created
#
#  Copyright 2014 Konstantin Androsov <konstantin.androsov@gmail.com>
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

if [ $# -ne 2 ] ; then
    echo "Usage: trees_path even_id"
    exit
fi

TREES_PATH=$1
EVENT_ID=$2

SCAN_FILE="RunTools/source/Scan.C"

find $TREES_PATH -type f -exec root -b -l -q $SCAN_FILE\($EVENT_ID,\"\{\}\"\) \; > /dev/null
