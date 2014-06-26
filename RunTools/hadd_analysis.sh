#!/bin/bash
#
#  \file hadd_analysis.sh
#  \brief Apply hadd command to all analysis output ROOT files and place the results into the output directory.
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

if [ $# -ne 2 ] ; then
    echo "Usage: analysis_path output_dir"
    exit
fi

ANALYSIS_PATH=$1
cd $ANALYSIS_PATH

FOLDER_PATH="anaMuTau anaETau anaTauTau"
OUTPUT_DIR=$2

if [ -d "$OUTPUT_DIR" ] ; then
	echo "ERROR: output dir '$OUTPUT_DIR' already exists."
	exit
fi

mkdir -p $OUTPUT_DIR

for FOLDER in $FOLDER_PATH ; do
    SUB_FOLDERS=$( find $FOLDER -maxdepth 1 -type d -printf "%f\n" )
    for SUB_FOLDER in $SUB_FOLDERS ; do
		if [ $SUB_FOLDER = $FOLDER ] ; then
			continue
		fi
        mkdir -p $OUTPUT_DIR/$FOLDER
        if [ $SUB_FOLDER = "Radion" ] ; then
            cp $FOLDER/$SUB_FOLDER/*.root $OUTPUT_DIR/$FOLDER/
        else
            hadd $OUTPUT_DIR/$FOLDER/${SUB_FOLDER}.root $FOLDER/$SUB_FOLDER/*.root
        fi
	done
done

tar cjvf ${OUTPUT_DIR}.tar.bz2 $OUTPUT_DIR
