#!/bin/bash
#
#  \file generate_file_list.sh
#  \brief Generate a set of text files that contain list of all the files from a chosen dataset.
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

if [ $# -le 2 -o $# -ge 5 ] ; then
    echo "Usage: input_path output_file_name n_jobs [reference_path]"
    exit
fi

INPUT_PATH=$1
mkdir -p $(dirname $2)
OUT_DIR=$(echo $(cd $(dirname $2); pwd))
OUT_FILE_NAME=$(basename $2)
OUT_FILE=$OUT_DIR/$OUT_FILE_NAME
TMP_DIR=$(mktemp -d)
TMP_OUT_FILE=$TMP_DIR/$OUT_FILE_NAME
N_JOBS=$3
REFERENCE_PATH=$4
PREFIX=""

if [ $# -ge 3 ] ; then
    cd "$REFERENCE_PATH"
    PREFIX="/"
    if [ "${INPUT_PATH:0:1}" == "/" ] ; then
        INPUT_PATH=${INPUT_PATH:1}
    fi
fi

find  $INPUT_PATH -name "*.root" -printf "${PREFIX}%p\n" > $TMP_OUT_FILE
cd $TMP_DIR
N_FILES=$( cat $TMP_OUT_FILE | wc -l )
N_SPLIT=$(( (N_FILES + N_JOBS - 1) / N_JOBS ))
echo "Total number of files: $N_FILES"
echo "Total number of jobs: $N_JOBS"
echo "Number of files per job: $N_SPLIT"

if [ $N_JOBS -eq 1 ] ; then
	mv "$TMP_OUT_FILE" "${OUT_FILE_NAME}.txt"
else
	split -a 3 -l $N_SPLIT $TMP_OUT_FILE
	find . -maxdepth 1 -type f -name "x*" -printf "%f ${OUT_FILE_NAME}_%f.txt\n" | xargs -n 2 mv
fi

mv *.txt $OUT_DIR
rm -rf $TMP_DIR
