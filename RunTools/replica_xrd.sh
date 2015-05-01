#!/bin/bash
#
#  \file replica_xrd.sh
#  \brief Run data replica between sites using xrdcp.
#  \author Konstantin Androsov (University of Siena, INFN Pisa)
#  \author Maria Teresa Grippo (University of Siena, INFN Pisa)
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

if [ $# -lt 2 -o $# -gt 3 ] ; then
    echo "Usage: file_list_path output_path [max_n_files]"
    exit
fi

FILE_LIST_PATH=$1
OUTPUT_PATH=$2

if [ $# -eq 3 ] ; then
    MAX_N_FILES=$3
else
    MAX_N_FILES=0
fi

SERVER="root://xrootd.ba.infn.it/"

if [ -d $FILE_LIST_PATH ] ; then
    JOB_FILES=$( find $FILE_LIST_PATH -maxdepth 1 -name "*.txt" )
elif [ -f $FILE_LIST_PATH ] ; then
    JOB_FILES=$FILE_LIST_PATH
else
    echo "ERROR: file list path '$FILE_LIST_PATH' does not exists."
    exit
fi

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' does not exist."
    exit
fi

N_JOB_FILES=$( echo "$JOB_FILES" | wc -l )
if [ $N_JOB_FILES -eq 0 ] ; then
    echo "ERROR: no job files found in '$FILE_LIST_PATH'"
    exit
fi

n=0

for JOB_FILE in $JOB_FILES ; do
    FILES_TO_TRANSFER=$( cat $JOB_FILE | sed '/^\s*$/d' )
    N_FILES_TO_TRANSFER=$( echo "$FILES_TO_TRANSFER" | sed '/^\s*$/d' | wc -l )
    if [ $N_FILES_TO_TRANSFER -eq 0 ] ; then
        echo "ERROR: input job file '$JOB_FILE' is empty."
        exit
    fi
    for FILE in $FILES_TO_TRANSFER ; do
        FILE_NAME=$( basename $FILE )
        if [ -f $OUTPUT_PATH/$FILE_NAME ] ; then
            echo "File '$FILE' is already exists. Skipping it."
            continue
        fi
        echo "Starting transfering file '$FILE'..."
        while [ 1 ] ; do
            xrdcp -force "${SERVER}${FILE}" $OUTPUT_PATH
            RESULT=$?
            if [ $RESULT -eq 0 ] ; then
                break
            fi
            echo "Transfer of '$FILE' failed. Will retry in 1 minute..."
            sleep 60
        done
        n=$(($n + 1))
        if [[ $MAX_N_FILES -ne 0 && $n -ge $MAX_N_FILES ]] ; then
            echo "Max number of files was transfered. Exiting..."
            exit
        fi
    done
done

echo "Output path size: " $( du -s $OUTPUT_PATH )
