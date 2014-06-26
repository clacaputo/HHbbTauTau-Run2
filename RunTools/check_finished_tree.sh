#!/bin/bash
#
#  \file check_finished_tree.sh
#  \brief Check which jobs were finished with success and copy jobs that are not successful in a new directory.
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

if [ $# -ne 3 ] ; then
    echo "Usage: file_list_path job_result new_file_list_path"
    exit
fi

FILE_LIST_PATH=$1
FILE_JOB_RESULT=$2
NEW_FILE_LIST_PATH=$3

if [ ! -d "$FILE_LIST_PATH" ] ; then
        echo "ERROR: file list path '$FILE_LIST_PATH' does not exist."
        exit
fi

JOBS=$( find $FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" | sort )

if [ "x$JOBS" = "x" ] ; then
        echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
        exit
fi

SUCCESSFULL_JOBS=$( cat $FILE_JOB_RESULT | sed -n "s/\(^0 \)\([^ ]*\)\(.*\)/\2/p" | sort )

echo -e "$JOBS"\\n"$SUCCESSFULL_JOBS" | sort | uniq -u | xargs -n 1 printf "$FILE_LIST_PATH/%b.txt $NEW_FILE_LIST_PATH/\n" \
     | xargs -n 2 cp 
