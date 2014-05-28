#!/bin/bash

if [ $# -lt 6 ] ; then
    echo "Usage: queue storage max_n_parallel_jobs file_list_path output_path analyzer_name args_for_make"
    exit
fi

QUEUE=$1
STORAGE=$2
MAX_N_PARALLEL_JOBS=$3
FILE_LIST_PATH=$4
OUTPUT_PATH=$5
ANALYZER_NAME=$6

WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
RUN_SCRIPT_PATH=$WORKING_PATH/RunTools/runAnalysis.sh
MAKE_PATH=$WORKING_PATH/RunTools/make_withFactory.sh

if [ $STORAGE = "Pisa" ] ; then
    PREFIX=/gpfs/ddn/cms/user/androsov
#    PREFIX="/gpfs/ddn/srm/cms"
elif [ $STORAGE = "Bari" ] ; then
    PREFIX="/lustre/cms"
elif [ $STORAGE = "Local" ] ; then
    PREFIX=$CMS_STORE
else
    echo "ERROR: unknown storage"
    exit
fi

if [ ! -d "$WORKING_PATH" ] ; then
	echo "ERROR: working path '$WORKING_PATH' does not exist."
	exit
fi

if [ ! -d "$WORKING_PATH/$FILE_LIST_PATH" ] ; then
	echo "ERROR: file list path '$WORKING_PATH/$FILE_LIST_PATH' does not exist."
	exit
fi

if [ ! -d "$OUTPUT_PATH" ] ; then
    echo "ERROR: output path '$OUTPUT_PATH' does not exist."
	exit
fi
OUTPUT_PATH=$( cd "$OUTPUT_PATH" ; pwd )

if [ ! -f "$RUN_SCRIPT_PATH" ] ; then
	echo "ERROR: script '$RUN_SCRIPT_PATH' does not exist."
	exit
fi

if [ ! -f "$MAKE_PATH" ] ; then
        echo "ERROR: script '$MAKE_PATH' does not exist."
        exit
fi

JOBS=$( find $WORKING_PATH/$FILE_LIST_PATH -maxdepth 1 -name "*.txt" -print0 | xargs -0 -n 1 basename | sed "s/\.txt//" )
#JOBS=$( find $WORKING_PATH/$FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" )

if [ "x$JOBS" = "x" ] ; then
	echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
	exit
fi

N_JOBS=$( echo "$JOBS" | wc -l )
echo "Following jobs will be submited:" $JOBS
echo "Total number of jobs to submit: $N_JOBS"

read -p "Compile these jobs and then submit (yes/no)? " -r REPLAY
if [ "$REPLAY" != "y" -a "$REPLAY" != "yes" -a "$REPLAY" != "Y" ] ; then
    echo "No jobs have been compiled or submitted."
    exit
fi

$MAKE_PATH $OUTPUT_PATH $ANALYZER_NAME
echo "Executable file is compiled."

i=0
n=0

if [ "$QUEUE" = "cms" -a "$STORAGE" = "Pisa" ] ; then
    for NAME in $JOBS ; do
        bsub -q $QUEUE -E /usr/local/lsf/work/infn-pisa/scripts/testq-preexec-cms.bash \
             -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH $OUTPUT_PATH/$ANALYZER_NAME "yes" "${@:7}"
    done
    echo "$N_JOBS have been submited in local in Pisa"
elif [ "$QUEUE" = "local" -a "$STORAGE" = "Bari" ] ; then
    for NAME in $JOBS ; do
        echo "$RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH  $OUTPUT_PATH/$ANALYZER_NAME yes ${@:7}" | \
            qsub -q $QUEUE -N $NAME -o $OUTPUT_PATH -e $OUTPUT_PATH -
    done
    echo "$N_JOBS have been submited in local in Bari"
elif [ "$QUEUE" = "fai5" ] ; then
    for NAME in $JOBS ; do
        bsub -Is -q $QUEUE -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH \
             $OUTPUT_PATH/$ANALYZER_NAME "yes" "${@:7}" &
        i=$(($i + 1))
		n=$(($n + 1))
		echo "job $n started"
        if [[ $i == $MAX_N_PARALLEL_JOBS ]] ; then
                wait
                i=0
        fi
    done
    wait
    echo "$N_JOBS finished on fai5"
elif [ "$QUEUE" = "fai" ] ; then
    for NAME in $JOBS ; do
        bsub -Is -q $QUEUE -R "select[defined(fai)]" -J $NAME \
                                                        $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH \
                                                        $OUTPUT_PATH/$ANALYZER_NAME "yes" "${@:7}" &
        i=$(($i + 1))
		n=$(($n + 1))
		echo "job $n started"
        if [[ $i == $MAX_N_PARALLEL_JOBS ]] ; then
                wait
                i=0
        fi
    done
    wait
    echo "$N_JOBS finished on fai"
elif [ $STORAGE = "Local" ] ; then
    for NAME in $JOBS ; do
        $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH $OUTPUT_PATH/$NAME "dont_set_cmsenv" \
                         $OUTPUT_PATH/$ANALYZER_NAME "yes" "${@:7}" &
        i=$(($i + 1))
        n=$(($n + 1))
        echo "job $n started"
        if [[ $i == $MAX_N_PARALLEL_JOBS ]] ; then
                wait
                i=0
        fi
    done
    wait
    echo "$N_JOBS finished on local"
else
    echo "unknow queue"
fi
