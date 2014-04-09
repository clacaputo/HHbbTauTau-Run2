#!/bin/bash

if [ $# -lt 6 ] ; then
    echo "Usage: queue storage max_n_parallel_jobs file_list_path output_path analyzer_name args_for_make"
    exit
fi

QUEUE=$1
STORAGE=$2
MAX_N_PARALLEL_JOBS=$3
FILE_LIST_PATH=$4
OUTPUT_PATH=$( cd "$5" ; pwd )
ANALYZER_NAME=$6

WORKING_PATH=$CMSSW_BASE/src/HHbbTauTau
RUN_SCRIPT_PATH=$WORKING_PATH/RunTools/runAnalysis.sh
MAKE_PATH=$WORKING_PATH/RunTools/make.sh

if [ $STORAGE = "Pisa" ] ; then
    PREFIX="/gpfs/ddn/srm/cms"
elif [ $STORAGE = "Bari" ] ; then
    PREFIX="/lustre/cms"
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

if [ ! -f "$RUN_SCRIPT_PATH" ] ; then
	echo "ERROR: script '$RUN_SCRIPT_PATH' does not exist."
	exit
fi

if [ ! -f "$MAKE_PATH" ] ; then
        echo "ERROR: script '$MAKE_PATH' does not exist."
        exit
fi

JOBS=$( find $WORKING_PATH/$FILE_LIST_PATH -maxdepth 1 -name "*.txt" -printf "%f\n" | sed "s/\.txt//" )

if [ "x$JOBS" = "x" ] ; then
	echo "ERROR: directory '$FILE_LIST_PATH' does not contains any job description."
	exit
fi

N_JOBS=$( echo "$JOBS" | wc -l )
echo "Following jobs will be submited:" $JOBS
echo "Total number of jobs to submit: $N_JOBS"

read -p "Compile these jobs (yes/no)? " -r REPLAY
if [ "$REPLAY" != "y" -a "$REPLAY" != "yes" -a "$REPLAY" != "Y" ] ; then
    echo "No jobs have been compiled."
    exit
fi

for NAME in $JOBS ; do
	echo "Compiling $NAME ..."
    $MAKE_PATH $OUTPUT_PATH $NAME $ANALYZER_NAME $FILE_LIST_PATH/${NAME}.txt $OUTPUT_PATH/${NAME}.root \
               $PREFIX @0 "${@:7}"
    MAKE_RESULT=$?
    if [ $MAKE_RESULT -ne 0 ] ; then
        echo "ERROR: failed to make job $NAME"
        exit
    fi
done
echo "All jobs compiled."
read -p "Submit these jobs (yes/no)? " -r REPLAY
if [ "$REPLAY" != "y" -a "$REPLAY" != "yes" -a "$REPLAY" != "Y" ] ; then 
	echo "No jobs have been submited."
	exit
fi

i=0
n=0

if [ "$QUEUE" = "local" -a "$STORAGE" = "Pisa" ] ; then
    for NAME in $JOBS ; do
        bsub -q $QUEUE -E /usr/local/lsf/work/infn-pisa/scripts/testq_pre-cms.bash \
             -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH $OUTPUT_PATH/$NAME
    done
    echo "$N_JOBS have been submited in local in Pisa"
elif [ "$QUEUE" = "local" -a "$STORAGE" = "Bari" ] ; then
    for NAME in $JOBS ; do
        echo "$RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH $OUTPUT_PATH/$NAME" | \
            qsub -q $QUEUE -N $NAME -o $OUTPUT_PATH -e $OUTPUT_PATH -
    done
    echo "$N_JOBS have been submited in local in Bari"
elif [ "$QUEUE" = "fai5" ] ; then
    for NAME in $JOBS ; do
        bsub -Is -q $QUEUE -J $NAME $RUN_SCRIPT_PATH $NAME $WORKING_PATH $OUTPUT_PATH $OUTPUT_PATH/$NAME &
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
                                                        $OUTPUT_PATH/$NAME &
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
else
    echo "unknow queue"
fi
