#!/bin/bash

if [ $# -ne 3 ] ; then
    echo "Usage: source destination directory"
    exit
fi

SOURCE_NAME=$1
DESTINATION_NAME=$2
DIRECTORY_NAME=$3

MIN_SPEED=5000000 # B/s
NUMBER_OF_STREAMS=4

if [ "$SOURCE_NAME" = "T2_IT_Pisa" ] ; then
    SOURCE_PREFIX="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN="
elif [ "$SOURCE_NAME" = "T2_IT_Bari" ] ; then
    SOURCE_PREFIX="srm://storm-se-01.ba.infn.it:8444/srm/managerv2?SFN="
else
    echo "ERROR: unknown source name '$SOURCE_NAME'."
    exit
fi

if [ "$DESTINATION_NAME" = "T2_IT_Pisa" ] ; then
    DESTINATION_PREFIX="srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN="
elif [ "$DESTINATION_NAME" = "T2_IT_Bari" ] ; then
    DESTINATION_PREFIX="srm://storm-se-01.ba.infn.it:8444/srm/managerv2?SFN="
else
    echo "ERROR: unknown destination name '$DESTINATION_NAME'."
    exit
fi

FULL_SOURCE_PATH=""
FULL_DESTINATION_PATH=""

LS_RESULT=$( lcg-ls -l $SOURCE_PREFIX/$DIRECTORY_NAME )
RESULT=$?
if [ $RESULT -ne 0 ] ; then
    echo "ERROR: directory '$DIRECTORY_NAME' not found in '$SOURCE_NAME'."
fi

echo "$LS_RESULT" | while read FILE_DESC ; do
    FILE_SIZE=$( echo $FILE_DESC | awk '{print $5}' )
    FILE_NAME=$( echo $FILE_DESC | awk '{print $7}' )
    TIMEOUT=$(( $FILE_SIZE / $MIN_SPEED ))

    WATCHDOG_KILLED=1
    while [ $WATCHDOG_KILLED -ne 0 ] ; do
        LS_RESULT_DST=$( lcg-ls --srm-timeout=$TIMEOUT --connect-timeout=$TIMEOUT -l $DESTINATION_PREFIX/$FILE_NAME )
        NEED_COPY=1
        ITERATION=0
        echo "$LS_RESULT_DST" | while read FILE_DESC_DST ; do
            if [ $ITERATION -ne 0 ] ; then
                echo "ERROR: more than one file in destination for input 'FILE_NAME'"
                exit
            fi
            ITERATION=$(($ITERATION + 1))
            FILE_SIZE_DST=$( echo $FILE_DESC_DST | awk '{print $5}' )
            if [ $FILE_SIZE -eq $FILE_SIZE_DST ] ; then
                echo "File $FILE_NAME is already transfered."
                NEED_COPY=0
                continue
            fi
            lcg-del --srm-timeout=$TIMEOUT --connect-timeout=$TIMEOUT $DESTINATION_PREFIX/$FILE_NAME
        done
        if [ $NEED_COPY -ne 1 ] ; then continue ; fi

        ( lcg-cp --srm-timeout=$TIMEOUT --connect-timeout=$TIMEOUT -b -D srmv2 -v -n $NUMBER_OF_STREAMS \
            $SOURCE_PREFIX/$FILE_NAME $DESTINATION_PREFIX/$FILE_NAME ) & LCG_PID=$!
        ( sleep $TIMEOUT && kill -9 $LCG_PID ) & WATCHDOG_PID=$!
        wait $LCG_PID
        kill -9 $WATCHDOG_PID &> /dev/null
        WATCHDOG_KILLED=$?
    done
done
