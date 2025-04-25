#!/usr/bin/env bash

# constants
readonly PYTHON_PATH=/home/skrhakv/CryptoBench/holo_venv/bin/activate
readonly OUTPUT_PATH=/home/skrhakv/CryptoBench/data/A-filter-ahojdb/holo_only_data
# CAUTION: if you have multiple folders (i.e. I got 3 folders - storage1, storage2, storage3),
# you need to run this for each of these folders separately 
readonly INPUT_PATH=/home/skrhakv/CryptoBench/data/ahoj-db/ahojdb_v2c/data

# vars
BASEDIR=$(realpath $(dirname $0))
TMP_DIR=$BASEDIR/tmp-data

# deactivate old venv and activate venv
deactivate

source $PYTHON_PATH
mkdir $TMP_DIR

for batch_path in $INPUT_PATH/*/; do
    BATCH=${batch_path: -4:-1}

    mkdir $TMP_DIR/$BATCH-pairs

    for directory in $INPUT_PATH/$BATCH/*/; do
        cd $directory
        if ! [ -f "apo_filtered_sorted_results.csv" ]; then
            continue
        fi
        if test "$(wc -l < apo_filtered_sorted_results.csv)" -gt 1; then 
            cp -r $directory $TMP_DIR/$BATCH-pairs
        fi
    done

    rm -rf $TMP_DIR/$BATCH

    # make pairs
    cd $BASEDIR
    python3 main.py $BATCH $TMP_DIR $OUTPUT_PATH 2>> error.txt
    rm -rf $TMP_DIR/$BATCH-pairs
    break
done

