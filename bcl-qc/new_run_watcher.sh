RUNS_DIR=$1

while true; do
# Watch `RUNS_DIR` for creation of new files
inotifywait -q -m -e create $RUNS_DIR | while read DIRECTORY EVENT FILE; do
    if [ $EVENT == 'CREATE*' ]; then
        # Check if the created file matches the run directory name format
        [[ $FILE =~ "\d{6}_[a-zA-Z0-9]{6}_\d{4}_[a-zA-Z0-9]{10}" ]] || continue
        # If it does, watch it for creation of "CopyComplete.txt"
        RUN_DIR=$FILE
        inotifywait -q -m -e create $RUN_DIR | while read DIRECTORY EVENT FILE; do
            if [ $EVENT == 'CREATE*' ] && [ $FILE == 'CopyComplete.txt' ]; then
                # Actually run the QC passes on RUN_DIR
                python3 bcl-qc.py $RUN_DIR
            fi
        done
    fi
done

done # while true
