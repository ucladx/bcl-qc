RUN_DIR=$1

while true; do
        inotifywait -q -m -e create $RUN_DIR | while read DIRECTORY EVENT FILE; do
            if [[ $EVENT == 'CREATE*' ]] && [[ $FILE == 'CopyComplete.txt' ]]; then
                python3 ../bcl-qc.py -m $RUN_DIR
            fi
        done
done # while true
