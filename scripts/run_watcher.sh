RUN_DIR=$1
inotifywait -q -e create $RUN_DIR | while read DIRECTORY EVENT FILE; do
    echo $EVENT
    echo $FILE
    if [[ $EVENT == 'CREATE' ]] && [[ $FILE == 'CopyComplete.txt' ]]; then
        parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
        python3 "$parent_path/../bcl-qc.py" $RUN_DIR
    fi
done
