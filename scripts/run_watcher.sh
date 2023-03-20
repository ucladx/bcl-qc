RUN_DIR=$1
inotifywait -q -m -e create $RUN_DIR | while read DIRECTORY EVENT FILE; do
    if [[ $EVENT == 'CREATE' ]] && [[ $FILE == 'CopyComplete.txt' ]]; then
        parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
        cd "$parent_path"
        python3 ../bcl-qc.py
    fi
done
