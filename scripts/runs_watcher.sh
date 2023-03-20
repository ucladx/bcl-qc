RUNS_DIR=$1
# # Watch `RUNS_DIR` for creation of new files
inotifywait -q -m -e create $RUNS_DIR | while read DIRECTORY EVENT FILE; do
    if [[ $EVENT == 'CREATE,ISDIR' ]]; then
        # Check if the created file matches the run directory name format
        # [[ $FILE =~ "\d{6}_[a-zA-Z0-9]{6}_\d{4}_[a-zA-Z0-9]{10}" ]] || continue
        parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")"; pwd -P )
        bash "$parent_path/run_watcher.sh" "$RUNS_DIR/$FILE"
    fi
done
