import time
import os
import sys
from bclqc import bclqc_run

def poll_directory(path):
    while True:
        if 'CopyComplete.txt' in os.listdir(path):
            bclqc_run(path.split('/')[-1], {})
            break
        else:
            time.sleep(600) # Sleep for 10 minutes

directory_path = sys.argv[-1]
if directory_path.endswith('/'):
    directory_path = directory_path[:-1]
poll_directory(directory_path)
