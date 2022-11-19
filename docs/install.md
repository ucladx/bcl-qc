## Setting up

### Install mamba
```shell=
curl -L https://github.com/conda-forge/miniforge/releases/download/22.9.0-1/Mambaforge-Linux-x86_64.sh -o mambaforge.sh
sh mambaforge.sh -bfp $HOME/mambaforge && rm -f mambaforge.sh
```

### setup environment
```shell=
ENV_NAME='bclqc'
conda init
conda create --name $ENV_NAME --file ~/config/bcl-qc/bclqc_conda_env.txt
conda activate $ENV_NAME
pip3 install -r ~/bcl-qc/config/requirements.txt
```

# start run watcher
sh ~/bcl-qc/new_run_watcher.sh <runs-directory>
