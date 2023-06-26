## Setting up

### Install Mamba
```shell
curl -L https://github.com/conda-forge/miniforge/releases/download/22.9.0-3/Mambaforge-Linux-x86_64.sh -o mambaforge.sh
sh mambaforge.sh -bfp $HOME/mambaforge && rm -f mambaforge.sh
$HOME/mambaforge/condabin/mamba init
```

### Setup environment
```shell
mamba env create -f config/conda_env.yml
conda activate bclqc
pip3 install -r config/py_requirements.txt
```
