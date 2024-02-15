# RmsdXNA

# Installation

Necessary packages should be installed to run the OnionNet model.

Dependecies:

    python >= 3.6
    numpy  
    scipy  
    pandas 
    scikit-learn
    xgboost 
    biopandas
    pymol

To install necessary environment, create a new env with conda commands
   
    # download the package and then enter the folder
    # Git Large File System usage: https://www.atlassian.com/git/tutorials/git-lfs   
    git lfs clone https://github.com/zhenglz/onionnet.git
    cd onionnet

    # create a new pearsonal conda environment
    conda create -n onionnet python=3.6
    conda activate onionnet

    # install necessary packages
    conda install -c anaconda scipy numpy pandas
    conda install tensorflow
    conda install -c omnia mdtraj
    conda install -c openbabel openbabel
    
    # do some tests now
    python generate_features.py -h
    python predict.py -h

Or alternatively, install the packages through the environment file.

    # create a new conda environment (name: onionnet)
    conda env create -f onet_env.yaml
    conda activate onionnet

# Usage
