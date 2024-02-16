# RmsdXNA

# Installation

Necessary software to be installed in your syste,:
    Openbable: https://ccsb.scripps.edu/adfr/downloads/
    

Necessary packages should be installed to run the OnionNet model.

Dependecies:

    python >= 3.7
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
    conda create -n RmsdXNA python=3.7
    conda activate RmsdXNA

    # install necessary packages
    conda install -c bioconda rdock
    conda install openbabel

    conda install -c conda-forge -c schrodinger pymol-bundle
    conda install conda-forge::biopanda
    conda install scipy pandas numpy
    conda install conda-forge::xgboost
    
    
    # do some tests now
    python generate_features.py -h
    python predict.py -h

Or alternatively, install the packages through the environment file.

    # create a new conda environment (name: onionnet)
    conda env create -f onet_env.yaml
    conda activate onionnet

# Usage
