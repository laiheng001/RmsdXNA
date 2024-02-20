# RmsdXNA

## Installation

Necessary packages should be installed to run the RmsdXNA model.

Dependecies:

python >= 3.7
    rdock
    openbabel
    pymol
    biopandas
    scipy
    numpy
    pandas
    scikit-learn
    xgboost

    # create a new pearsonal conda environment
    conda create -n RmsdXNA python=3.7
    conda activate RmsdXNA

    # install necessary packages
    conda install -c bioconda rdock
    conda install openbabel
    conda install -c conda-forge -c schrodinger pymol-bundle
    conda install conda-forge::biopanda
    conda install scipy pandas numpy
    conda install scikit-learn
    conda install conda-forge::xgboost
    

## Usage

Files required: Protonated receptor .mol2 and .pdb files in the same directory, protonated screening compounds in .sdf format in the same directory.\
Code can be executed without any input to use the example receptor and ligand.

Create a biopandas .csv file for receptor that contains the label of the receptor atoms based on its residue and atom_name.

    python 1_process_receptor.py -receptor path/to/receptor

Perform docking of ligands onto receptor at a position using rDock and obtain .csv file containing the path to the ligand and receptor and the rDock score of the poses.

    python 2_local_dock.py -receptor path/to/receptor -folder_lig directory/to/docking/library/ -folder_dock output/directory/ -x xcenter -y ycenter -z zcenter -n_poses 100 -ncpus no_of_core

Generate features of the poses. Output file name is pose_feature.csv.

    python 3_compile_feature.py -rec_csv path/to/receptor.csv -folder_dock directory/to/docked/poses/

Obtain RmsdXNA score of the poses. Output file name is pose_score.csv.

    python 4_getscore.py -folder_dock directory/to/docked/poses/
