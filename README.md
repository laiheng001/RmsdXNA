# RmsdXNA

## Installation

> [!Note]
> The trained model can be obtained upon request

Necessary packages should be installed to run the RmsdXNA model.

Dependecies:

python >= 3.8
    pymol
    biopandas
    scipy
    numpy
    pandas
    scikit-learn
    xgboost
    rxdock

    # create a new pearsonal conda environment
    conda create -n RmsdXNA python=3.8
    conda activate RmsdXNA

    # install necessary packages
    
    conda install biopandas xgboost scipy pandas numpy scikit-learn -c conda-forge
    conda install -c conda-forge pymol-open-source
    conda install -c bioconda rxdock




## Usage

Files required: Protonated receptor .mol2 and .pdb files in the same directory, protonated screening compounds in .sdf format in the same directory.\
Code can be executed without any input to use the example receptor and ligand.

1. Create a biopandas .csv file for receptor that contains the label of the receptor atoms based on its residue and atom_name.

        python 1_process_receptor.py -receptor path/to/receptor

2. a) Perform docking of ligands onto receptor at a position using rDock and obtain .csv file containing the path to the ligand and receptor and the rDock score of the poses.

        python 2_local_dock.py -receptor path/to/receptor -folder_lig directory/to/docking/library/ -folder_dock output/directory/ -x xcenter -y ycenter -z zcenter -n_poses 100 -ncpus no_of_core

   b) If reference docking is preferred, use the command below

        python 2_ref_dock.py -receptor path/to/receptor -ref path/to/reference_ligand -folder_dock output/directory/ -n_poses 100 -ncpus no_of_core
   
   c) If there is no need for docking, use the command below

        python 2_ref_dock.py -receptor path/to/receptor -ref path/to/reference_ligand -folder_dock output/directory/ -n_poses 0 -ncpus no_of_core -prm rdock_parameter/score.prm

4. Generate features of the poses. Output file name is pose_feature.csv.

        python 3_compile_feature.py -receptor path/to/receptor -folder_dock directory/to/docked/poses/

5. Obtain RmsdXNA score of the poses. Output file name is pose_score.csv in docked_poses folder

        python 4_getscore.py -folder_dock directory/to/docked/poses/

6. Finally, select the best pose from each ligand based on the RmsdXNA score, and compile them in a .csv file

        for f in directory/to/docked/poses/*_score.csv; do sort -t "," -k 2 -n $f | head -n 2 | tail -n 1 >> compile_score.csv; done
