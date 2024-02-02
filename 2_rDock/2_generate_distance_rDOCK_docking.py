import numpy as np
import os
import argparse
import pandas as pd
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
from multiprocessing import Pool
from scipy.spatial import distance
from glob import glob

"""
lig in docking folder, .mol2 format
nohup sh -c 'for f in 1*/*.sd; do fn=$(basename -s ".sd" $f); dir=$(dirname $f); cd $dir; obabel -isd "$fn".sd -omol2 -O "$fn".mol2 -m; cd ../; done' >/dev/null 2>&1 &

na in pdb_splitchain folder, needs .mol2 and .pdb.
If not script will convert 
"""

error_log = "error2_feature2.log"
if os.path.exists(error_log):
    with open(error_log, 'r') as f:
        error_list = [l.strip().split()[0] for l in f.readlines()]
        f.close()
else:
    error_list = []
    
def process(na, lig, filename):
    ### Read lig.mol2 file ###
    try:
        lig_mol2 = PandasMol2().read_mol2(lig)
        df_lig_mol2 = lig_mol2.df[["x","y","z","atom_type"]]
        df_lig_mol2_noH = df_lig_mol2[~df_lig_mol2['atom_type'].isin(['H', 'D'])] ###
        df_lig_mol2_noH.reset_index(inplace=True)
    except:
        with open(error_log, "a") as f:
            f.write("Error reading ligand file: {}\n".format(lig))
        f.close()
        return
    
    ### Read na.pdb and na.mol2 file. Obabel conversion to na.mol2 if it does not exist ###
    dir, fn = os.path.dirname(na), os.path.basename(na).split(".")[0]
    na_pdb_fn, na_mol2_fn = os.path.join(dir, fn + ".pdb"), os.path.join(dir, fn + ".mol2")
    try:
        if os.path.exists(na_pdb_fn) and os.path.exists(na_mol2_fn):
            na_pdb =  PandasPdb().read_pdb(na_pdb_fn)
            na_mol2 = PandasMol2().read_mol2(na_mol2_fn)
        elif os.path.exists(na_pdb_fn) and (os.path.exists(na_mol2_fn) == False): #missing .mol2 file
            na_pdb =  PandasPdb().read_pdb(na_pdb_fn)
            os.system("obabel -ipdb {} -omol2 -O {}".format(na_pdb_fn, na_mol2_fn))
            try:
                na_mol2 = PandasMol2().read_mol2(na_mol2_fn)
            except ValueError:
                with open(error_log, "a") as f:
                    f.write("Obabel conversion of nucleic acid fail: {}".format(na_mol2_fn))
                f.close()
                return
        else:
            with open(error_log, "a") as f:
                f.write("Missing na.pdb: {}".format(na_pdb_fn))
            f.close()
            return
    except:
        with open(error_log, "a") as f:
            f.write("Error reading receptor file: {}\n".format(fn))
        f.close()
        return
    
    try:
        ### load pymol to select those within cutoff distance
        df_na_mol2 = na_mol2.df[["x","y","z","atom_type"]] ###
        df_na_mol2.columns = ["x_coord","y_coord","z_coord","atom_type"]
        df_na_mol2 = df_na_mol2[~df_na_mol2['atom_type'].isin(['H', 'D'])]
        df_na_pdb = pd.concat([na_pdb.df["ATOM"], na_pdb.df["HETATM"]])[["atom_name","residue_name","x_coord","y_coord","z_coord"]] ###
        df_na_join = pd.merge(df_na_pdb, df_na_mol2, on = ["x_coord","y_coord","z_coord"])
        df_na_join = df_na_join.dropna().reset_index(drop=True)
        df_na_join['Label'] = df_na_join.apply(lambda x: ",".join([x["residue_name"], x["atom_name"], x["atom_type"]]), axis = 1)
        
        output_folder_lig = os.path.join(args.output_folder, filename.split("_out")[0])
        os.makedirs(output_folder_lig, exist_ok=True)
        fileout = os.path.join(output_folder_lig, filename + ".out")
        distance_matrix = np.array(distance.cdist(df_na_join[['x_coord','y_coord','z_coord']].values, df_lig_mol2_noH[['x','y','z']].values))
        indices = np.where(distance_matrix[:,:]<=args.cutoff)
        if len(indices[0]) == 0:
            with open(error_log, "a") as f:
                f.write("{} no nearby interaction\n".format(filename))
            f.close()
        else:
            with open(fileout, "w") as f:
                for x,y in zip(indices[0], indices[1]):
                    na_label = df_na_join.loc[x,'Label']
                    lig_lab = df_lig_mol2_noH.loc[y,'atom_type']
                    distance_value = distance_matrix[x,y]
                    f.write("{},{},{:.10f}\n".format(na_label, lig_lab, distance_value))
                f.close()
    except Exception as e:
        with open(error_log, "a") as f:
            f.write("{} problem processing, {}\n".format(fn, e))
        f.close()
        return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="python generate_distance_feature.py -inp rdock_allresult_all.csv -cutoff 10 -output_folder 3_distance_rdock_docking -ncpus 8")
    parser.add_argument("-inp", type=str, default="2_rDock_allresult_filter.csv")
    parser.add_argument("-cutoff", type=float, default=10)
    parser.add_argument("-output_folder", type=str, default="2_distance_rdock_docking")
    parser.add_argument("-ncpus", type=int, default=20)
    args = parser.parse_args()
    
    os.makedirs(args.output_folder, exist_ok = True)
    ### load input.dat file ###
    """
    For pdb_splitchain, na.pdb lig.pdb
    For docking, na.pdb lig.mol2
    """

    if os.path.exists(args.inp) == False:
        df1 = pd.read_csv("1_rDock_allresult.csv")
        df2 = pd.read_csv("1_rDock_allresult2.csv")
        df_all = pd.concat([df1,df2], ignore_index=True)
        df_all = df_all[~df_all["Name"].isin(error_list)]
        df_all = df_all[df_all["rmsd"] < 10]
        df_all["Lig"] = df_all.apply(lambda x: x["Name"].split("_out")[0], axis = 1)
        count_cutoff = df_all.groupby(["Lig"]).size().reset_index()
        count_cutoff = count_cutoff[count_cutoff[0] >= 100]

        df_all = df_all[df_all["Lig"].isin(list(count_cutoff["Lig"].unique()))]
        df = df_all.groupby('Lig').head(100)

        df["Ligand_pose"] = df.apply(lambda x: os.path.join(os.path.dirname(x['Ligand']), os.path.basename(x['Ligand']).replace(".sd", "{}.mol2".format(x['pose_no']))), axis=1)
        df["Lig_out"] = df.apply(lambda x: os.path.join(args.output_folder, x['Lig'], x['Name'] + ".out"), axis=1)
        df.to_csv(args.inp, index = False)
    else:
        df = pd.read_csv(args.inp)
    
    # Remove those that are finished already
    df = df[(~df["Lig_out"].isin(glob(args.output_folder + "/*/*.out")))]
    inputs = df[['Receptor', 'Ligand_pose', "Name"]].values
    # print(inputs)
    pool = Pool(args.ncpus)
    pool.starmap(process, inputs)
    # process(*inputs[0])
    
### Remove empty file ###
# for file in glob("3_distance_rdock_docking/*"):
#      if os.stat(file).st_size == 0:
#          os.remove(file)


