import numpy as np
import os
import pandas as pd
from biopandas.mol2 import PandasMol2
from multiprocessing import Pool
from scipy.spatial import distance
from glob import glob
import json
from pymol import cmd
import argparse


# for ligand_sd in ligand_list[:2]:
def process(ligand_csv):
    df = pd.read_csv(ligand_csv)
    df_all = pd.DataFrame()
    feature_error_list = []
    
    for index, row in df.iterrows():
        # try:
            ligand_mol2 = os.path.join(os.path.dirname(row['Ligand']), os.path.basename(row['Ligand']).replace(".sd", str(row["pose_no"]) + ".mol2"))
            cmd.reinitialize()
            cmd.load(ligand_mol2, "lig")
            cmd.load(row['Receptor'], "rec")
            cmd.remove("elem H")
            cmd.select("nearbylig", f"(lig around {args.cutoff}) and rec")
            receptor_nearby_id = cmd.identify("nearbylig")
            df_na_nearby = df_na_join[df_na_join["atom_id"].isin(receptor_nearby_id)]
            
            
            df_lig_mol2 = PandasMol2().read_mol2(ligand_mol2).df[["x","y","z","atom_type"]]
            df_lig_mol2_noH = df_lig_mol2[~df_lig_mol2['atom_type'].isin(['H', 'D'])] ###
            df_na_nearby.reset_index(inplace=True)
            df_lig_mol2_noH.reset_index(inplace=True)
            
            label_value_current = label_value_dict.copy()
            label_value_current["Name"] = row['Name']
            for score_label in score_columns:
                try:
                    label_value_current[score_label] = row[score_label]
                except:
                    continue
            
            distance_matrix = np.array(distance.cdist(df_na_nearby[['x_coord','y_coord','z_coord']].values, df_lig_mol2_noH[['x','y','z']].values))
            indices = np.where((distance_matrix[:,:] <= args.cutoff)) #& (distance_matrix[:,:] > 2)
            if len(indices[0]) == 0:
                with open(args.log, "a") as f:
                    f.write("{} no nearby interaction\n".format(os.path.basename(ligand_mol2).replace(".mol2","")))
                f.close()
            else:
                for x,y in zip(indices[0], indices[1]):
                    na_label = df_na_nearby.loc[x,'Label']
                    lig_label = df_lig_mol2_noH.loc[y,'atom_type']
                    try:
                        Lig_label = lig_atom_dict[lig_label]
                    except:
                        Lig_label = lig_label
                    distance_value = distance_matrix[x,y]
                    feature1, feature2 = pow(distance_value, -1), pow(distance_value, -6)
                    try:
                        label_value_current["{},{}:1/r".format(na_label, Lig_label)] += feature1
                        label_value_current["{},{}:1/r6".format(na_label, Lig_label)] += feature2
                    except:
                        feature_error_list.append("{},{}".format(na_label, Lig_label))
                df_out = pd.DataFrame([label_value_current])
                df_out = df_out[['Name'] + feature_col]
                df_all = pd.concat([df_all, df_out], ignore_index=True)
        # except Exception as e:
        #     with open(args.log, "a") as f:
        #         f.write("{} {}\n".format(ligand_mol2, e))
        #     f.close()
        #     continue
    if len(feature_error_list) > 0:
        with open(args.log, "a") as f:
            f.write("{} unseen feature: {}\n".format(os.path.basename(row['Ligand']).replace(".sd",""), set(feature_error_list)))
        f.close()
    if df_all.empty == False:
        df_all.to_csv(ligand_csv.replace(".csv", '_feature.csv'), index = False)
    # os.system("rm {}/*.mol2".format(os.path.dirname(ligand_sd)))

# process("docking_CB/core_581950/core_581950_out.sd")      
# pool = Pool(16)
# pool.map(process, ligand_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate features of the poses")
    parser.add_argument("-receptor", type=str, default="example/example_6x5n.pdb", help = "Input receptor file")
    parser.add_argument("-folder_dock", type=str, default="example/docking", help = "Folder containing ligand docking information .csv files")
    parser.add_argument("-cutoff", type=float, default=8, help = "Cutoff distance for feature extraction")
    parser.add_argument("-log", type=str, default="error_feature.log", help = ".log file to record error")
    parser.add_argument("-cla", type=str, default="conversion_files/dict_ligatom_convert.csv", help = "Converted ligand atom_name")
    parser.add_argument("-column", type=str, default="feature_column.json", help = "json file containing dataset columns")
    parser.add_argument("-ncpus", type=int, default=1, help= "no. of CPU used to run jobs in parallel")
    args = parser.parse_args()
    
    
    
    df_na_join = pd.read_csv(args.receptor.split(".")[0] + ".csv")

    try:
        df_lig=pd.read_csv(args.cla)
        lig_atom_dict = {k:v for k,v in df_lig.values}
    except Exception as e:
        print("Error reading residue generalization file. {}".format(e))
        print("Ignoring ligand atom conversion")
        lig_atom_dict={}
        
    with open(args.column) as f_dict:
        label_value_dict = json.load(f_dict)
        f_dict.close()
    feature_col = sorted(label_value_dict.keys())
    score_columns = list(filter(lambda x: x.startswith("SCORE"), label_value_dict.keys()))
    
    ligand_list = glob(args.folder_dock + "/*.csv")
    pool = Pool(args.ncpus)
    pool.map(process, ligand_list)
    
