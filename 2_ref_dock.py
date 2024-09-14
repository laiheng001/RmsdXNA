import os
import pandas as pd
from glob import glob
from multiprocessing import Pool
import argparse
from pymol import cmd

def parameter_text(receptor, ligand, prm_input_path, prm_output_path):
    with open(prm_input_path, 'r') as f:
        prm_text = f.read()
        f.close()
    prm_text = prm_text.replace("{ligand_sd}",ligand)
    prm_text = prm_text.replace("{receptor_mol2}",str(receptor))
    with open(prm_output_path, "w") as f_prm:
        f_prm.write(prm_text)
        f_prm.close()

def sd_to_dict(sd_file):
    with open(sd_file, "r") as f:
        all_poses = f.read().split("$$$$")
    full_dict = {}
    for num, pose in enumerate(all_poses):
        data = list(map(lambda l: l.strip(), pose.split("\n")))
        try:
            data_score = data[data.index(">  <SCORE>"):]
            values_dict = {}
            for i in range(0, len(data_score), 3):
                try:
                    current = data_score[i:i+3]
                    header, value = current[0][4:-1], float(current[1])
                    values_dict[header] = value
                except Exception:
                    break
            full_dict[num+1] = values_dict
        except ValueError:
            break
    return full_dict
   
def process(ligand_sdf):
    if (ligand_sdf.endswith(".sd") == False) and (ligand_sdf.endswith(".sdf") == False):
        print("convert ligand file to .sdf format".format(ligand_sdf))
        return
    # os.chdir(pwd)
    ligand_name = os.path.basename(ligand_sdf).split(".")[0]
    ligand_docking_folder = os.path.join(args.folder_dock, ligand_name)
    os.makedirs(ligand_docking_folder, exist_ok=True)
    # os.chdir(ligand_docking_folder)
    os.system("rbdock -i {ligand} -o {output} -r {prm}  -p {dock_prm} -n {poses_no}".format(ligand = ligand_sdf, output = os.path.join(ligand_docking_folder, ligand_name + "_out"),
                                                                                              prm = prm_text_filepath, dock_prm = args.prm, poses_no = args.n_poses))
    # rbdock = "path/to/rbdock"
    # os.system("{rbdock} -i {ligand} -o {output} -r {prm}  -p {dock_prm} -n {poses_no}".format(rbdock = rbdock, ligand = ligand_sdf, output = os.path.join(ligand_docking_folder, ligand_name + "_out"),
    #                                                                                           prm = prm_text_filepath, dock_prm = args.prm, poses_no = poses_no))
    dock_out_sd = os.path.join(ligand_docking_folder, ligand_name + "_out.sd")
    os.system("obabel -isd {} -omol2 -O {} -m".format(dock_out_sd, dock_out_sd.replace(".sd", '.mol2')))

    sd_dict = sd_to_dict(dock_out_sd)
    df_all = pd.DataFrame.from_dict(sd_dict, orient='index')
    score_columns = list(df_all.columns)
    df_all['pose_no'] = df_all.index
    df_all["Receptor"] = args.receptor
    df_all["Ligand"] = dock_out_sd
    df_all["Name"] = df_all.apply(lambda x: "{}_{}".format(ligand_name, x["pose_no"]), axis = 1)
    df_all = df_all[["Receptor", "Ligand", "Name", "pose_no"] + list(score_columns)]
    df_all.reset_index(drop=True, inplace = True)
    df_all.to_csv(os.path.join(args.folder_dock, ligand_name + ".csv"), index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform docking of ligands onto receptor at a position using rDock and obtain .csv file containing the path to the ligand and receptor and the rDock score of the poses\
                                                  Usage: python 2_local_dock.py ")
    parser.add_argument("-receptor", type=str, default="example/example_6x5n.mol2", help = "Input protonated receptor .mol2 receptor filename.")
    parser.add_argument("-folder_lig", type=str, default="example/ligand", help = "Folder containing ligand .sdf files for docking")
    parser.add_argument("-folder_dock", type=str, default="example/docking", help = "Output folder for docked ligand files")
    parser.add_argument("-ref", type=str, help = "reference ligand")
    parser.add_argument("-cav", type=str, default="rdock_parameter/cavity_reference_ligand.txt", help = "cavity definition file")
    parser.add_argument("-prm", type=str, default="rdock_parameter/dock.prm", help = "dock prm file")
    parser.add_argument("-n_poses", type=int, default=100, help= "no. poses generated")
    parser.add_argument("-ncpus", type=int, default=1, help= "no. of CPU used to run jobs in parallel")
    args = parser.parse_args()
    
    pwd = os.getcwd()
    os.makedirs(args.folder_dock, exist_ok=True)
    prm_text_filepath = args.receptor.split(".")[0] + ".prm"
    center = (args.x,args.y,args.z)
    receptor_mol2 = args.receptor.split(".")[0] + ".mol2"
    
    ligand_sd = args.ref.split(".")[0] + ".sdf"
    cmd.load(args.ref)
    cmd.save(ligand_sd)
    
    parameter_text(receptor_mol2, ligand_sd, args.cav, prm_text_filepath)
    os.system("rbcavity -was -d -r {}".format(prm_text_filepath))
    
    ligand_list = glob(f"{args.folder_lig}/*.sdf")
    pool = Pool(args.ncpus)
    pool.map(process, ligand_list)
