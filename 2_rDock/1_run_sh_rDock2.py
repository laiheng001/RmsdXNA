import os, sys
import pandas as pd
from glob import glob
from multiprocessing import Pool
from spyrmsd import io, rmsd
from pymol import cmd



splitchain_folder = "/home/laiheng/RmsdXNA/1_pdb_files/split_file_filter"
docking_folder = "/home/laiheng/RmsdXNA/2_rDock/1_docking2"
os.makedirs(docking_folder, exist_ok=True)
error_log = "/home/laiheng/RmsdXNA/2_rDock/error2_docking.log"
dock_prm = "/home/laiheng/SOFTWARE/rDock-main/dock.prm"
rbcavity = "/home/laiheng/SOFTWARE/rDock-main/bin/rbcavity"
rbdock = "/home/laiheng/SOFTWARE/rDock-main/bin/rbdock"
poses_no = 100
repeat_dock = 4
rmsd_cutoff = 10
rmsdcutoff_count = 100

def parameter_text(receptor, ligand):
    text = """RBT_PARAMETER_FILE_V1.00
TITLE gart_DUD

RECEPTOR_FILE {receptor_mol2}
RECEPTOR_FLEX 3.0

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {ligand_sd}
    RADIUS 10
    SMALL_SPHERE 1.5
    MIN_VOLUME 100
    MAX_CAVITIES 99
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION        RbtCavityGridSF
    WEIGHT                  1.0
    RMAX                    0.1
    QUADRATIC               FALSE
END_SECTION
""".format(receptor_mol2 = receptor, ligand_sd = ligand)
    return text

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

def get_rmsd(sd_ligand, dock_out_sd):
    ref = io.loadmol(sd_ligand)
    mols = io.loadallmols(dock_out_sd)
    ref.strip()
    for mol in mols:
        mol.strip()
    RMSD = rmsd.symmrmsd(ref.coordinates,
                         [mol.coordinates for mol in mols],
                         ref.atomicnums,
                         mols[0].atomicnums,
                         ref.adjacency_matrix,
                         mols[0].adjacency_matrix)
    return RMSD


# for pdb_model_folder in os.listdir(splitchain_folder):
def process(pdb_model_ligand, rmsd_count):
    pdb, model, ligandno = pdb_model_ligand.split("_")
    mol2_na = os.path.join(splitchain_folder, "{pdb}_{model}/{pdb}_{model}_na.mol2".format(pdb = pdb, model = model))
    sd_ligand = os.path.join(splitchain_folder, "{pdb}_{model}/{pdb}_{model}_{ligand}.sd".format(pdb = pdb, model = model, ligand = ligandno))
    mol2_na_nearby = mol2_na.replace(".mol2", "_nearby{}.mol2".format(ligandno))
    if os.path.exists(mol2_na_nearby) == False:
        mol2_na_nearby = mol2_na
    counter = 2
    while counter <=repeat_dock:
        if os.path.exists(os.path.join(docking_folder, pdb_model_ligand + "_{}.csv".format(counter))):
            continue
        ligand_docking_folder = os.path.join(docking_folder, "{pdb}_{model}_{ligand}_{counter}".format(pdb = pdb, model = model, ligand = ligandno, counter=counter))
        os.makedirs(ligand_docking_folder,exist_ok=True)
        # prm_text = parameter_text(mol2_na, sd_ligand) # (if dont want create nearby)
        prm_text = parameter_text(mol2_na_nearby, sd_ligand)
        prm_text_filepath = os.path.join(ligand_docking_folder, pdb_model_ligand + ".prm")
        with open(prm_text_filepath, "w") as f_prm:
            f_prm.write(prm_text)
            f_prm.close()
        os.system("{} -was -d -r {}".format(rbcavity, prm_text_filepath))
        os.system("{rbdock} -i {ligand} -o {output} -r {prm}  -p {dock_prm} -n {poses_no}".format(rbdock = rbdock, ligand = sd_ligand, output = os.path.join(ligand_docking_folder, pdb_model_ligand + "_out"),
                                                                                                  prm = prm_text_filepath, dock_prm = dock_prm, poses_no = poses_no))
        try:
            dock_out_sd = os.path.join(ligand_docking_folder, pdb_model_ligand + "_out.sd")
            os.system("obabel -isd {} -omol2 -O {} -m".format(dock_out_sd, dock_out_sd.replace(".sd", '.mol2')))
            rmsd_list = get_rmsd(sd_ligand, dock_out_sd)
            sd_dict = sd_to_dict(dock_out_sd)
            df_all = pd.DataFrame.from_dict(sd_dict, orient='index')
            df_all['pose_no'] = df_all.index
            df_all['rmsd'] = rmsd_list
            df_all["Receptor"] = mol2_na
            df_all["Ligand"] = dock_out_sd
            df_all["Name"] = df_all.apply(lambda x: os.path.basename(x["Ligand"]).replace(".sd", "{}_dock{}".format(x['pose_no'],counter)), axis = 1)
            df_all = df_all[["Receptor", "Ligand", "Name", "pose_no", "rmsd", "SCORE", "SCORE.norm"]]
            df_all.reset_index(drop=True, inplace = True)
            df_all.to_csv(os.path.join(docking_folder, pdb_model_ligand + "_{}.csv".format(counter)), index=False)
            rmsd_count = rmsd_count + len(df_all[df_all["rmsd"]<=rmsd_cutoff])
            if rmsd_count >= rmsdcutoff_count:
                break
            else:
                counter += 1
        except:
            with open(error_log, "a") as flog:
                flog.write("{} unable to get rmsd\n".format(pdb_model_ligand))
                flog.close()
            
    if (counter == repeat_dock) & (rmsd_count < rmsdcutoff_count):
        with open(error_log, 'a') as f:
            f.write("{} unable to complete docking. Count = {}\n".format(pdb_model_ligand, rmsd_count))

try:
    df_finished1 = pd.read_csv('1_rDock_allresult.csv')
    df_finished2 = pd.read_csv('1_rDock_allresult2.csv')
    df_finished = pd.concat([df_finished1, df_finished2], ignore_index=True)
except FileNotFoundError:
    df_finished = pd.read_csv('1_rDock_allresult.csv')
df_finished["Lig"] = df_finished.apply(lambda x: x["Name"].split("_out")[0], axis = 1)
# Exclude those that docked for more than 400 times
count_finished = df_finished.groupby(["Lig"]).size().reset_index()
count_finished = count_finished[count_finished[0] >= repeat_dock*poses_no]
df_finished = df_finished[~df_finished["Lig"].isin(count_finished["Lig"].unique())]
# Count those with rmsd <= cutoff
df_finished = df_finished[df_finished["rmsd"] < rmsd_cutoff]
count_rmsd10_finished = df_finished.groupby(["Lig"]).size().reset_index()
count_rmsd10_finished_list = count_rmsd10_finished[count_rmsd10_finished[0] <= rmsdcutoff_count].values

pool = Pool(20)
pool.starmap(process, count_rmsd10_finished_list)

os.system("awk '(NR == 1) || (FNR > 1)' {}/*.csv > 1_rDock_allresult2.csv".format(docking_folder))