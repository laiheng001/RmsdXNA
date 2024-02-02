import os
import pandas as pd
from glob import glob
from multiprocessing import Pool
from spyrmsd import io, rmsd
import numpy as np
import signal
from pymol import cmd
import math

splitchain_folder = "/home/laiheng/RmsdXNA/1_pdb_files/split_file_filter_protonate"
docking_folder = "/home/laiheng/RmsdXNA/2_rDock/1_docking"
os.makedirs(docking_folder, exist_ok=True)
error_log = "/home/laiheng/RmsdXNA/2_rDock/error1_docking.log"
dock_prm = "/home/laiheng/SOFTWARE/rDock-main/dock.prm"
rbcavity = "/home/laiheng/SOFTWARE/rDock-main/bin/rbcavity"
rbdock = "/home/laiheng/SOFTWARE/rDock-main/bin/rbdock"
poses_no = 100

def rgyrate(selection='(all)', quiet=1):
    try:
        from itertools import izip
    except ImportError:
        izip = zip
    quiet = int(quiet)
    model = cmd.get_model(selection).atom
    x = [i.coord for i in model]
    mass = [i.get_mass() for i in model]
    xm = [(m*i,m*j,m*k) for (i,j,k),m in izip(x,mass)]
    tmass = sum(mass)
    rr = sum(mi*i+mj*j+mk*k for (i,j,k),(mi,mj,mk) in izip(x,xm))
    mm = sum((sum(i)/tmass)**2 for i in izip(*xm))
    rg = math.sqrt(rr/tmass - mm)
    if not quiet:
        print("Radius of gyration: %.2f" % (rg))
    return rg

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
    RMSD = rmsd.symmrmsd(ref.coordinates,[mol.coordinates for mol in mols],ref.atomicnums,mols[0].atomicnums,ref.adjacency_matrix,mols[0].adjacency_matrix)
    return RMSD

def get_nearby(mol2_na, sd_ligand):
    ligand_no = int(os.path.basename(sd_ligand).split("ligand")[1].replace(".sd", ""))
    cmd.reinitialize()
    cmd.load(mol2_na, object = "na")
    cmd.load(sd_ligand, object = "lig")
    box_size = max(int(3*rgyrate("lig")), 20)
    cmd.pseudoatom(pos = cmd.centerofmass("lig"), object = "lig_center")
    cmd.select("nearby_na" ,"(lig_center around {}) and na".format(box_size))
    na_nearby = mol2_na.replace(".mol2", "_nearby{}.mol2".format(ligand_no))
    cmd.save(na_nearby, "nearby_na")
    return na_nearby

def timeout_handler(signum, frame):
        raise TimeoutError
 
def process(pdb_model_ligand):
    pdb, model, ligandno = pdb_model_ligand.split("_")
    mol2_na = os.path.join(splitchain_folder, "{pdb}_{model}/{pdb}_{model}_na.mol2".format(pdb = pdb, model = model))
    sd_ligand = os.path.join(splitchain_folder, "{pdb}_{model}/{pdb}_{model}_{ligand}.sd".format(pdb = pdb, model = model, ligand = ligandno))
    ligand_docking_folder = os.path.join(docking_folder, "{pdb}_{model}_{ligand}".format(pdb = pdb, model = model, ligand = ligandno))
    os.makedirs(ligand_docking_folder, exist_ok=True)
    prm_text_filepath = os.path.join(ligand_docking_folder, pdb_model_ligand + ".prm")
    os.chdir(ligand_docking_folder)
    
    ### Perform docking ###
    def run_rdock1():
        prm_text = parameter_text(mol2_na, sd_ligand)
        with open(prm_text_filepath, "w") as f_prm:
            f_prm.write(prm_text)
            f_prm.close()
        os.system("{} -was -d -r {}".format(rbcavity, prm_text_filepath))
        os.system("{rbdock} -i {ligand} -o {output} -r {prm}  -p {dock_prm} -n {poses_no}".format(rbdock = rbdock, ligand = sd_ligand, output = os.path.join(ligand_docking_folder, pdb_model_ligand + "_out"),
                                                                                                  prm = prm_text_filepath, dock_prm = dock_prm, poses_no = poses_no))
        
    def run_rdock2():
        mol2_na_nearby = get_nearby(mol2_na, sd_ligand)
        prm_text = parameter_text(mol2_na_nearby, sd_ligand)
        with open(prm_text_filepath, "w") as f_prm:
            f_prm.write(prm_text)
            f_prm.close()
        os.system("{} -was -d -r {}".format(rbcavity, prm_text_filepath))
        os.system("{rbdock} -i {ligand} -o {output} -r {prm}  -p {dock_prm} -n {poses_no}".format(rbdock = rbdock, ligand = sd_ligand, output = os.path.join(ligand_docking_folder, pdb_model_ligand + "_out"),
                                                                                                  prm = prm_text_filepath, dock_prm = dock_prm, poses_no = poses_no))
    completed = False
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(3600)
    try:
        run_rdock1()
        completed = True
    except TimeoutError:
        signal.alarm(0)
        with open(error_log, "a") as f:
            f.write(pdb_model_ligand + " unable to dock within 1h for attempt1.\n")
            f.close()
        signal.signal(signal.SIGALRM, timeout_handler)
        signal.alarm(3600)
        try:
            run_rdock2()
            completed = True
        except TimeoutError:
            with open(error_log, "a") as f:
                f.write(pdb_model_ligand + " unable to dock within 1h for attempt2.\n")
                f.close()
        finally:
            signal.alarm(0)
    finally:
        signal.alarm(0)
        
    if completed == True:
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
            df_all["Name"] = df_all.apply(lambda x: os.path.basename(x["Ligand"]).replace(".sd", str(x["pose_no"])), axis = 1)
            df_all = df_all[["Receptor", "Ligand", "Name", "pose_no", "rmsd", "SCORE", "SCORE.norm"]]
            df_all.reset_index(drop=True, inplace = True)
            df_all.to_csv(os.path.join(docking_folder, pdb_model_ligand + ".csv"), index=False)
        except Exception as e:
            with open(error_log, "a") as flog:
                flog.write("{} unable to get rmsd, {}\n".format(pdb_model_ligand, e))
                flog.close()

# process(sys.argv[2])
input_list = list(map(lambda x: os.path.basename(x).replace(".sd", ""), glob(splitchain_folder + "/*/*ligand*.sd")))
finished = list(map(lambda x: os.path.basename(x).replace(".csv", ""), glob(docking_folder+"/*.csv")))
if os.path.exists(error_log):
    with open(error_log, "r") as f:
        problem_list = [l.strip().split()[0] for l in f.readlines() if "attempt2" in l]
        f.close()
    finished = finished + problem_list
left = list(set(input_list) - set(finished))

print(left)
# pool = Pool(min(20, len(left)))
# pool.map(process, left)

os.system("awk '(NR == 1) || (FNR > 1)' {}/*.csv > 1_rDock_allresult.csv".format(docking_folder))


