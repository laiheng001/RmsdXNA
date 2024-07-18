
import pandas as pd
from pymol import cmd
from biopandas.mol2 import PandasMol2
from biopandas.pdb import PandasPdb
import json, argparse

def load_dict(json_file, json_label, alt_dict = {}):
    try:
        with open(json_file, 'r') as f:
            output_dictionary = json.load(f)  
            f.close()
    except json.decoder.JSONDecodeError:
        print("Error parsing dictionary in {} json file.".format(json_label))
        output_dictionary=alt_dict
    except FileNotFoundError:
        print("{} json file, {} not found.".format(json_label, json_file))
        output_dictionary=alt_dict
    return output_dictionary

def get_receptor_column(residue_name, atom_name, atom_type, element_symbol):
    if atom_name in metal_list:
        return 'Metal,Metal'
    else:
        Residue_name = res_dict[residue_name]
        # convert atomname
        if atom_name in atomname_convert.keys():
            atom_name = atomname_convert[atom_name]
        try:
            atom_name = res_atom_dict[residue_name][atom_name]
        except KeyError:
            atom_name = atom_name
        # Check if in standard atomname
        if atom_name in atomname["N"].keys():
            return "N," + atomname["N"][atom_name]
        elif atom_name in atomname[Residue_name]:
            return Residue_name + "," + atomname[Residue_name][atom_name]
        elif (atom_type == "O.co2") or (atom_type == "P.3"):
            return  "N," + atom_type
        else:
            return "OTH," + element_symbol
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a biopandas .csv file for receptor that contains the label of the receptor atoms based on its residue and atom_name.\
                                                  Usage: python 1_process_receptor.py -receptor path/to/receptor")
    parser.add_argument("-receptor", type=str, default="example/example_6x5n", help = "Input receptor filename ending with .pdb or .mol2 or without extension. The filename will also be used as the output .csv file.")
    parser.add_argument("-an", type=str, default='conversion_files/dict_atomname.json', help = "atom_name of residues")
    parser.add_argument("-crn", type=str, default="conversion_files/dict_residuename_convert.csv", help = "Converted residue_name")
    parser.add_argument("-cran", type=str, default="conversion_files/dict_residueatom_convert.json", help = "Converted residue atom_name")
    parser.add_argument("-can", type=str, default="conversion_files/dict_atomname_convert.json", help = "Converted atom_name")
    parser.add_argument("-m", type=str, default="conversion_files/list_metal.json", help = "Metal atom_name list")
    args = parser.parse_args()

    atomname = load_dict(args.an, 'atomname', alt_dict = {})
    res_atom_dict = load_dict(args.cran, 'residueconvert', alt_dict = {})
    atomname_convert = load_dict(args.can, 'atomnameconvert', alt_dict = {})

    try:
        df_res=pd.read_csv(args.crn)
        res_dict = {k:v for k,v in df_res.values}
    except Exception as e:
        print("Error reading residue generalization file. {}".format(e))
        print("Accept only canonical residues residue_generalization")
        res_dict={'A':'A', 'U':'U', 'T':'U', 'G':'G', 'C':'C', 'N':'N'}

    try:
        with open(args.m, "r") as f:
            metal_list= json.load(f)
            f.close()
    except Exception as e:
        print("Error reading metal_list file. {}".format(e))
        print("Setting default metal list")
        metal_list = ['Ag', 'Au', 'Ba', 'Ca', 'Cd', 'Co.oh', 'Cr.oh', 'Cs', 'Cu', 'Fe', 'K', 'Hg', 'Mg', 'Mn', 'Na', 'Ni', 'Pt', 'Rb', 'Rh', 'Ru', 'Sr', 'Tl', 'Zn']
        
    receptor_basename = args.receptor.split(".")[0]
    
    # Save receptor file as .pdb and .mol2 file
    cmd.reinitialize()
    cmd.load(args.receptor)
    cmd.save(receptor_basename + ".pdb")
    cmd.save(receptor_basename + ".mol2")
    
    na_pdb =  PandasPdb().read_pdb(receptor_basename + ".pdb")
    na_mol2 = PandasMol2().read_mol2(receptor_basename + ".mol2")

    df_na_mol2 = na_mol2.df[["atom_id", "x","y","z","atom_type"]]
    df_na_mol2.columns = ["atom_id","x_coord","y_coord","z_coord","atom_type"]
    df_na_mol2 = df_na_mol2[~df_na_mol2['atom_type'].isin(['H', 'D'])]
    df_na_pdb = pd.concat([na_pdb.df["ATOM"], na_pdb.df["HETATM"]])[["atom_name","residue_name","x_coord","y_coord","z_coord","element_symbol"]]
    df_na_join = pd.merge(df_na_pdb, df_na_mol2, on = ["x_coord","y_coord","z_coord"])
    df_na_join = df_na_join.dropna().reset_index(drop=True)
    keep_columns = df_na_join.columns
    df_na_join['Label'] = df_na_join.apply(lambda x: get_receptor_column(x["residue_name"], x["atom_name"], x["atom_type"], x["element_symbol"]), axis = 1)
    df_na_join.to_csv(receptor_basename + ".csv", index=False)
    print("Receptor label csv saved in {}".format(receptor_basename + ".csv"))
