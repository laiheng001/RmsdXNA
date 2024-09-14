import xgboost as xgb
import pandas as pd
from glob import glob
from multiprocessing import Pool
import os, json
import argparse

def process(csv_file):
    try:
        df = pd.read_csv(csv_file)
        df = df[["Name"] + feature_col]
        X_test = df.set_index(["Name"]).fillna(0)
        df["Lig"] = df.apply(lambda x: os.path.basename(x["Name"]).split("_")[0], axis = 1)
        df["prediction"] = model.predict(X_test)
        df[['Name', 'prediction', "SCORE", "Lig"]].to_csv(csv_file.replace("_feature.csv", "_score.csv"), index = False)
    except Exception as e:
        with open(args.log, "a") as f_log:
            f_log.write(csv_file + "{} \n".format(e))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Obtain RmsdXNA score of the poses.")
    parser.add_argument("-folder_dock", type=str, default="example/docking", help = "Folder containing ligand features .csv files for prediction")
    parser.add_argument("-cutoff", type=float, default=8, help = "Cutoff distance for feature extraction")
    parser.add_argument("-log", type=str, default="example/error_getscore.log", help = ".log file to record error")
    parser.add_argument("-o", type=str, default="allscores.csv", help = "compiled output filename for scores")
    parser.add_argument("-model", type=str, default="xgboostfull.model", help = "Converted ligand atom_name")
    parser.add_argument("-column", type=str, default="feature_column.json", help = "json file containing dataset columns")
    parser.add_argument("-ncpus", type=int, default=1, help= "no. of CPU used to run jobs in parallel")
    args = parser.parse_args()
    
    model = xgb.XGBRegressor()
    os.environ["OMP_THREAD_LIMIT"] = "1"
    model.load_model(args.model)

    with open(args.column, "r") as f_dict:
        label_value_dict = json.load(f_dict)
        f_dict.close()
    feature_col = sorted(label_value_dict.keys())
    
    ligand_csv = glob(args.folder_dock + "/*_feature.csv")
    pool = Pool(args.ncpus)
    pool.map(process, ligand_csv)
    
    os.system("for f in {}/*_score.csv; do sort -t ',' -k 2 -n $f | head -n 2 | tail -n 1 >> {}; done".format(args.folder_dock, args.o))

    

