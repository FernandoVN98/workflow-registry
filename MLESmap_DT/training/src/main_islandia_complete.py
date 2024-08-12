import json
from dislib.preprocessing import MinMaxScaler

from modules.simulation_phase import Simulation
from modules.cybershake_workflows import cybershake_workflow
from modules.data_manager.base import DataManager
from modules.data_target import Disk
from modules.preprocessing import Preprocessing
from modules.DigitalTwin import DigitalTwin
from modules.training import ModelSelection
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.tree import DecisionTreeRegressor
from pycompss.api.api import compss_wait_on, compss_wait_on_directory, compss_wait_on_file, compss_barrier
from auxiliar_function import _copy_output_cybershake, _correct_output_cybershake
from modules.preprocessing_workflow import merged_steps_preprocess_mlesmap
from os import walk
import sys

def read_configuration_file(configuration_file):
    with open(configuration_file, 'r') as f:
        config_file = json.load(f)
    cybershake_input_file = config_file["cybershake_file"]
    path_to_folders = config_file["output_simulation"]
    output_folder = config_file["output_folder"]
    table_CS_Ruptures = config_file["CS_Ruptures_file"]
    return cybershake_input_file, path_to_folders, output_folder, table_CS_Ruptures

def main(configuration_file):
    cybershake_input_file, path_to_folders, output_folder, table_CS_Ruptures = read_configuration_file(configuration_file)
    sim = Simulation(cybershake_workflow, format_result=merged_steps_preprocess_mlesmap)
    args_format_results = {"path_to_folders": path_to_folders, 
            "output_folder": output_folder, 
            "table_CS_Ruptures": table_CS_Ruptures}
    sim.set_args_format_results(args_format_results)
    ml_train = sim.execute({cybershake_input_file})
    # Disk and Data Manager Reads input Data
    disk = Disk("/gpfs/scratch/bsc19/bsc019756/Iceland_geoinq/Files/output_folder", read_file="ml_train_dataset.csv")
    dm = DataManager(disk)
    #data_readed = dm.get_data()
    minmax_scaler = MinMaxScaler(feature_range=(0, 1))
    preprocess = Preprocessing(method=minmax_scaler)
    model_sel = ModelSelection(scoring=r2_score)
    model_sel.set_models([DecisionTreeRegressor()])
    model_sel.set_paramaters_models(
        [{"criterion": ["squared_error", "friedman_mse"], "max_depth": [2, 3, 5, 7], "random_state": [0]}])
    # Digital Twin Object
    dt = DigitalTwin(simulator=sim, data_manager=dm, data_preprocessing=preprocess, training=model_sel)
    best_model_instance = dt.execute()
    dt.visualize_results_training()
    return best_model_instance

if __name__ == '__main__':
    configuration_file = sys.argv[1]
    main(configuration_file)

