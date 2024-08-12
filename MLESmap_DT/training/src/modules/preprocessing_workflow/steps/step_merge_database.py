from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import os
import pandas as pd


def step_merge_database(path_to_folders, psa_1s, psa_2s, psa_3s, psa_5s, psa_7s, psa_10s, output_folder=None):
    if output_folder is not None:
        folder_Out = output_folder + "/step_4"
    else:
        folder_Out = path_to_folders + "/step_4"
    if not os.path.exists(folder_Out):
        os.mkdir(folder_Out)
    PSA_1s = folder_Out + '/PSA_1_SS.csv'
    PSA_2s = folder_Out + '/PSA_2_SS.csv'
    PSA_3s = folder_Out + '/PSA_3_SS.csv'
    PSA_5s = folder_Out + '/PSA_5_SS.csv'
    PSA_7s = folder_Out + '/PSA_7_SS.csv'
    PSA_10s = folder_Out + '/PSA_10_SS.csv'
    merge_psas_step_4(psa_1s, PSA_1s)
    merge_psas_step_4(psa_2s, PSA_2s)
    merge_psas_step_4(psa_3s, PSA_3s)
    merge_psas_step_4(psa_5s, PSA_5s)
    merge_psas_step_4(psa_7s, PSA_7s)
    merge_psas_step_4(psa_10s, PSA_10s)
    return PSA_1s, PSA_2s, PSA_3s, PSA_5s, PSA_7s, PSA_10s


@constraint(computing_units="${ComputingUnits}")
@task(psas=COLLECTION_FILE_IN, name_out=FILE_OUT, returns=1)
def merge_psas_step_4(psas, name_out):
    PSA = pd.DataFrame()
    print("PSAAAS")
    print(psas)
    for file in psas:
        f1 = pd.read_csv(file)
        print("FILLLEEEEEEE")
        print(f1, flush=True)
        PSA = pd.concat([PSA, f1])
    print("CONCATTED DATAFRAME")
    print(PSA, flush=True)
    PSA.to_csv(name_out)
