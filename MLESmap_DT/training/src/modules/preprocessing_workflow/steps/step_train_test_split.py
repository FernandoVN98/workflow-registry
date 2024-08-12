from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import pandas as pd


@constraint(computing_units="${ComputingUnits}")
@task(path_to_source_rupvar=FILE_IN, psa_to_split=FILE_IN, data_test=FILE_OUT, data_train=FILE_OUT)
def train_test_split_tasks(path_to_source_rupvar, psa_to_split, data_test, data_train):
    # PASO 1 TEST
    df10Selected = pd.read_csv(path_to_source_rupvar)
    Database = pd.read_csv(psa_to_split)

    dfEQ = pd.merge(Database, df10Selected, on=['Source_ID', 'Rupture_ID', 'Rupture_Variation_ID', 'Magnitude'],
                    how='inner')
    dfEQ.to_csv(data_test, index=False)

    # PASO 2 TRAIN
    chunk_size = 600000
    chunks = pd.read_csv(psa_to_split, chunksize=chunk_size)
    DataOut1 = pd.read_csv(data_test)

    TrainChunks = []
    for chunk in chunks:
        chunk_df = chunk
        pd.merge(chunk_df, DataOut1, on=['Source_ID', 'Rupture_ID', 'Rupture_Variation_ID', 'Magnitude'], how='left')
        merged = chunk_df.merge(DataOut1, how='left', indicator=True)

        dfNotEQ = merged[merged['_merge'] == 'left_only']  # in the merge column, only those listed as left_only
        dfNotEQ = dfNotEQ.drop('_merge', axis=1)

        TrainChunks.append(dfNotEQ)

    Train = pd.concat(TrainChunks, ignore_index=True)
    Train.to_csv(data_train, index=False)

    return data_train, data_test


def train_test_split(path_to_source_rupvar, path_to_folder, PSA, name_train_out, name_test_out):
    name_out_train = path_to_folder + name_train_out
    name_out_test = path_to_folder + name_test_out
    train_test_split_tasks(path_to_source_rupvar, PSA, name_out_test, name_out_train)
    return name_out_train, name_out_test
