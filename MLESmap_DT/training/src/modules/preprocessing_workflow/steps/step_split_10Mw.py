from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import pandas as pd


@constraint(computing_units="${ComputingUnits}")
@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_10Mw_compss(PSA, name_out, seed=42):
    Database = pd.read_csv(PSA)
    MwData = pd.DataFrame()
    unique_magnitudes = Database['Magnitude'].unique()  # Magnitude unique values
    for magnitude in unique_magnitudes:  # Iterate through each unique value of 'Magnitude'

        subset = Database[
            Database['Magnitude'] == magnitude]  # Filter Database to get the rows with the current magnitude
        num_samples = int(0.20 * len(subset))
        random_subset = subset.sample(n=num_samples, random_state=seed)  # Seed

        MwData = pd.concat([MwData, random_subset])

    MwData.to_csv(name_out, index=False)


@constraint(computing_units="${ComputingUnits}")
@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_10Random_compss(PSA, name_out, seed=42):
    DataIn = pd.read_csv(PSA)
    dfRandom = DataIn.sample(frac=0.10, random_state=seed)

    dfRandom.to_csv(name_out, index=False)


def step_split_10Mw(path_to_folder, PSA, random=False):
    if random:
        name_out = path_to_folder + "/10_source_PSA_" + PSA.split("_", -1)[-2] + "_rup_rupvar_rand.csv"
    else:
        name_out = path_to_folder + "/10_source_PSA_" + PSA.split("_", -1)[-2] + "_rup_rupvar_10Mw.csv"
    if random:
        split_10Random_compss(PSA, name_out)
    else:
        split_10Mw_compss(PSA, name_out)
    return name_out
