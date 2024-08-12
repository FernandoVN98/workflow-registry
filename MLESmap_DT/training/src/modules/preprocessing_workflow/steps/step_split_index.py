from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import pandas as pd


@constraint(computing_units="${ComputingUnits}")
@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_index_compss_random(PSA, name_out):
    Site_Rup_RupVar = pd.DataFrame()
    Database = pd.read_csv(PSA)

    Site_Rup_RupVar['Source_ID'] = Database['Source_ID']
    Site_Rup_RupVar['Rupture_ID'] = Database['Rupture_ID']
    Site_Rup_RupVar['Rupture_Variation_ID'] = Database['Rupture_Variation_ID']

    Site_Rup_RupVar.to_csv(name_out, index=False)


@constraint(computing_units="${ComputingUnits}")
@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_index_compss(PSA, name_out):
    Site_Rup_RupVar = pd.DataFrame()
    Database = pd.read_csv(PSA)

    Site_Rup_RupVar['Source_ID'] = Database['Source_ID']
    Site_Rup_RupVar['Rupture_ID'] = Database['Rupture_ID']
    Site_Rup_RupVar['Rupture_Variation_ID'] = Database['Rupture_Variation_ID']
    Site_Rup_RupVar['Magnitude'] = Database['Magnitude']

    Site_Rup_RupVar.to_csv(name_out, index=False)


def step_split_index(output_folder, PSA, random=False):
    name_out = output_folder + "/10_Source_Rup_RupVar_Mw.csv"
    if random:
        split_index_compss_random(PSA, name_out)
    else:
        split_index_compss(PSA, name_out)
    return name_out
