from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import pandas as pd
import numpy as np


@constraint(computing_units="${ComputingUnits}")
@task(input_file=FILE_IN, output_file=FILE_OUT)
def adapt_dataset_to_dislib(input_file, output_file):
    dfDataSet0 = pd.read_csv(input_file)

    # Generate features of the model
    dfFeatures = pd.DataFrame()
    dfFeatures["Site_Lat"] = dfDataSet0["Site_Lat"]
    dfFeatures["Site_Lon"] = dfDataSet0["Site_Lon"]
    dfFeatures["Magnitude"] = dfDataSet0["Magnitude"]
    dfFeatures["Hypocenter_Lat"] = dfDataSet0["Hypocenter_Lat"]
    dfFeatures["Hypocenter_Lon"] = dfDataSet0["Hypocenter_Lon"]
    dfFeatures["Hypocenter Depth"] = dfDataSet0["Hypocenter_Depth"]
    dfFeatures["EuclideanDistance"] = dfDataSet0["Geodesic_distance"]#dfDataSet0["EuclideanDistance_EQ-Sta"]
    dfFeatures["Azimuth"] = dfDataSet0["Azimuth"]#dfDataSet0["Azim_EQ-Sta"]
    #dfFeatures["Lower Intensity Value 10s"] = np.log10(dfDataSet0["PSA"])#np.log10(dfDataSet0["Intensity Value"])
    #dfFeatures["Target Intensity Value 3s"] = np.log10(dfDataSet0["PSA_max"])#np.log10(dfDataSet0["Intensity Value 3s"])
    dfFeatures["Intensity Value"] = np.log10(dfDataSet0["PSA_max"])
    dfFeatures.to_csv(output_file, index=True)


def generate_dislib_dataset(input_file, output_file):
    adapt_dataset_to_dislib(input_file, output_file)
    return output_file
