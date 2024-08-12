from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import os
from os import walk
import pandas as pd
import utm


def step_srf_extraction(path_to_folders, output_folder=None, save_intermediate_files=False):
    f_Folders = []
    f_files_out = []
    for (dirpath, dirnames, filenames) in walk(path_to_folders):
        f_Folders.extend(dirnames)
        break
    for nk in f_Folders[:]:
        f_files_out.append([])
        path_to_file = path_to_folders + "/" + nk + '/'
        f_Stations = []
        for (dirpath, dirnames, filenames) in walk(path_to_file):
            f_Stations.extend(dirnames)
            break
        for i in f_Stations:  # Loop that goes through the stations
            srf_in_folder = []
            stations = []
            stations_compss = []
            path_Stat = path_to_file + i + '/'
            if output_folder is not None:
                folder_Out = output_folder + "/step_2"
                if not os.path.exists(output_folder):
                    os.mkdir(output_folder)
            else:
                folder_Out = path_Stat + "step_2"
                if not os.path.exists(path_Stat):
                    os.mkdir(path_Stat)
            if not os.path.exists(folder_Out):
                os.mkdir(folder_Out)
            f_stat = i.split('_')
            stat = f_stat[1]
            name_Out = folder_Out + '/srf_extraction_' + stat + '_src.csv'
            name_compss = 'file_srf_extraction_' + stat + '_src.csv'
            flag = os.path.exists(path_Stat + 'post-processing')
            srf_out = None
            if flag == True:
                srf_out = compute_srf(path_Stat, folder_Out, name_Out)
            srf_in_folder.append(srf_out)
            stations.append(name_Out)
            stations_compss.append(name_compss)
            for out_station, name_station, name_compss in zip(srf_in_folder, stations, stations_compss):
                output_step_srf_extraction(out_station, name_Out)
                f_files_out[-1].append([])
                f_files_out[-1][-1].append(name_Out)
    return f_files_out

def compute_srf(path_Stat,folder_Out,name_Out):
    f_Files = []    # Contains the seismogram files
    for (dirpath,dirnames, filenames) in walk(path_Stat+'post-processing/'):
        f_Files.extend(filenames)
        break
    f_SRF = []
    for ik in f_Files:
        try:
            f2 = ik.split(".")
            if f2[3] == 'srf':
                f_SRF.append(path_Stat+'post-processing/'+ik)
        except:
            R = 1+1
    if f_SRF != []:
        srf_out_station = srf_station_extraction(f_SRF)
    else:
        print('station-no-post-pro',folder_Out)
        srf_out_station = pd.DataFrame()
    return srf_out_station

def srf_station_extraction(f_SRF):
    srf_out_station = []
    for srf_file in f_SRF:
        srf_out = srf_extraction(srf_file)
        srf_out_station.append(srf_out)
    return srf_out_station


@constraint(computing_units="${ComputingUnits}")
@task(srf_file2=FILE_IN, returns=1)
def srf_extraction(srf_file2):
    srf_temp =  os.path.basename(srf_file2)
    srf_file = srf_temp
    my_file = pd.read_table(srf_file2)
    df1 = my_file.iloc[1]              # first row data extraction
    df2 = df1.str.split( )          # split first row data by space
    df3 = my_file.iloc[2]              # second row data extraction
    df4 = df3.str.split( )
    data = {'ELON': float(df2.iloc[0][0]),
            'ELAT': float(df2.iloc[0][1]),
            'DTOP': float(df4.iloc[0][2]),
            'SHYP': float(df4.iloc[0][3]),
            'DHYP': float(df4.iloc[0][4])}
    x_hypo, y_hypo, utm_number, utm_letter = utm.from_latlon((data['ELAT']), (data['ELON']))
    y_hypo += 1000.0 * data['SHYP']
    Lat_hypo, Lon_hypo = utm.to_latlon(x_hypo, y_hypo,utm_number, utm_letter)
    Depth_hypo = data['DTOP'] + data['DHYP']
    f = srf_file.split("_")
    f2 = f[3].split(".")
    f3 = f[4].split("s")
    data_final = {'Source_ID': int(f[2]),
                  'Rupture_ID': int(f2[0]),
                  'Rupture_Variation_ID':int(f3[1]),
                  'Lat_hypo': float(Lat_hypo),
                  'Lon_hypo': float(Lon_hypo),
                  'Depth_hypo': float(Depth_hypo)}
    data_final.update(data)
    return data_final


@constraint(computing_units="${ComputingUnits}")
@task(out_station=COLLECTION_IN, outFile=FILE_OUT)
def output_step_srf_extraction(out_station, outFile):
    pd_to_store = pd.DataFrame(out_station)
    pd_to_store.to_csv(outFile)
