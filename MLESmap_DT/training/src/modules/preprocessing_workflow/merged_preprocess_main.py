
from modules.preprocessing_workflow.rotinv_maxavg import *
from modules.preprocessing_workflow.read_seis_TEST import *
from pycompss.api.parameter import *
import numpy as np
import pandas as pd
import utm
import os
from os import walk
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import sys

import numpy as np
import pandas as pd
from geopy.distance import geodesic
from geopy import Point
import math
import os
from os import walk
from pycompss.api.task import task
from pycompss.api.parameter import *
from pycompss.api.api import compss_wait_on, compss_wait_on_file, compss_barrier
import sys

from dislib.preprocessing import MinMaxScaler

from modules.data_manager.base import DataManager
from modules.data_target import Disk
from modules.preprocessing import Preprocessing
from modules.DigitalTwin import DigitalTwin
from modules.training import ModelSelection
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error
from sklearn.tree import DecisionTreeRegressor
import sys


def obtain_stations_name(path_to_folders, output_folder):#TODO CREO QUE HAY QUE BORRARLO
    f_Folders = []
    stations = []
    for (dirpath, dirnames, filenames) in walk(path_to_folders):
        f_Folders.extend(dirnames)
        break
    for nk in f_Folders[:]:
        path_to_file = path_to_folders + "/" + nk + '/'
        f_Stations = []
        for (dirpath, dirnames, filenames) in walk(path_to_file):
            f_Stations.extend(dirnames)
            break
        stations.append([])
        for i in f_Stations:
            path_Stat = path_to_file + i + '/'
            folder_Out = path_Stat + output_folder
            f_stat = i.split('_')
            stat = f_stat[1]
            name_Out = folder_Out + '/srf_extraction_' + stat + '_src.csv'
            stations[-1].append(name_Out)
    return stations

#MAIN
def psa_start(path_to_folders,output_folder):
    f_Folders = []
    for (dirpath, dirnames, filenames) in walk(path_to_folders):
        f_Folders.extend(dirnames)
        break
    for nk in f_Folders[0:1]:
        path_to_file = path_to_folders+'/'+nk+'/'
        f_Stations = []
        for (dirpath, dirnames,filenames) in walk(path_to_file):
            f_Stations.extend(dirnames)
            break
        for i in f_Stations:
            path_Stat = path_to_file+i+'/PSAoutputExtraction/'
            path_srf = path_to_file+i+'/PSAsrfExtraction/'
            if not os.path.exists(path_Stat):
                print('There is not a PSAoutputExtraction folder in :' + path_to_file+i)
                continue

            folder_Out = path_to_file + i + '/' + output_folder
            if not os.path.exists(folder_Out):
                os.mkdir(folder_Out)

            f_Files = []
            for (dirpath, dirnames, filenames) in walk(path_Stat):
                f_Files.extend(filenames)
                break
            CyberS_sites_temp = (
                    "./inputFiles/CyberShake_Sites.csv"
                )
            f_stat = i.split("_")
            stat = f_stat[1]

            name_Out1 = folder_Out + '/PSA_1_SS_' + stat + '.csv'
            name_Out2 = folder_Out + '/PSA_2_SS_' + stat + '.csv'
            name_Out3 = folder_Out + '/PSA_3_SS_' + stat + '.csv'
            name_Out5 = folder_Out + '/PSA_5_SS_' + stat + '.csv'
            name_Out7 = folder_Out + '/PSA_7_SS_' + stat + '.csv'
            name_Out10 = folder_Out + '/PSA_10_SS_' + stat + '.csv'

            SS_psa(path_Stat,path_srf,f_Files,CyberS_sites_temp,name_Out1,name_Out2,name_Out3,name_Out5,name_Out7,name_Out10)


@task(db=COLLECTION_IN, name_Out=FILE_OUT)
def PSAs(db, name_Out):
    PSA = pd.DataFrame(dtype=np.float64)
    for db_ss in db:
        PSA = pd.concat([PSA, db_ss])
    PSA['Haversine_distance'] = 0
    PSA['Azimuth'] = 0
    PSA['Geodesic_distance'] = 0
    for r in range(len(PSA)):
        x1, y1, z1 = PSA['Hypocenter_Lat'].iloc[r], PSA['Hypocenter_Lon'].iloc[r], \
        PSA['Hypocenter_Depth'].iloc[r]
        start = Point(PSA['Hypocenter_Lat'].iloc[r], PSA['Hypocenter_Lon'].iloc[r])
        x2, y2, z2 = PSA['Site_Lat'].iloc[r], PSA['Site_Lon'].iloc[r], 0
        end = Point(PSA['Site_Lat'].iloc[r], PSA['Site_Lon'].iloc[r])

        distance = (geodesic(start, end).km)
        azimuth = (math.degrees(math.atan2((x2 - x1), (y2 - y1))))
        R = 6371
        lat1 = math.radians(x1)
        lon1 = math.radians(y1)
        lat2 = math.radians(x2)
        lon2 = math.radians(y2)
        dlat = lat2 - lat1
        dlon = lon2 - lon1
        ddepth = z2 - z1
        a = math.sin(dlat / 2) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        distance_surface = R * c
        distance1 = math.sqrt(distance_surface ** 2 + ddepth ** 2)
        PSA.loc[r, 'Haversine_distance'] = distance1
        PSA.loc[r, 'Azimuth'] = azimuth
        PSA.loc[r, 'Geodesic_distance'] = distance
    PSA = PSA.dropna()
    PSA.to_csv(name_Out)
    #return PSA


@task(path_Stat=FILE_IN, path_srf=FILE_IN, CyberS_sites_temp=FILE_IN)
def task_f_Files(file, path_Stat, CyberS_sites_temp, path_srf):
    ij = file
    ij_split = ij.split('_')
    SS_csv = int(ij_split[4])
    src_csv = int(ij_split[6])
    psa_stat_temp = path_Stat + file
    f = pd.read_csv(psa_stat_temp)
    psa_db = pd.DataFrame()

    for k in range(len(f['Unnamed: 0'])):
        l = f['Unnamed: 0'].iloc[k]
        ik_split = l.split('_')
        psa_number = int(ik_split[2])

        psa_db = pd.concat([psa_db, pd.DataFrame({'SS': SS_csv, 'Source_ID': src_csv, 'PSA': psa_number,
                                                  'PSA_max': f['PSA_max'].iloc[k],
                                                  'Rupture_Variation_ID': f['rv'].iloc[k]}, index=[
            0])])

        f3 = pd.read_csv(CyberS_sites_temp)
        data_site = {
            'SS': (f3['CS_Site_ID'] - 1000),
            'Site_Lat': (f3['CS_Site_Lat']),
            'Site_Lon': (f3['CS_Site_Lon'])}
        site_db = pd.DataFrame(data_site)

        db_ss2 = pd.merge(psa_db, site_db, on='SS')

        psa_1_db = db_ss2.query('PSA == 1')
        psa_2_db = db_ss2.query('PSA == 2')
        psa_3_db = db_ss2.query('PSA == 3')
        psa_5_db = db_ss2.query('PSA == 5')
        psa_7_db = db_ss2.query('PSA == 7')
        psa_10_db = db_ss2.query('PSA == 10')

        f_Test = []
        for (dirpath, dirnames, filenames) in walk(path_srf):
            f_Test.extend(filenames)

        for m in f_Test:
            im = m
            im_split = im.split('_')
            SS_Test = int(im_split[2])
            f2 = pd.read_csv(path_srf + m)
            data_srf = {
                'SS': SS_Test,
                'Source_ID': (f2['Source_ID']),
                'Rupture_ID': (f2['Rupture_ID']),
                'Rupture_Variation_ID': (f2['Rupture_Variation_ID']),  # ,Rupture_Variation_ID    (+1)
                'Hypocenter_Lat': (f2['Lat_hypo']),
                'Hypocenter_Lon': (f2['Lon_hypo']),
                'Hypocenter_Depth': (f2['Depth_hypo']).round(2),
                'Magnitude': (f2['Rupture_ID'] * 0.1).round(2)
            }
            srf_db = pd.DataFrame(data_srf)

    db_ss1_1s = pd.merge(psa_1_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_2s = pd.merge(psa_2_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_3s = pd.merge(psa_3_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_5s = pd.merge(psa_5_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_7s = pd.merge(psa_7_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_10s = pd.merge(psa_10_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    return db_ss1_1s, db_ss1_2s, db_ss1_3s, db_ss1_5s, db_ss1_7s, db_ss1_10s

def SS_psa(path_Stat,path_srf,f_Files,CyberS_sites_temp,name_Out1,name_Out2,name_Out3,name_Out5,name_Out7,name_Out10):
    PSA_1s = pd.DataFrame(dtype=np.float64)
    PSA_2s = pd.DataFrame(dtype=np.float64)
    PSA_3s = pd.DataFrame(dtype=np.float64)
    PSA_5s = pd.DataFrame(dtype=np.float64)
    PSA_7s = pd.DataFrame(dtype=np.float64)
    PSA_10s = pd.DataFrame(dtype=np.float64)
    srf_db = pd.DataFrame(dtype=np.float64)
    f_Files = []
    for (dirpath, dirnames, filenames) in walk(path_Stat):
        f_Files.extend(filenames)
        break
    num_file = 0
    db_1s = []
    db_2s = []
    db_3s = []
    db_5s = []
    db_7s = []
    db_10s = []
    for j in f_Files:
        num_file += 1
        db_ss1_1s, db_ss1_2s, db_ss1_3s, db_ss1_5s, db_ss1_7s, db_ss1_10s = task_f_Files(j, path_Stat, CyberS_sites_temp, path_srf)
        db_1s.append(db_ss1_1s)
        db_2s.append(db_ss1_2s)
        db_3s.append(db_ss1_3s)
        db_5s.append(db_ss1_5s)
        db_7s.append(db_ss1_7s)
        db_10s.append(db_ss1_10s)
    '''PSAs(db_1s, name_Out1)
    PSAs(db_2s, name_Out2)
    PSAs(db_3s, name_Out3)
    PSAs(db_5s, name_Out5)
    PSAs(db_7s, name_Out7)
    PSAs(db_10s, name_Out10)'''

def srf_station_extraction(f_SRF):
    srf_out_station = []
    for srf_file in f_SRF:
        srf_out = srf_extraction(srf_file)
        srf_out_station.append(srf_out)
    return srf_out_station

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
def step_1(path_to_folders, output_folder=None, save_intermediate_files=False):
    f_Folders = []
    f_files_out = []
    if output_folder is not None:
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)
    for (dirpath, dirnames, filenames) in walk(path_to_folders):
        f_Folders.extend(dirnames)
        break
    for nk in f_Folders:
        f_files_out.append([])
        path_to_file = path_to_folders+"/"+nk+'/'
        try:
            f_Stations = []
            for (dirpath, dirnames, filenames) in walk(path_to_file):
                f_Stations.extend(dirnames)
                break
            for i in f_Stations:  # Loop that goes through the stations
                f_files_out[-1].append([])
                path_Stat = path_to_file + i +'/'
                if output_folder is not None:
                    folder_Out = output_folder+"/step_1"
                else:
                    folder_Out = path_Stat+"step_1"
                try:
                    if not os.path.exists(folder_Out):
                        os.mkdir(folder_Out)
                    f_Files = []  # Contains the seismogram files
                    for (dirpath, dirnames, filenames) in walk(path_Stat+'post-processing/'):
                        f_Files.extend(filenames)
                        break
                    f_Seismograms = []
                    for ik in f_Files:
                        try:
                            f2 = ik.split(".")
                            if len(f2) > 1:
                                if f2[1] == 'grm':
                                    f_Seismograms.append(path_Stat+'post-processing/' + ik)
                        except Exception as excep:
                            print('Exception: ' + str(excep))
                    for seismogram_file in f_Seismograms:
                        j = os.path.basename(seismogram_file)
                        f1 = j.split("_")
                        src = f1[4]
                        stat = f1[2]
                        f3 = f1[5].split(".")
                        rup = f3[0]
                        compss_file = folder_Out+'/RotI_PGA_PSA_stat_'+stat+'_src_'+src+'_rup_'+rup+'.csv'
                        seismogram_station_extraction_only_one_file( seismogram_file, compss_file, write_pandas_df=save_intermediate_files)
                        f_files_out[-1][-1].append(compss_file)
                except Exception as excep:
                    print('ERROR_STATION', i)
                    print(excep)
        except:
            print('Error_folder', path_to_file)
    return f_files_out

@task(f_Seismogram=FILE_IN, file_out=FILE_OUT)
def seismogram_station_extraction_only_one_file(f_Seismogram, file_out, write_pandas_df=False):
    j = os.path.basename(f_Seismogram)
    f1 = j.split("_")
    src = f1[4]
    stat = f1[2]
    f3 = f1[5].split(".")
    rup = f3[0]
    Mw = int(rup)*0.1
    g = 0.01
    damp = 0.05

    #nameOut = folder_Out+'/RotI_PGA_PSA_stat_'+stat+'_src_'+src+'_rup_'+rup+'.csv'
    seismogram_data_extraction(f_Seismogram,Mw,rup,src,stat,g,damp, file_out)

def seismogram_data_extraction(grm_file, Mw,rup,src,stat,g,damp, file_out):
    periods  = np.asarray((1, 2, 3, 5, 7.5, 10))
    seis = Seismogram(grm_file)
    seis.readData()
    dfUmax = pd.DataFrame()
    index_labels = []

    for k in range(0, seis.num_rvs):
        nt = seis.nt
        dt = seis.dt
        (x_data, y_data) = seis.data_dict[k]
        x_data = np.asarray(x_data)
        y_data = np.asarray(y_data)
        Ax = (1/dt)*(x_data[1:nt] - x_data[0:nt-1]); Ay = (1/dt)*(y_data[1:nt] - y_data[0:nt-1])
        Umax = rotinv_maxavg(Ax*g,Ay*g,dt,periods,damp)
        dir1 = 'rv_'+str(k+1)+'_1_s'
        dir2 = 'rv_'+str(k+1)+'_2_s'
        dir3 = 'rv_'+str(k+1)+'_3_s'
        dir4 = 'rv_'+str(k+1)+'_5_s'
        dir5 = 'rv_'+str(k+1)+'_7_s'
        dir6 = 'rv_'+str(k+1)+'_10_s'
        index_labels.append(dir1)
        index_labels.append(dir2)
        index_labels.append(dir3)
        index_labels.append(dir4)
        index_labels.append(dir5)
        index_labels.append(dir6)
        Umax['Mw'] = Mw*np.ones(len(Umax))
        Umax['src'] = int(src)*np.ones(len(Umax))
        Umax['rv'] = k*np.ones(len(Umax))+1
        Umax = Umax.set_axis(index_labels, axis='index')
        dfUmax = pd.concat([dfUmax, Umax])
        index_labels = []
    dfUmax.to_csv(file_out)

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

#@task(out_station=COLLECTION_IN)
#def store_df_step_2a(out_station, name_station):
#    pd_to_store = pd.DataFrame(out_station)
#    pd_to_store.to_csv(name_station)

@task(out_station=COLLECTION_IN, outFile=FILE_OUT)
def output_step_2(out_station, outFile):
    pd_to_store = pd.DataFrame(out_station)
    pd_to_store.to_csv(outFile)

def step_2(path_to_folders, output_folder=None, save_intermediate_files=False):
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
        srf_in_folder = []
        stations = []
        stations_compss = []
        for i in f_Stations:  # Loop that goes through the stations
            srf_in_folder = []
            stations = []
            stations_compss = []
            #f_files_out[-1].append([])
            path_Stat = path_to_file + i + '/'
            if output_folder is not None:
                folder_Out = output_folder + "/step_2"
            else:
                folder_Out = path_Stat + "step_2"
            if not os.path.exists(folder_Out):
                os.mkdir(folder_Out)
            f_stat = i.split('_')
            print("IIII")
            print(i)
            print(f_stat)
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
                output_step_2(out_station, name_Out)
                f_files_out[-1].append([])
                f_files_out[-1][-1].append(name_Out)
        '''if save_intermediate_files:
            for out_station, name_station, name_compss in zip(srf_in_folder, stations, stations_compss):
                store_df_step_2a(out_station, name_station)
                #out_file = "file.csv"
                output_step_2(out_station, name_compss)
                f_files_out[-1].append(name_compss)#f_files_out[-1][-1].append(name_station)
        else:
            for out_station, name_station, name_compss in zip(srf_in_folder, stations, stations_compss):
                #out_file = "file.csv"
                output_step_2(out_station, name_station)
                f_files_out[-1].append(name_station)#f_files_out[-1][-1].append(name_station)'''
    return f_files_out


@task(file=FILE_IN, CyberS_sites_temp=FILE_IN, path_srf=COLLECTION_FILE_IN, returns=6)
def task_f_Files_compss(file, CyberS_sites_temp, path_srf):
    ij = file.split('/')[-1]
    ij_split = ij.split('_')
    SS_csv = int(ij_split[4])
    src_csv = int(ij_split[6])
    #psa_stat_temp = path_Stat + file
    f = pd.read_csv(file)
    psa_db = pd.DataFrame()

    for k in range(len(f['Unnamed: 0'])):
        l = f['Unnamed: 0'].iloc[k]
        ik_split = l.split('_')
        psa_number = int(ik_split[2])

        psa_db = pd.concat([psa_db, pd.DataFrame({'SS': SS_csv, 'Source_ID': src_csv, 'PSA': psa_number,
                                                  'PSA_max': f['PSA_max'].iloc[k],
                                                  'Rupture_Variation_ID': f['rv'].iloc[k]}, index=[
            0])])

        f3 = pd.read_csv(CyberS_sites_temp)
        data_site = {
            'SS': 179,#(f3['CS_Site_ID'] - 1000),
            'Site_Lat': (f3['CS_Site_Lat']),
            'Site_Lon': (f3['CS_Site_Lon'])}
        site_db = pd.DataFrame(data_site)

        db_ss2 = pd.merge(psa_db, site_db, on='SS')# ESTO ESTABA BIEN NO SÉ QUE PASA

        psa_1_db = db_ss2.query('PSA == 1')
        psa_2_db = db_ss2.query('PSA == 2')
        psa_3_db = db_ss2.query('PSA == 3')
        psa_5_db = db_ss2.query('PSA == 5')
        psa_7_db = db_ss2.query('PSA == 7')
        psa_10_db = db_ss2.query('PSA == 10')

        f_Test = []

        for m in path_srf:
            im = m
            im_split = im.split('_')
            SS_Test = int(im_split[-2])
            f2 = pd.read_csv(m)
            data_srf = {
                'SS': SS_Test,
                'Source_ID': (f2['Source_ID']),
                'Rupture_ID': (f2['Rupture_ID']),
                'Rupture_Variation_ID': (f2['Rupture_Variation_ID']),  # ,Rupture_Variation_ID    (+1)
                'Hypocenter_Lat': (f2['Lat_hypo']),
                'Hypocenter_Lon': (f2['Lon_hypo']),
                'Hypocenter_Depth': (f2['Depth_hypo']).round(2),
                'Magnitude': (f2['Rupture_ID'] * 0.1).round(2)
            }
            srf_db = pd.DataFrame(data_srf)
    db_ss1_1s = pd.merge(psa_1_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_2s = pd.merge(psa_2_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_3s = pd.merge(psa_3_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_5s = pd.merge(psa_5_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_7s = pd.merge(psa_7_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    db_ss1_10s = pd.merge(psa_10_db, srf_db, on=['Rupture_Variation_ID', 'Source_ID'])
    return db_ss1_1s, db_ss1_2s, db_ss1_3s, db_ss1_5s, db_ss1_7s, db_ss1_10s

def SS_psa_compss(path_Stat,path_srf,CyberS_sites_temp,name_Out1,name_Out2,name_Out3,name_Out5,name_Out7,name_Out10):
    db_1s = []
    db_2s = []
    db_3s = []
    db_5s = []
    db_7s = []
    db_10s = []
    #real_path_srf = path_srf
    xj = path_srf.split('/')[:-1]
    real_path_srf = "/" + "/".join(xj)
    for j in path_Stat:
        db_ss1_1s, db_ss1_2s, db_ss1_3s, db_ss1_5s, db_ss1_7s, db_ss1_10s = task_f_Files_compss(j, CyberS_sites_temp, real_path_srf)
        db_1s.append(db_ss1_1s)
        db_2s.append(db_ss1_2s)
        db_3s.append(db_ss1_3s)
        db_5s.append(db_ss1_5s)
        db_7s.append(db_ss1_7s)
        db_10s.append(db_ss1_10s)
    PSA_1s = PSAs(db_1s, name_Out1)
    PSA_2s = PSAs(db_2s, name_Out2)
    PSA_3s = PSAs(db_3s, name_Out3)
    PSA_5s = PSAs(db_5s, name_Out5)
    PSA_7s = PSAs(db_7s, name_Out7)
    PSA_10s = PSAs(db_10s, name_Out10)


    PSA_1s.to_csv(name_Out1)
    PSA_2s.to_csv(name_Out2)
    PSA_3s.to_csv(name_Out3)
    PSA_5s.to_csv(name_Out5)
    PSA_7s.to_csv(name_Out7)
    PSA_10s.to_csv(name_Out10)

def step_3(first_step_files, second_step_files, path_fo_file, save_intermediate_files=False):

    for file_folder_output_PSA, file_folder_output_srf in zip(first_step_files, second_step_files):
        for file_station_psa, file_station_output in zip(file_folder_output_PSA, file_folder_output_srf):
            path_out = file_station_output.split("/", -1)[:-2]
            i_split = file_station_output.split("/", -1)[-1]
            f_stat = i_split.split("_")
            stat = f_stat[1]
            folder_Out = "/" + "/".join(path_out) + "/step_3"
            if not os.path.exists(folder_Out):
                os.mkdir(folder_Out)
            name_Out1 = folder_Out + '/PSA_1_SS_' + stat + '.csv'
            name_Out2 = folder_Out + '/PSA_2_SS_' + stat + '.csv'
            name_Out3 = folder_Out + '/PSA_3_SS_' + stat + '.csv'
            name_Out5 = folder_Out + '/PSA_5_SS_' + stat + '.csv'
            name_Out7 = folder_Out + '/PSA_7_SS_' + stat + '.csv'
            name_Out10 = folder_Out + '/PSA_10_SS_' + stat + '.csv'
            CyberS_sites_temp = (
                "./inputFiles/CyberShake_Sites.csv"
            )
            SS_psa_compss(file_station_psa, file_station_output, CyberS_sites_temp, name_Out1, name_Out2, name_Out3, name_Out5,
                   name_Out7, name_Out10)


def SS_psa_compss_nuevo(path_Stat, path_srf, CyberS_sites_temp, name_Out1, name_Out2, name_Out3, name_Out5, name_Out7,
                  name_Out10):
    db_1s = []
    db_2s = []
    db_3s = []
    db_5s = []
    db_7s = []
    db_10s = []
    # real_path_srf = path_srf
    #xj = path_srf.split('/')[:-1]
    #real_path_srf = "/" + "/".join(xj)
    for j in path_Stat:
        db_ss1_1s, db_ss1_2s, db_ss1_3s, db_ss1_5s, db_ss1_7s, db_ss1_10s = task_f_Files_compss(j, CyberS_sites_temp,
                                                                                                path_srf)#real_path_srf)
        db_1s.append(db_ss1_1s)
        db_2s.append(db_ss1_2s)
        db_3s.append(db_ss1_3s)
        db_5s.append(db_ss1_5s)
        db_7s.append(db_ss1_7s)
        db_10s.append(db_ss1_10s)
    PSA_1s = PSAs(db_1s, name_Out1)
    PSA_2s = PSAs(db_2s, name_Out2)
    PSA_3s = PSAs(db_3s, name_Out3)
    PSA_5s = PSAs(db_5s, name_Out5)
    PSA_7s = PSAs(db_7s, name_Out7)
    PSA_10s = PSAs(db_10s, name_Out10)

    return name_Out1, name_Out2, name_Out3, name_Out5, name_Out7, name_Out10

    '''PSA_1s.to_csv(name_Out1)
    PSA_2s.to_csv(name_Out2)
    PSA_3s.to_csv(name_Out3)
    PSA_5s.to_csv(name_Out5)
    PSA_7s.to_csv(name_Out7)
    PSA_10s.to_csv(name_Out10)'''


def step_3_nuevo(first_step_files, second_step_files, path_fo_file, CyberS_sites_temp, output_folder=None, save_intermediate_files=False):
    one_second_PSAs = []
    two_second_PSAs = []
    three_second_PSAs = []
    five_second_PSAs = []
    seven_second_PSAs = []
    ten_second_PSAs = []
    for file_folder_output_PSA, file_folder_output_srf in zip(first_step_files, second_step_files):
        for file_station_psa, file_station_output in zip(file_folder_output_PSA, file_folder_output_srf):
            path_out = file_station_output[0].split("/", -1)[:-2]
            i_split = file_station_output[0].split("/", -1)[-1]
            f_stat = i_split.split("_")
            stat = f_stat[1]
            if output_folder is not None:
                folder_Out = output_folder + "/step_3"
            else:
                folder_Out = "/" + "/".join(path_out) + "/step_3"
            if not os.path.exists(folder_Out):
                os.mkdir(folder_Out)
            name_Out1 = folder_Out + '/PSA_1_SS_' + stat + '.csv'
            name_Out2 = folder_Out + '/PSA_2_SS_' + stat + '.csv'
            name_Out3 = folder_Out + '/PSA_3_SS_' + stat + '.csv'
            name_Out5 = folder_Out + '/PSA_5_SS_' + stat + '.csv'
            name_Out7 = folder_Out + '/PSA_7_SS_' + stat + '.csv'
            name_Out10 = folder_Out + '/PSA_10_SS_' + stat + '.csv'
            #CyberS_sites_temp = (
            #    "./inputFiles/CyberShake_Sites.csv"
            #)
            one_second_PSA, two_second_PSA, three_second_PSA, five_second_PSA, seven_second_PSA, ten_second_PSA = SS_psa_compss_nuevo(file_station_psa, file_station_output, CyberS_sites_temp, name_Out1, name_Out2, name_Out3,
                          name_Out5,
                          name_Out7, name_Out10)
            one_second_PSAs.append(one_second_PSA)
            two_second_PSAs.append(two_second_PSA)
            three_second_PSAs.append(three_second_PSA)
            five_second_PSAs.append(five_second_PSA)
            seven_second_PSAs.append(seven_second_PSA)
            ten_second_PSAs.append(ten_second_PSA)
    return one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs


@task(psas=COLLECTION_FILE_IN, name_out=FILE_OUT, returns=1)
def merge_psas_step_4(psas, name_out):
    PSA = pd.DataFrame()
    for file in psas:
        f1 = pd.read_csv(file)
        PSA = pd.concat([PSA, f1])
    PSA.to_csv(name_out)

def step_4_compss(path_to_folders, psa_1s, psa_2s, psa_3s, psa_5s, psa_7s, psa_10s, output_folder=None):
    #path_out = psa_1s.split("/", -1)[:-3]
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

def step_4(path_to_folders, *args):
    PSA_1s = pd.DataFrame()
    PSA_2s = pd.DataFrame()
    PSA_3s = pd.DataFrame()
    PSA_5s = pd.DataFrame()
    PSA_7s = pd.DataFrame()
    PSA_10s = pd.DataFrame()
    f_Folders = []
    for (dirpath, dirnames, filenames) in walk(path_to_folders):
        f_Folders.extend(dirnames)
        break
    for nk in f_Folders:
        folder_out = path_to_folders + "/step_4"
        if not os.path.exists(folder_out):
            os.mkdir(folder_out)
        f_Stations = []
        for (dirpath, dirnames, filenames) in walk(path_to_folders + "/" +nk):
            f_Stations.extend(dirnames)
            break
        for i in f_Stations:
            path_Stat = path_to_folders + "/" +nk +"/" + i + "/step_3" #"/outputExtraction/PSA_output"
            f_Files = []
            for (dirpath, dirnames, filenames) in walk(path_Stat):
                f_Files.extend(filenames)
                break
            for j in f_Files:
                ij = j
                ij_split = ij.split('_')
                PSA_num = int(ij_split[1])
                SS = ij_split[3]
                if PSA_num == 1:
                    f1 = pd.read_csv(path_Stat + "/" + j)
                    PSA_1s = pd.concat([PSA_1s, f1])
                if PSA_num == 2:
                    f1 = pd.read_csv(path_Stat + "/" +j)
                    PSA_2s = pd.concat([PSA_2s, f1])
                if PSA_num == 3:
                    f1 = pd.read_csv(path_Stat + "/" + j)
                    PSA_3s = pd.concat([PSA_3s, f1])
                if PSA_num == 5:
                    f1 = pd.read_csv(path_Stat + "/" + j)
                    PSA_5s = pd.concat([PSA_5s, f1])
                if PSA_num == 7:
                    f1 = pd.read_csv(path_Stat + "/" + j)
                    PSA_7s = pd.concat([PSA_7s, f1])
                if PSA_num == 10:
                    f1 = pd.read_csv(path_Stat + "/" + j)
                    PSA_10s = pd.concat([PSA_10s, f1])
            PSA_1s.to_csv(folder_out + '/PSA_1.csv')
            PSA_2s.to_csv(folder_out + '/PSA_2.csv')
            PSA_3s.to_csv(folder_out + '/PSA_3.csv')
            PSA_5s.to_csv(folder_out + '/PSA_5.csv')
            PSA_7s.to_csv(folder_out + '/PSA_7.csv')
            PSA_10s.to_csv(folder_out + '/PSA_10.csv')


@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_index_compss_random(PSA, name_out):
    Site_Rup_RupVar = pd.DataFrame()
    Database = pd.read_csv(PSA)

    Site_Rup_RupVar['Source_ID'] = Database['Source_ID']
    Site_Rup_RupVar['Rupture_ID'] = Database['Rupture_ID']
    Site_Rup_RupVar['Rupture_Variation_ID'] = Database['Rupture_Variation_ID']

    Site_Rup_RupVar.to_csv(name_out, index=False)

@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_index_compss(PSA, name_out):
    Site_Rup_RupVar = pd.DataFrame()
    Database = pd.read_csv(PSA)

    Site_Rup_RupVar['Source_ID'] = Database['Source_ID']
    Site_Rup_RupVar['Rupture_ID'] = Database['Rupture_ID']
    Site_Rup_RupVar['Rupture_Variation_ID'] = Database['Rupture_Variation_ID']
    Site_Rup_RupVar['Magnitude'] = Database['Magnitude']

    Site_Rup_RupVar.to_csv(name_out, index=False)

def step_split_index(output_folder, PSA):
    name_out = output_folder + "/10_Source_Rup_RupVar_Mw.csv"#"/MwData_PSA_" + PSA.split("_", -1)[-2] + ".csv"
    split_index_compss(PSA, name_out)
    '''name_out2 = output_folder + "/Source_Rup_RupVar.csv"
    split_index_compss_random(PSA, name_out2)'''
    return name_out

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


def step_split_10Mw(path_to_folder, PSA):
    name_out = path_to_folder + "/10_source_PSA_" + PSA.split("_", -1)[-2] + "_rup_rupvar_10Mw.csv"
    split_10Mw_compss(PSA, name_out)
    return name_out

@task(PSA=FILE_IN, name_out=FILE_OUT)
def split_10Random_compss(PSA, name_out, seed=42):
    DataIn = pd.read_csv(PSA)
    dfRandom = DataIn.sample(frac=0.10, random_state=seed)

    dfRandom.to_csv(name_out, index=False)

def step_split_10Random(path_to_folder, PSA):
    name_out = path_to_folder + "/10_source_PSA_" + PSA.split("_", -1)[-2] + "_rup_rupvar.csv"
    split_10Random_compss(PSA, name_out)
    return name_out

def train_test_split_random(path_to_folder, PSA):
    pass

@task(path_to_source_rupvar=FILE_IN, psa_to_split=FILE_IN, data_test=FILE_OUT, data_train=FILE_OUT)
def train_test_split_tasks(path_to_source_rupvar, psa_to_split, data_test, data_train):
    # PASO 1 TEST
    df10Selected = pd.read_csv(path_to_source_rupvar)#pd.read_csv(path_to_folders + '/10_Source_Rup_RupVar_Mw.csv')#'/path/to/10%selected/10_Source_Rup_RupVar_Mw.csv')
    Database = pd.read_csv(psa_to_split)#'/path/to/Database.csv')

    DataOut1 = data_test#path_to_folders + '/Test.csv'
    DataOut2 = data_train#path_to_folders + '/Train.csv'

    dfEQ = pd.merge(Database, df10Selected, on=['Source_ID', 'Rupture_ID', 'Rupture_Variation_ID', 'Magnitude'],
                    how='inner')
    dfEQ.to_csv(data_test, index=False)

    # PASO 2 TRAIN
    chunk_size = 600000
    chunks = pd.read_csv(psa_to_split, chunksize=chunk_size)#'/path/to/Database.csv', chunksize=chunk_size)
    DataOut1 = pd.read_csv(data_test)#path_to_folders + '/Test.csv')

    TrainChunks = []
    for chunk in chunks:
        chunk_df = chunk
        df = pd.merge(chunk_df, DataOut1, on=['Source_ID', 'Rupture_ID', 'Rupture_Variation_ID', 'Magnitude'], how='left')
        merged = chunk_df.merge(DataOut1, how='left', indicator=True)  # left,right,both

        dfNotEQ = merged[merged['_merge'] == 'left_only']  # in the merge column, only those listed as left_only
        dfNotEQ = dfNotEQ.drop('_merge', axis=1)

        TrainChunks.append(dfNotEQ)

    Train = pd.concat(TrainChunks, ignore_index=True)
    Train.to_csv(data_train, index=False)

    return data_train, data_test#path_to_folders + '/Train.csv', data_test#path_to_folders + '/Test.csv'

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

def train_test_split(path_to_source_rupvar, path_to_folder, PSA):
    name_out_train = path_to_folder + "/Train.csv"
    name_out_test = path_to_folder + "/Test.csv"
    train_test_split_tasks(path_to_source_rupvar, PSA, name_out_test, name_out_train)
    return name_out_train, name_out_test

def generate_dislib_dataset(input_file, output_file):
    adapt_dataset_to_dislib(input_file, output_file)
    return output_file

def merged_steps_preprocess_mlesmap(path_to_folders, output_folder, table_CS_Ruptures, save_intermediate_files=False):
    CyberS_sites_temp = table_CS_Ruptures
    files_out_first_step = step_1(path_to_folders, output_folder=output_folder, save_intermediate_files=save_intermediate_files)
    compss_wait_on_file(files_out_first_step)
    files_output_second_step = step_2(path_to_folders, output_folder, save_intermediate_files=save_intermediate_files)
    compss_wait_on_file(files_output_second_step)
    one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs = step_3_nuevo(files_out_first_step, files_output_second_step, path_to_folders, CyberS_sites_temp, output_folder=output_folder, save_intermediate_files=save_intermediate_files)
    PSA_1s, PSA_2s, PSA_3s, PSA_5s, PSA_7s, PSA_10s = step_4_compss(path_to_folders, one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs, output_folder=output_folder)
    PSA_1s_sindex_Mw = step_split_index(output_folder, one_second_PSAs[-1])
    path_to_source_rupvar = step_split_10Mw(output_folder, PSA_1s_sindex_Mw)
    train, test = train_test_split(path_to_source_rupvar, output_folder, PSA_3s)
    ml_train = output_folder + "/ml_train_dataset.csv"
    ml_test = output_folder + "/ml_test_dataset.csv"
    ml_train = generate_dislib_dataset(train, ml_train)
    ml_test = generate_dislib_dataset(test, ml_test)
    compss_wait_on_file(ml_train)
    return ml_train


def main(path_to_folders, CyberS_sites_temp, output_folder=None, random_split=False, save_intermediate_files=False):
    files_out_first_step = step_1(path_to_folders, output_folder=output_folder, save_intermediate_files=save_intermediate_files)
    compss_wait_on_file(files_out_first_step)
    files_output_second_step = step_2(path_to_folders, output_folder, save_intermediate_files=save_intermediate_files)
    one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs = step_3_nuevo(files_out_first_step, files_output_second_step, path_to_folders, CyberS_sites_temp, output_folder=output_folder, save_intermediate_files=save_intermediate_files)
    PSA_1s, PSA_2s, PSA_3s, PSA_5s, PSA_7s, PSA_10s = step_4_compss(path_to_folders, one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs, output_folder=output_folder)
    PSA_1s_sindex_Mw = step_split_index(output_folder, one_second_PSAs[-1])
    path_to_source_rupvar = step_split_10Mw(output_folder, PSA_1s_sindex_Mw)
    train, test = train_test_split(path_to_source_rupvar, output_folder, PSA_3s)
    ml_train = output_folder + "/ml_train_dataset.csv"
    ml_test = output_folder + "/ml_test_dataset.csv"
    ml_train = generate_dislib_dataset(train, ml_train)
    ml_test = generate_dislib_dataset(test, ml_test)
    compss_wait_on_file(ml_train)
    return ml_train

