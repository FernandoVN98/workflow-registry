from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *
import os
import pandas as pd
import numpy as np
import math
from geopy import Point
from geopy.distance import geodesic


def step_srf_database(first_step_files, second_step_files, path_fo_file, CyberS_sites_temp, output_folder=None, save_intermediate_files=False):
    one_second_PSAs = []
    two_second_PSAs = []
    three_second_PSAs = []
    five_second_PSAs = []
    seven_second_PSAs = []
    ten_second_PSAs = []
    files_station_psas = []
    files_station_outputs = []
    print("FIRST STEP FILES")
    print(first_step_files)
    print("SECOND STEP FILES")
    print(second_step_files)
    for file_folder_output_PSA, file_folder_output_srf in zip(first_step_files, second_step_files):
        for file_station_psa, file_station_output in zip(file_folder_output_PSA, file_folder_output_srf):
            path_out = file_station_output[0].split("/", -1)[:-2]
            i_split = file_station_output[0].split("/", -1)[-1]
            f_stat = i_split.split("_")
            stat = f_stat[2]
            if output_folder is not None:
                folder_Out = output_folder + "/step_3/Station_" + str(f_stat[2])
            else:
                folder_Out = "/" + "/".join(path_out) + "/step_3/Station_" + str(f_stat[2])
            if not os.path.exists(folder_Out):
                if not os.path.exists(output_folder + "/step_3"):
                    os.mkdir(output_folder + "/step_3")
                os.mkdir(folder_Out)
            one_second_PSAs.append(folder_Out + '/PSA_1_SS_' + stat + '.csv')
            two_second_PSAs.append(folder_Out + '/PSA_2_SS_' + stat + '.csv')
            three_second_PSAs.append(folder_Out + '/PSA_3_SS_' + stat + '.csv')
            five_second_PSAs.append(folder_Out + '/PSA_5_SS_' + stat + '.csv')
            seven_second_PSAs.append(folder_Out + '/PSA_7_SS_' + stat + '.csv')
            ten_second_PSAs.append(folder_Out + '/PSA_10_SS_' + stat + '.csv')
            files_station_psas.append(file_station_psa)
            files_station_outputs.append(file_station_output)
    for i in range(len(files_station_psas)):
        SS_psa_compss_nuevo(files_station_psas[i], files_station_outputs[i],
                CyberS_sites_temp,
                one_second_PSAs[i], two_second_PSAs[i], three_second_PSAs[i],
                five_second_PSAs[i], seven_second_PSAs[i], ten_second_PSAs[i])
    return one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs



def step_srf_database_orig(first_step_files, second_step_files, path_fo_file, CyberS_sites_temp, output_folder=None, save_intermediate_files=False):
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
                folder_Out = output_folder + "/step_3/Station_" + str(f_stat[2])
            else:
                folder_Out = "/" + "/".join(path_out) + "/step_3/Station_" + str(f_stat[2])
            if not os.path.exists(folder_Out):
                if not os.path.exists(output_folder + "/step_3"):
                    os.mkdir(output_folder + "/step_3")
                os.mkdir(folder_Out)
            name_Out1 = folder_Out + '/PSA_1_SS_' + stat + '.csv'
            name_Out2 = folder_Out + '/PSA_2_SS_' + stat + '.csv'
            name_Out3 = folder_Out + '/PSA_3_SS_' + stat + '.csv'
            name_Out5 = folder_Out + '/PSA_5_SS_' + stat + '.csv'
            name_Out7 = folder_Out + '/PSA_7_SS_' + stat + '.csv'
            name_Out10 = folder_Out + '/PSA_10_SS_' + stat + '.csv'
            one_second_PSA, two_second_PSA, three_second_PSA, five_second_PSA, \
            seven_second_PSA, ten_second_PSA = SS_psa_compss_nuevo(
                file_station_psa, file_station_output, CyberS_sites_temp, name_Out1,
                name_Out2, name_Out3, name_Out5,
                name_Out7, name_Out10)
            one_second_PSAs.append(one_second_PSA)
            two_second_PSAs.append(two_second_PSA)
            three_second_PSAs.append(three_second_PSA)
            five_second_PSAs.append(five_second_PSA)
            seven_second_PSAs.append(seven_second_PSA)
            ten_second_PSAs.append(ten_second_PSA)
    return one_second_PSAs, two_second_PSAs, three_second_PSAs, \
           five_second_PSAs, seven_second_PSAs, ten_second_PSAs

def SS_psa_compss_nuevo(path_Stat, path_srf, CyberS_sites_temp, name_Out1,
                        name_Out2, name_Out3, name_Out5, name_Out7,
                        name_Out10):
    db_1s = []
    db_2s = []
    db_3s = []
    db_5s = []
    db_7s = []
    db_10s = []
    for j in path_Stat:
        print("FILE SENT TO TASK F FILES")
        print(j)
        db_ss1_1s, db_ss1_2s, db_ss1_3s, db_ss1_5s, \
        db_ss1_7s, db_ss1_10s = task_f_Files_compss(j, CyberS_sites_temp,
                                                    path_srf)
        db_1s.append(db_ss1_1s)
        db_2s.append(db_ss1_2s)
        db_3s.append(db_ss1_3s)
        db_5s.append(db_ss1_5s)
        db_7s.append(db_ss1_7s)
        db_10s.append(db_ss1_10s)
    PSAs(db_1s, name_Out1)
    PSAs(db_2s, name_Out2)
    PSAs(db_3s, name_Out3)
    PSAs(db_5s, name_Out5)
    PSAs(db_7s, name_Out7)
    PSAs(db_10s, name_Out10)

    return name_Out1, name_Out2, name_Out3, name_Out5, name_Out7, name_Out10


@constraint(computing_units="${ComputingUnits}")
@task(file=FILE_IN, CyberS_sites_temp=FILE_IN, path_srf=COLLECTION_FILE_IN, returns=6)
def task_f_Files_compss(file, CyberS_sites_temp, path_srf):
    ij = file.split('/')[-1]
    ij_split = ij.split('_')
    print("IJ SPLIT")
    print(ij_split, flush=True)
    SS_csv = int(ij_split[4])
    src_csv = int(ij_split[6])
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


@constraint(computing_units="${ComputingUnits}")
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
