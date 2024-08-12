import os
import numpy as np
from os import walk
from modules.preprocessing_workflow.read_seis_TEST import *
from modules.preprocessing_workflow.rotinv_maxavg import *
import pandas as pd
from pycompss.api.constraint import constraint
from pycompss.api.task import task
from pycompss.api.api import compss_wait_on
from pycompss.api.parameter import *


def step_seismogram_extraction(path_to_folders, output_folder=None, save_intermediate_files=False):
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

@constraint(computing_units="${ComputingUnits}")
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
