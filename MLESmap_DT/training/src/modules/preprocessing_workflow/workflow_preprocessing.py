from modules.preprocessing_workflow.steps import step_srf_extraction, step_srf_extraction_future_object, step_seismogram_extraction_future_object, step_seismogram_extraction, step_srf_database, step_merge_database, step_split_index, step_split_10Mw, train_test_split, generate_dislib_dataset
from pycompss.api.api import compss_wait_on_file


def obtain_files_from_directory(directory_path):
    seismogram_files = []
    srf_files = []
    for (dirpath, dirnames, filenames) in walk(directory_path):
        for filename in filenames:
            if filename.endswith(".grm"):
                seismogram_files.append(dirpath + "/" + filename)
            elif filename.endswith(".srf"):
                srf_files.append(dirpath + "/" + filename)
    return seismogram_files, srf_files



def load_seismogram_files_directories(directories):
    numbers_station = []
    seismogram_files_out = []
    srf_files_out = []
    for nk in directories:
        nk_number = nk.split("_")[-2]
        try:
            index_file_out = numbers_station.index(nk_number)
            seismogram_files, srf_files = obtain_files_from_directory(nk)
            seismogram_files_out[index_file_out].append([seismogram_files])
            srf_files_out[index_file_out].append([srf_files])
        except ValueError:
            seismogram_files, srf_files = obtain_files_from_directory(nk)
            seismogram_files_out.append([seismogram_files])
            srf_files_out.append([srf_files])
            numbers_station.append(nk_number)
    return seismogram_files_out, srf_files_out



def merged_steps_preprocess_mlesmap(output_directories, path_to_folders, output_folder, table_CS_Ruptures, random_split=False, intensity_split=True):
    if len(output_directories) == 1 and isinstance(output_directories, list):
        output_directories = output_directories[0]
    for directory_out in output_directories:
        compss_wait_on_directory(directory_out)
    seismogram_files, srf_files = load_seismogram_files_directories(output_directories) 
    CyberS_sites_temp = table_CS_Ruptures
    files_out_first_step = step_seismogram_extraction_future_object(path_to_folders, seismogram_files, output_folder)
    files_output_second_step = step_srf_extraction_future_object(path_to_folders, srf_files, output_folder=output_folder)
    one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs = step_srf_database(files_out_first_step, files_output_second_step, path_to_folders, CyberS_sites_temp, output_folder=output_folder)
    PSA_1s, PSA_2s, PSA_3s, PSA_5s, PSA_7s, PSA_10s = step_merge_database(path_to_folders, one_second_PSAs, two_second_PSAs, three_second_PSAs, five_second_PSAs, seven_second_PSAs, ten_second_PSAs, output_folder=output_folder)
    PSA_1s_sindex_Mw = step_split_index(output_folder, one_second_PSAs[-1])
    if intensity_split:
        path_to_source_rupvar = step_split_10Mw(output_folder, PSA_1s_sindex_Mw)
        train, test = train_test_split(path_to_source_rupvar, output_folder, PSA_3s,
                "/Train_10MW.csv", "/Test_10MW.csv")
        ml_train = output_folder + "/ml_train_dataset_MW.csv"
        ml_test = output_folder + "/ml_test_dataset_MW.csv"
        ml_train = generate_dislib_dataset(train, ml_train)
        ml_test = generate_dislib_dataset(test, ml_test)
    if random_split:
        path_to_source_rupvar = step_split_10Mw(output_folder, PSA_1s_sindex_Mw, random=True)
        train_rand, test_rand = train_test_split(path_to_source_rupvar, output_folder, PSA_3s, "/Train_rand.csv",
                "/Test_rand.csv")
        ml_train_rand = output_folder + "/ml_train_dataset_random.csv"
        ml_test_rand = output_folder + "/ml_test_dataset_random.csv"
        ml_train_rand = generate_dislib_dataset(train_rand, ml_train_rand)
        ml_test_rand = generate_dislib_dataset(test_rand, ml_test_rand)
    if intensity_split and random_split:
        compss_wait_on_file(ml_train)
        compss_wait_on_file(ml_train_rand)
        return ml_train, ml_train_rand
    elif intensity_split:
        compss_wait_on_file(ml_train)
        return ml_train
    elif random_split:
        compss_wait_on_file(ml_train_rand)
        return ml_train_rand
