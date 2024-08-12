from pycompss.api.task import task
from pycompss.api.mpi import mpi
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *
import os
from modules.cybershake_workflows.steps.utils import getProperty
from pycompss.api.multinode import multinode
from modules.cybershake_workflows.utils.DAL import SQLiteHandler
import subprocess
import sys


@constraint(computing_units="${ComputingUnits}")
@task(rupture_file=FILE_IN, rupture_file_list=FILE_IN, returns=1)
def make_lns_orig(rupture_file, ruptures_erf, rupture_file_list, run_path):
    rupture_root = getProperty("RUPTURE_ROOT", rupture_file)
    rupture_root = rupture_root + ruptures_erf
    if len(sys.argv)<3:
        print("Usage: %s <rupture file list> <rupture root directory>" % sys.argv[0])
        sys.exit(1)
    with open(rupture_file_list, "r") as fp_in:
        data = fp_in.readlines()
        fp_in.close()
        for line in data[1:]:
            lfn = line.split()[0]
            #e36_rv6_39_5.txt
            pieces = lfn.split("_")
            src = int(pieces[2])
            rup = int(pieces[3].split(".")[0])
            try:
                print("COMMMAAAAND", flush=True)
                print("cp %s/%d/%d/%d_%d.txt " % (rupture_root, src, rup, src, rup) + run_path + "/post-processing/" + lfn, flush=True)
                os.system("cp %s/%d/%d/%d_%d.txt " % (rupture_root, src, rup, src, rup) + run_path + "/post-processing/" + lfn)
                #os.symlink("%s/%d/%d/%d_%d.txt" % (rupture_root, src, rup, src, rup), run_path + "/post-processing/" + lfn)
            except OSError as e:
                print("Exception")
                print(e)
    return 1


@constraint(computing_units="${ComputingUnits}")
@task(rupture_file=FILE_IN, rupture_file_list=FILE_IN, working_directory=DIRECTORY_INOUT)
def make_lns(rupture_file, ruptures_erf, rupture_file_list, run_path, working_directory):
    rupture_root = getProperty("RUPTURE_ROOT", rupture_file)
    rupture_root = rupture_root + ruptures_erf
    if len(sys.argv)<3:
        print("Usage: %s <rupture file list> <rupture root directory>" % sys.argv[0])
        sys.exit(1)
    with open(rupture_file_list, "r") as fp_in:
        data = fp_in.readlines()
        fp_in.close()
        for line in data[1:]:
            lfn = line.split()[0]
            #e36_rv6_39_5.txt
            pieces = lfn.split("_")
            src = int(pieces[2])
            rup = int(pieces[3].split(".")[0])
            try:
                os.symlink("%s/%d/%d/%d_%d.txt" % (rupture_root, src, rup, src, rup), working_directory + "/" + lfn)
            except OSError as e:
                print("Exception")
                print(e)


@constraint(computing_units="${DirectSint_CU}")
@multinode(computing_nodes="${DirectSint_NODES}", processes_per_node="${DirectSint_PPN}")
@task(rupt_file=FILE_IN, work_dir=DIRECTORY_INOUT, returns=1)
def direct_synt(site, cybershake_path, database_path, run_path, graves_pitarka, sgt_handlers, runID, debug, max_buf, ntout, rupt_file, x_sgt_file, y_sgt_file, x_header, y_header, setup_frequency, stoch_freq, psa, rotd, dur, dtouts, out_pointsX, out_pointsY, timesamples, timeskip, highHZ, dep_pre, work_dir):
    os.environ["SLURM_CPUS_ON_NODE"]="112"
    os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x" + str(os.environ["DirectSint_NODES"]) + ")"
    os.environ["SLURM_NPROCS"]=str(int(os.environ["DirectSint_PPN"])*int(os.environ["DirectSint_NODES"]))#"168"
    os.environ["SLURM_NTASKS"]=str(int(os.environ["DirectSint_PPN"])*int(os.environ["DirectSint_NODES"]))#"168"
    os.environ["SLURM_TASKS_PER_NODE"]=str(os.environ["DirectSint_PPN"]) + "x(" + str(os.environ["DirectSint_NODES"]) + ")"
    os.environ["SLURM_THREADS_PER_CORE"]="1"
    os.environ["SLURM_CPUS_PER_TASK"]="2"
    os.environ["SLURM_MEM_PER_CPU"]="2000"
    os.environ["SLURM_NTASKS_PER_NODE"]=str(os.environ["DirectSint_PPN"])#"56"
    os.environ["COMPSS_WORKING_DIR"]= work_dir
    os.chdir(work_dir)
    file_sh = work_dir + "/run_DS_.sh"
    file_python = work_dir + "/script.py"
    f_in = open(file_sh, "w")
    f_in.write("#!/bin/bash\n\n")
    f_in.write("module load python/3.12.1\n")
    f_in.write("module purge\n")
    f_in.write("module load impi\n")
    f_in.write("module load gcc\n")
    f_in.write("module load intel\n")
    f_in.write("module load fftw\n\n")
    f_in.write("ulimit -c unlimited\n\n")
    f_in.write("export LD_LIBRARY_PATH=/gpfs/scratch/bsc19/bsc019756/CyberShake/software/AWP-ODC-SGT/bin:$LD_LIBRARY_PATH\n\n")
    f_in.write("cd " + run_path + "/post-processing/\n\n")
    f_in.write("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(graves_pitarka) + "\n")
    f_in.write("%s\n" % file_python)
    f_in.flush()
    f_in.close()
    dal = SQLiteHandler(database_path, True)
    lat, lon = dal.getSiteLocation(site)
    f_python = open(file_python, "w")
    f_python.write("#!/usr/bin/env python\n\n")
    f_python.write("import sys\n\n")
    f_python.write("import os\n\n")
    f_python.write("ds_path = \"" + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(graves_pitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + work_dir + "/" + str(x_sgt_file) + " sgt_yfile=" + work_dir + "/" + str(y_sgt_file) + " x_header=" + work_dir + "/" + str(x_header) + " y_header=" + work_dir + "/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filter_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_byteswap=no\"\n\n")
    f_python.write("cmd = \"srun \" + str(ds_path)\n\n")
    f_python.write("os.system(cmd)\n")
    f_python.flush()
    f_python.close()
    os.system("chmod 755 " + work_dir + "/run_DS_.sh")
    os.system("chmod 755 " + work_dir + "/script.py")
    os.system("module load python/3.12.1")
    os.system("module purge")
    os.system("module load impi gcc intel fftw")
    cmd = "cp " + run_path + "/" + site + "_fx_" + str(runID) + ".sgt " +work_dir + "/"+ site+ "_fx_" + str(runID) + ".sgt"
    os.system(cmd)
    cmd = "cp " + run_path + "/" + site + "_fy_" + str(runID) + ".sgt " +work_dir+ "/"+ site+ "_fy_" + str(runID) + ".sgt"
    os.system(cmd)
    cmd = "cp " + run_path + "/" + site + "_fx_" + str(runID) + ".sgthead " + work_dir+ "/"+ site + "_fx_" + str(runID) + ".sgthead"
    os.system(cmd)
    cmd = "cp " + run_path + "/" +site + "_fy_" + str(runID) + ".sgthead " + work_dir + "/"+ site + "_fy_" + str(runID) + ".sgthead"
    os.system(cmd)

    command = [work_dir + "/run_DS_.sh"]
    process = subprocess.Popen(command, env=os.environ, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Obtén la salida y los errores
    stdout, stderr = process.communicate()

    # Imprime la salida
    print("Salida:", flush=True)
    print(stdout.decode(), flush=True)

    # Imprime los errores
    print("Errores:")
    print(stderr.decode())
    return stderr.decode()


@constraint(computing_units="${DirectSint_CU}")
@multinode(computing_nodes="${DirectSint_NODES}", processes_per_node="${DirectSint_PPN}")
@task(rupt_file=FILE_IN,  returns=1)
def direct_synt_orig(site, cybershake_path, database_path, run_path, graves_pitarka, sgt_handlers, runID, debug, max_buf, ntout, rupt_file, x_sgt_file, y_sgt_file, x_header, y_header, setup_frequency, stoch_freq, psa, rotd, dur, dtouts, out_pointsX, out_pointsY, timesamples, timeskip, highHZ, dependency):
    os.environ["SLURM_CPUS_ON_NODE"]="112"
    os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x" + str(os.environ["DirectSint_NODES"]) + ")"
    os.environ["SLURM_NPROCS"]=str(int(os.environ["DirectSint_PPN"])*int(os.environ["DirectSint_NODES"]))#"168"
    os.environ["SLURM_NTASKS"]=str(int(os.environ["DirectSint_PPN"])*int(os.environ["DirectSint_NODES"]))#"168"
    os.environ["SLURM_TASKS_PER_NODE"]=str(os.environ["DirectSint_PPN"]) + "x(" + str(os.environ["DirectSint_NODES"]) + ")"
    os.environ["SLURM_THREADS_PER_CORE"]="1"
    os.environ["SLURM_CPUS_PER_TASK"]="2"
    os.environ["SLURM_MEM_PER_CPU"]="2000"
    os.environ["SLURM_NTASKS_PER_NODE"]=str(os.environ["DirectSint_PPN"])#"56"
    os.environ["COMPSS_WORKING_DIR"]= run_path + "/post-processing/"
    os.chdir(run_path + "/post-processing/")
    file_sh = run_path + "/post-processing/run_DS_.sh"
    file_python = run_path + "/post-processing/script.py"
    f_in = open(file_sh, "w")
    f_in.write("#!/bin/bash\n\n")
    f_in.write("module load python/3.12.1\n")
    f_in.write("module purge\n")
    f_in.write("module load impi\n")
    f_in.write("module load gcc\n")
    f_in.write("module load intel\n")
    f_in.write("module load fftw\n\n")
    f_in.write("ulimit -c unlimited\n\n")
    f_in.write("export LD_LIBRARY_PATH=/gpfs/scratch/bsc19/bsc019756/CyberShake/software/AWP-ODC-SGT/bin:$LD_LIBRARY_PATH\n\n")
    f_in.write("cd " + run_path + "/post-processing/\n\n")
    f_in.write("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(graves_pitarka) + "\n")
    f_in.write("%s\n" % file_python)
    f_in.flush()
    f_in.close()
    dal = SQLiteHandler(database_path, True)
    lat, lon = dal.getSiteLocation(site)
    f_python = open(file_python, "w")
    f_python.write("#!/usr/bin/env python\n\n")
    f_python.write("import sys\n\n")
    f_python.write("import os\n\n")
    f_python.write("ds_path = \"" + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(graves_pitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + run_path + "/post-processing/" + str(x_sgt_file) + " sgt_yfile=" + run_path + "/post-processing/" + str(y_sgt_file) + " x_header=" + run_path + "/post-processing/" + str(x_header) + " y_header=" + run_path + "/post-processing/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filter_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_byteswap=no\"\n\n")
    f_python.write("cmd = \"srun \" + str(ds_path)\n\n")
    f_python.write("os.system(cmd)\n")
    f_python.flush()
    f_python.close()
    os.system("chmod 755 " + run_path + "/post-processing/run_DS_.sh")
    os.system("chmod 755 " + run_path + "/post-processing/script.py")
    os.system("module load python/3.12.1")
    os.system("module purge")
    os.system("module load impi gcc intel fftw")
    cmd = "ln -s " + run_path + "/" + site + "_fx_" + str(runID) + ".sgt " +run_path + "/"  + "/post-processing/"+ site+ "_fx_" + str(runID) + ".sgt"
    os.system(cmd)
    cmd = "ln -s " + run_path + "/" + site + "_fy_" + str(runID) + ".sgt " +run_path + "/"  + "/post-processing/"+ site+ "_fy_" + str(runID) + ".sgt"
    os.system(cmd)
    cmd = "ln -s " + run_path + "/" + site + "_fx_" + str(runID) + ".sgthead " + run_path + "/post-processing/"+ site + "_fx_" + str(runID) + ".sgthead"
    os.system(cmd)
    cmd = "ln -s " + run_path + "/" +site + "_fy_" + str(runID) + ".sgthead " + run_path + "/post-processing/"+ site + "_fy_" + str(runID) + ".sgthead"
    os.system(cmd)

    command = [run_path + "/post-processing/run_DS_.sh"]
    process = subprocess.Popen(command, env=os.environ, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Obtén la salida y los errores
    stdout, stderr = process.communicate()

    # Imprime la salida
    print("Salida:", flush=True)
    print(stdout.decode(), flush=True)

    # Imprime los errores
    print("Errores:")
    print(stderr.decode())
    return stderr.decode()


def step_runDS(cybershake_path, database_path, rupture_file, rupture_file_list, run_path, site, runID, input_region, input_erf_ruptures, input_erf_id, graves_pitarka, compute_setup_frequency, config):
    working_directory_out = run_path + "/post-processing"
    dependency = make_lns(rupture_file, "Ruptures_erf" + str(input_erf_id) + "/", rupture_file_list, run_path, working_directory_out)
    direct_synt(site, cybershake_path, database_path, run_path, graves_pitarka, os.environ["num_sgt_handlers"], runID, 1, 512, 3000, rupture_file_list, site + "_fx_" + str(runID) + ".sgt", site + "_fy_" + str(runID) + ".sgt", site + "_fx_" + str(runID) + ".sgthead", site + "_fy_" + str(runID) +".sgthead", config["compute"]["setup"]["frequency"], -1.0, 1, 1, 1, 0.1, 2, 1, 3000, 0.1, 5.0, dependency, working_directory_out)
    return working_directory_out 

