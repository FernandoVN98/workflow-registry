import os
from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *
from pycompss.api.multinode import multinode


@constraint(computing_units="${ComputingUnits}")
@task(cordfile=FILE_IN, in3d_file=FILE_IN, returns=6)
def read_previous_files(cordfile, in3d_file, dependency):
    fp_in = open(in3d_file, "r")
    data = fp_in.readlines()
    fp_in.close()
    fp_in = open(cordfile, "r")
    fp_in.readline()
    fp_in.readline()
    fp_in.readline()
    fp_in.readline()
    num_sgt_pts = int(fp_in.readline())
    fp_in.close()
    params = dict()
    for line in data:
        try:
            [value, key] = line.split()
            params[key] = value
        except Exception as e:
            continue

    total_ts = int(float(params["TMAX"])/float(params["DT"])+0.5)
    decimation = int(params["NTISKP_SGT"])
    spacing = float(params["DH"])/1000.0
    params_dt = float(params["DT"])
    total_ts_decimation = total_ts/decimation
    return total_ts, decimation, spacing, total_ts_decimation, num_sgt_pts, params_dt


@constraint(computing_units="${ComputingUnitspostAWP}")
@multinode(computing_nodes="${ComputingNodespostAWP}")
@task(input_filename=FILE_IN, returns=1)#output_file=FILE_OUT, returns=1)
def reformat(cybershake_path, input_filename, timesteps, num_pts, output_filename, comp):
    #On Titan, need aprun to execute this on a compute node, not the aprun node
    command = "srun %s/SgtHead/bin/reformat_awp_mpi %s %d %d %s" % (cybershake_path, input_filename, timesteps, num_pts, output_filename)
    if comp=="z":
        #Add z flag to apply additional factor of two
        command = "%s -z" % command
    exitcode = os.system(command)
    #fp_in = open(output_filename, "rb")
    #fp_out = open(output_file, "wb")
    #data = fp_in.read()
    #fp_out.write(data)
    #fp_out.close()
    #fp_in.close()
    return exitcode



@constraint(computing_units="${ComputingUnitspostAWP}")
@multinode(computing_nodes="${ComputingNodespostAWP}")
@task(cordfile=FILE_IN, fdloc=FILE_IN, gridout=FILE_IN, modelbox=FILE_IN, returns=1)# header_file_out=FILE_OUT, returns=1)
def write_head(cybershake_path, modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, source_freq, media, header_name, header_file_out, dep):
    #On Titan, need aprun to execute this on a compute node, not the aprun node
    command = "srun %s/SgtHead/bin/write_head %s %s %s %s %f %d %f %d %s %s %f %s %s" % (cybershake_path, modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, source_freq, media, header_name)
    if not comp=="z":
        #Add c flag to use corrected mu
        command = "%s -c " % command
    exitcode = os.system(command)
    return exitcode

def step_postAWP(cybershake_path, run_path, variable, site, model_box, model_gridout, runID, frequency, source_frequency, gravesPitarka, erf_id, database_path, cordfile, fdloc, media, in3d_file, dependency, skip_md5=False):
    awp_sgt_filename = run_path + "/comp_"+variable+"/output_sgt/awp-strain-" + site+ "-f"+variable
    rwg_comp = "z"
    if variable=="x":
        rwg_comp = "y"
    elif variable=="y":
        rwg_comp = "x"
    awp_reformat_sgt_filename = run_path + "/" + site + "_f" + rwg_comp + "_" + str(runID) + ".sgt"
    awp_reformat_sgt_file = run_path + "/" + site + "_f" + rwg_comp + "_" + str(runID) + ".sgt"
    run_id = str(runID)
    MOMENT = "1.0e20"
    header_out = run_path + "/" + site + "_f" + variable + "_" + str(runID) + ".sgthead"
    total_ts, decimation, spacing, total_ts_decimation, num_sgt_pts, params_dt = read_previous_files(cordfile, in3d_file, dependency)
    rc = reformat(cybershake_path, awp_sgt_filename, total_ts_decimation, num_sgt_pts, awp_reformat_sgt_filename, variable)
    rwg_comp = "z"
    if variable=="x":
        rwg_comp = "y"
    elif variable=="y":
        rwg_comp = "x"
    header_name = run_path + "/%s_f%s_%s.sgthead" % (site, rwg_comp, run_id)
    header_file_out = run_path + "/%s_f%s_%s.sgthead" % (site, rwg_comp, run_id)
    rc = write_head(cybershake_path, model_box, cordfile, fdloc, model_gridout, spacing, total_ts, params_dt, decimation, rwg_comp, MOMENT, source_frequency, media, header_name, header_file_out, rc)
    if not skip_md5:
        os.system("md5sum %s > %s.md5" % (run_path + "/" + site + "_f" + variable + "_" + str(runID) + ".sgt", run_path + "/" + site + "_f" + variable + "_" + str(runID) + ".sgt"))
    return rc, awp_reformat_sgt_filename, header_file_out

