from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *
from pycompss.api.on_failure import on_failure
from modules.cybershake_workflows.steps.utils import getProperty
from modules.cybershake_workflows.utils.DAL import SQLiteHandler
import os
import sqlite3
from pycompss.api.multinode import multinode
import subprocess


@task(rupture_out=FILE_OUT)
def generate_new_rupture_out(rupture_out, database_path, dependency):
    dal = SQLiteHandler(database_path, True)#ESTO ESTA BIEN???
    dal.generateRuptureFile(rupture_out)#ESTO ESTA BIEN???



@on_failure(management='IGNORE', returns=0)
@constraint(computing_units="${ComputingUnitsrupVar}")
@multinode(computing_nodes="${ComputingNodesrupVar}")
@task(rupture_file=FILE_IN, returns=1)
def populate_rvs(cybershake_path, graves_pitarka, erf_id, database_path, run_path, input_region, rupture_file, dependency, dependency2):
    conn = sqlite3.connect(database_path)
    cur = conn.cursor()
    RUPTURE_ROOT = getProperty("RUPTURE_ROOT", rupture_file)
    query_prefix = "insert into Rupture_Variations (Source_ID, Rupture_ID, Rup_Var_ID, ERF_ID, Rup_Var_Scenario_ID, Rup_Var_LFN) values "
    rup_var_scenario_id = 1 #ESTO ESTABA ASíANTERIORMENTE, NO SE SI ESTA BIEN
    root = "%s/Ruptures_erf%d" % (RUPTURE_ROOT, erf_id)
    for src_dir in os.listdir(root):
        src_entry = os.path.join(root, src_dir)
        if os.path.isdir(src_entry):
            src = int(src_dir)
            print("Processing directory %s." % src_entry)
            for rup_dir in os.listdir(src_entry):
                rup_entry = os.path.join(src_entry, rup_dir)
                if os.path.isdir(rup_entry):
                    rup = int(rup_dir)
                    print("Processing directory %s." % rup_entry)
                    for f in os.listdir(rup_entry):
                        if f.find("%d_%d.txt" % (src, rup))==0:
                            cmd =  cybershake_path + "/RuptureCodes/RupGen-api-%s/utils/num_rvs %s" % (str(graves_pitarka), os.path.join(rup_entry, f))
                            result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
                            num_rvs = int(result.stdout.split()[-2])
                            for i in range(0, num_rvs):
                                query = "%s (%d, %d, %d, %d, %d, '%s')" % (query_prefix, src, rup, i, erf_id, rup_var_scenario_id, "e%d_rv%d_%d_%d.txt.variation-r%06d" % (erf_id, rup_var_scenario_id, src, rup, i))
                                print("QUERYYYYY", flush=True)
                                print(query, flush=True)
                                cur.execute(query)
    conn.commit()
    conn.close()
    return 1

def rupVar(cybershake_path, graves_pitarka, erf_id, database_path, run_path, input_region, rupture_file, dependency, dependency2):
    cmd = "mkdir " + run_path +"/post-processing"
    print("CMDDDDD")
    print(cmd, flush=True)
    os.system(cmd)
    rupture_out = run_path + "/post-processing/rupture_file_list_" + str(input_region)
    dependency_pop_rvs = populate_rvs(cybershake_path, graves_pitarka, erf_id, database_path, run_path, input_region, rupture_file, dependency, dependency2)
    generate_new_rupture_out(rupture_out, database_path, dependency_pop_rvs)
    return rupture_out

