from pycompss.api.task import task
from pycompss.api.mpi import mpi
from pycompss.api.parameter import *


@mpi(runner="srun", binary="$PMCL3D_BIN", args="{{input_file}}", processes="$CYBERSHAKE_PROCS", processes_per_node="$CYBERSHAKE_PPN", working_dir="work_dir}}")
@task(input_file=FILE_IN, auxiliar_files=COLLECTION_FILE_IN, returns=1)
def AWP_variable(input_file, auxiliar_files, work_dir):
    pass
