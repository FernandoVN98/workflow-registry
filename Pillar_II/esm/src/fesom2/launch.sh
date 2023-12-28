#!/bin/bash

# Code used:
# - https://stackoverflow.com/a/29754866 CC-BY-SA 4.0
# - https://stackoverflow.com/a/246128 CC-BY-SA 4.0
set -o errexit -o pipefail -o noclobber -o nounset

# This script is used to launch experiments with FESOM2.

MODEL="fesom2"

HPC=""
DEBUG=""
CORES=0
CORES_PER_NODE=48
MEMBERS=1
QOS=debug

# Parse options. Note that options may be followed by one colon to indicate
# they have a required argument.
if ! options=$(getopt --name "$(basename "$0")" --options dc:m:q: --longoptions debug,cores:,members:,qos:,hpc:,cores_per_node: -- "$@"); then
  # Error, getopt will put out a message for us
  exit 1
fi

eval set -- "${options}"

while [ $# -gt 0 ]; do
  # Consume next (1st) argument
  case "$1" in
  -d | --debug)
    DEBUG="--debug"
    ;;
  # Options with required arguments, an additional shift is required
  -c | --cores)
    CORES="$2"
    shift
    ;;
  -m | --members)
    MEMBERS="$2"
    shift
    ;;
  -q | --qos)
    QOS="$2"
    shift
    ;;
  --hpc)
    HPC="${2,,}"
    shift
    ;;
  --cores_per_node)
    CORES_PER_NODE="${2}"
    shift
    ;;
  --)
    shift
    break
    ;;
  -*)
    echo "$0: error - unrecognized option $1" 1>&2
    exit 1
    ;;
  *)
    break
    ;;
  esac
  # Fetch next argument as 1st
  shift
done

if [ -z "${HPC}" ]; then
  echo "Please provide a valid HPC environment name"
  exit 1
fi

if [ "${CORES}" -le 0 ]; then
  echo "Cores must be equal or greater than 1"
  exit 1
fi

if [ "${CORES_PER_NODE}" -le 0 ]; then
  echo "Cores per node must be equal or greater than 1"
  exit 1
fi

if [ "${MEMBERS}" -le 0 ]; then
  echo "Members must be equal or greater than 1"
  exit 1
fi

printf "Launching %s eFlows4HPC ESM experiment...\U1F680\n" "${MODEL}"

# Hecuba configuration
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
HECUBA_CONFIGURATION="$(realpath -e -- "${SCRIPT_DIR}/storage_props.cfg")"

# Experiment configuration. The variables exported here are used
# by PyCOMPSs (some Python decorators use values like ="${FESOM_CORES}").
FESOM_CORES=${CORES}

# math.ceil(fesom_cores / cores_per_node)
NODE_ALLOCATION="$(((FESOM_CORES + CORES_PER_NODE - 1) / CORES_PER_NODE))"
# Now multiply by number of ensemble members...
NODE_ALLOCATION="$((NODE_ALLOCATION * MEMBERS))"

if [ "${CORES}" -lt "${CORES_PER_NODE}" ]; then
  echo "WARNING: You have ${CORES_PER_NODE} cores per node, but requested less: ${CORES}"
fi
if [ $((CORES % CORES_PER_NODE)) -ne 0 ]; then
  echo "WARNING: You are not using all the cores of your nodes (${CORES_PER_NODE}), you requested: ${CORES}"
fi

export FESOM_CORES
export QOS
export MEMBERS
export NODE_ALLOCATION

echo -e "\nLaunch arguments:\n"

echo "MODEL           : ${MODEL}"
echo "HPC             : ${HPC}"
echo "CORES           : ${FESOM_CORES}"
echo "CORES PER NODE  : ${CORES_PER_NODE}"
echo "NODES           : ${NODE_ALLOCATION}"
echo "QOS             : ${QOS}"
echo "MEMBERS         : ${MEMBERS}"

echo -e "\nLoading ${HPC} configurations..."
HPC_ENV_FILE="$(realpath -e -- "${SCRIPT_DIR}/env/${HPC}.sh")"
# Use an example so shellcheck can at least check that one when parsing
# this file (you can lint all files independently from this).
# shellcheck source=env/mn4.sh
source "${HPC_ENV_FILE}"
echo -e "Done! ${HPC} environment loaded correctly!\n"

# Sample invocation of this script:
#
# ./launch.sh --hpc mn4 --cores 288 --qos debug --members 1

# EXP_ID=$(printf "%06d\n" $((1 + $RANDOM % 100000)))
# --expid is now optional. Python does the same thing.

# Launch the ESM ensemble simulation with Hecuba infrastructure using COMPSs.
# N.B.: HECUBA_ROOT is defined when you load a Hecuba HPC Module (or manually).
enqueue_compss -t -g -d --sc_cfg=mn.cfg \
  --qos="${QOS}" \
  --storage_props="${HECUBA_CONFIGURATION}" \
  --storage_home="${HECUBA_ROOT}/compss" \
  --job_name=esm_workflow \
  --exec_time=120 \
  --keep_workingdir \
  --worker_working_dir="${PWD}" \
  --worker_in_master_cpus="${CORES_PER_NODE}" \
  --num_nodes="${NODE_ALLOCATION}" \
  --pythonpath="${PWD}":"${HECUBA_ROOT}/compss" esm_simulation.py \
  --model "${MODEL}" \
  "${DEBUG}"