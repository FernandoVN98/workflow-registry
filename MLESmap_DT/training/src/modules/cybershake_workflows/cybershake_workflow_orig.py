from modules.cybershake_workflows.utils import SQLiteHandler, Ruptures
from modules.cybershake_workflows.steps import step_preSGT, step_preAWP, AWP_variable, checkAWP, step_postAWP, rupVar, step_runDS
import json
import os
import traceback
import sys
import glob


def simulation_site(config, runID, site):
    try:
        id = (runID % 2 + 1) * 2 - 1
        dal = SQLiteHandler(config["input"]["database"]["path"])
        # Create the directory for the current run
        cRunPath = config["output"]["path"] + "/" + site + "_" + str(runID)
        # Check if current site has started
        config["compute"]["restartFile"] = "stage.txt"
        os.makedirs(cRunPath, exist_ok=True)

        os.chdir(cRunPath)

        # Set the current path
        config["output"]["cRunpath"] = cRunPath

        # Obtain the Site ID for a given site name
        siteId = dal.getSiteID(site)

        # Check if the site do exist
        if not siteId:
            print("[ERROR] Skipping site '" + site + "': It was not found in the database", file=sys.stderr)
            return

        # Just a shortcut
        csetup = config["compute"]["setup"]

        # Insert the current run information onto the Database
        dal.addRunInfo(runID, siteId, config["input"]["ERF"]["id"],
                       config["input"]["model"]["id"],
                       csetup["sourceFrequency"], csetup["frequency"])

        # Store the ERF ID
        if not 'ruptures' in config["input"]["ERF"].keys():
            config["input"]["ERF"]["ruptures"] = cRunPath + "/ruptures/"
            if "focalMechanism" in config["compute"]["setup"].keys(): 
                # Generate rupture file
                rupt = Ruptures(config["input"]["ERF"]["ruptures"] + "Ruptures_erf"
                        + str(config["input"]["ERF"]["id"]) + "/",
                        config["input"]["database"]["path"])
                fm = config["compute"]["setup"]["focalMechanism"]
                rupt.generateRuptures(id, site, fm)
            else:
                # Generate rupture file
                rupt = Ruptures(config["input"]["ERF"]["ruptures"] + "Ruptures_erf"
                        + str(config["input"]["ERF"]["id"]) + "/",
                        config["input"]["database"]["path"])
                rupt.generateRuptures(id, site)

        # Write the CyberShake CFG file
        file_cfg = config["output"]["cRunpath"] + "/cybershake.cfg"
        with open(file_cfg, 'w') as f:
            f.write("CS_PATH = %s\n" % config["input"]["cyberShake"]["path"])
            f.write("SCRATCH_PATH = %s\n" % ("/scratch"))
            f.write("TMP_PATH = %s\n" % ("/tmp"))
            f.write("RUPTURE_ROOT = %s\n" % config["input"]["ERF"]["ruptures"])
            f.write("MPI_CMD = %s\n" % "srun")
            f.write("LOG_PATH = %s\n" % ("/logs"))
        generateCyberShakeCFG_compss(config, file_cfg)
        # Write the CyberShake CFG file
        siteSN = dal.getSiteShortName(site)
        sgtcordFileName, faultlistFileName, radiusFileName, fdlocFileName = step_preSGT(config, id, runID, siteSN, config["output"]["cRunpath"], 
                str(config["input"]["ERF"]["id"]), str(config["input"]["model"]["box"]),
               config["input"]["model"]["gridOut"], config["input"]["model"]["coords"], config["input"]["database"]["path"], 
               str(config["compute"]["setup"]["spacing"]), str(config["compute"]["setup"]["frequency"]), file_cfg,
               config["output"]["cRunpath"] + "/" + str(siteSN) + ".fdloc", config["output"]["cRunpath"] + "/" + str(siteSN) + ".faultlist", 
               config["output"]["cRunpath"] + "/" + str(siteSN) + ".radiusfile", config["output"]["cRunpath"] + "/" + str(siteSN) + ".coordfile")
        out_IN3D_x, out_IN3D_y, out_src_x, out_src_y, out_cordfile, out_media = step_preAWP(config["output"]["cRunpath"], 
                config["input"]["cyberShake"]["path"], config["input"]["model"]["path"], str(site) + "_" + str(runID),
               config["input"]["cyberShake"]["path"] + "/AWP-ODC-SGT/utils/build_awp_inputs.py", siteSN,
               config["input"]["model"]["gridOut"], config["compute"]["setup"]["frequency"],
               config["compute"]["decomposition"]["x"], config["compute"]["decomposition"]["y"],
               config["compute"]["decomposition"]["z"], config["compute"]["setup"]["sourceFrequency"],
               runID, config["input"]["model"]["path"], sgtcordFileName, fdlocFileName)

        dep1 = AWP_variable(out_IN3D_x, [out_src_x, out_cordfile, out_media])
        dep2 = AWP_variable(out_IN3D_y, [out_src_y, out_cordfile, out_media])

        dep1, awp_reformat_file, header_file = step_postAWP(config["input"]["cyberShake"]["path"], config["output"]["cRunpath"], "x", siteSN,
                config["input"]["model"]["box"], config["input"]["model"]["gridOut"], runID, config["compute"]["setup"]["frequency"],
                config["compute"]["setup"]["sourceFrequency"], config["input"]["cyberShake"]["gravesPitarka"], config["input"]["ERF"]["id"],
                config["input"]["database"]["path"], sgtcordFileName, fdlocFileName, out_media,  out_IN3D_x, dep1)
        dep2, awp_reformat_file, header_file = step_postAWP(config["input"]["cyberShake"]["path"], config["output"]["cRunpath"], "y", siteSN,
                config["input"]["model"]["box"], config["input"]["model"]["gridOut"], runID, config["compute"]["setup"]["frequency"],
                config["compute"]["setup"]["sourceFrequency"], config["input"]["cyberShake"]["gravesPitarka"], config["input"]["ERF"]["id"],
                config["input"]["database"]["path"], sgtcordFileName, fdlocFileName, out_media, out_IN3D_y, dep2)
        siteSN = dal.getSiteShortName(site)
        rupture_file_list = rupVar(config["input"]["cyberShake"]["path"], config["input"]["cyberShake"]["gravesPitarka"], 
                config["input"]["ERF"]["id"], config["input"]["database"]["path"],
                config["output"]["cRunpath"], config["input"]["region"], file_cfg, dep1, dep2)
        step_runDS(config["input"]["cyberShake"]["path"], config["input"]["database"]["path"], file_cfg, rupture_file_list, 
                config["output"]["cRunpath"], siteSN, runID, config["input"]["region"], config["input"]["ERF"]["ruptures"], 
                config["input"]["ERF"]["id"], config["input"]["cyberShake"]["gravesPitarka"], config["compute"]["setup"]["frequency"], 
                config)

    except Exception as error:
        print("Exception in code:")
        print('-'*80)
        traceback.print_exc(file=sys.stdout)
        print('-'*80)

def cybershake_workflow(configuration_file, **kwargs):
    try:
        # Read the configuration file
        with open(configuration_file, 'r') as f:
            config = json.load(f)

        # Create the base directory for the run and enter in it
        workSpace = os.path.abspath(config["output"]["path"])
        config["output"]["path"] = workSpace

        os.makedirs(workSpace, exist_ok=True)
        os.chdir(workSpace)

        # DataBase handling
        db = config["input"]["database"]

        # Create connection
        dal = SQLiteHandler(db["path"])

        # Check if the DB must be populated
        if db["populate"]:
            dal.importData(db["importFrom"])

        # Obtain a valid starting point RunID from DataBase
        runID = dal.getValidRunId()
        if config["compute"]["restart"]:
            tmpids = [i.split("_")[-1] for i in glob.glob(config["output"]["path"]+"/*[0-9]")]
            if tmpids:
              runID = int(min([i.split("_")[-1] for i in glob.glob(config["output"]["path"]+"/*[0-9]")]))

        # Obtain both model and ERF IDs
        modelID = dal.getModelID(config["input"]["model"]["name"])
        erfID = dal.getERFID(config["input"]["ERF"]["name"])

        # Check if both erfID and modelID were found
        foo =  "Please, check input JSON file '" + configuration_file + "' "
        if erfID == None:
            raise Exception("ERF not found." + foo)

        if modelID == None:
            raise Exception("Model not found." + foo)

        # Store the ERF ID and Model ID
        config["input"]["ERF"]["id"] = erfID
        config["input"]["model"]["id"] = modelID
        
        # For each Site defined as input
        runIDs = []
        i = 0
        for site in config["input"]["sites"]:
            runIDs.append(runID+i)
            i += 1
        for runID, site in zip(runIDs, config["input"]["sites"]):
            simulation_site(config, runID, site)
    except Exception as error:
        print("Exception in code:")
        print('-'*80)
        traceback.print_exc(file=sys.stdout)
        print('-'*80)
    
