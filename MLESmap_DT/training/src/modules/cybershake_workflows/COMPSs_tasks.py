from pycompss.api.task import task
from pycompss.api.mpi import mpi
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *
import os
from pycompss.api.multinode import multinode
import sys
import time
import sqlite3
import pymysql
import subprocess

from unifiedCSWFlow import DAL

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as ex:
        if ex.errno==errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def readCfg(file_input):
    varscfg={}
    try:
        filename =file_input#'%s/cybershake.cfg' % os.path.dirname(__file__)
        cfg = open(filename)
    except IOError:
        print("%s not found.\n" % filename)
        sys.exit(-2)
    cfgContents = cfg.readlines()
    for line in cfgContents:
        if line[0]!='#':
            pieces = line.split('=')
            varscfg[pieces[0].strip()] = pieces[1].strip()
    return varscfg

def getProperty(property, file_input):
    #if len(vars)==0:
    varscfg = readCfg(file_input)
    try:
        propertyVal = varscfg[property]
    except KeyError:
        print("No %s found in cybershake.cfg.\n" % property)
        sys.exit(-1)
    return propertyVal


@task(returns=2)
def getSiteCoords(site, db_file):
    connection = sqlite3.connect(db_file)
    cursor = connection.cursor()
    sql_string = 'select CS_Site_Lat, CS_Site_Lon from CyberShake_Sites where CS_Short_Name="%s"' % site
    cursor.execute(sql_string)
    siteCoords = cursor.fetchone()
    return float(siteCoords[0]), float(siteCoords[1])

@task(outputName=FILE_OUT, cordfileName=FILE_IN, returns=1)
def genFdloc(outputName, site, mlon, mlat, cordfileName):
    cordfile = open(cordfileName)
    minDistance = 1000.0
    minLoc = [-1, -1]
    i = 0
    for line in cordfile:
        i = i+1
        if i % 100000==0:
            print('%d' % i)
        pieces = line.split()
        distance = (float(pieces[0])-mlon)*(float(pieces[0])-mlon)+(float(pieces[1])-mlat)*(float(pieces[1])-mlat)
        if distance < minDistance:
            minDistance = distance
            minLoc[0] = pieces[2]
            minLoc[1] = pieces[3]
    cordfile.close()
    fdlocFile = open(outputName, 'w')
    fdlocFile.write('%s %s\n' % (minLoc[0], minLoc[1]))
    fdlocFile.flush()
    fdlocFile.close()
    return [int(minLoc[0]), int(minLoc[1])]

def fault_list_script(site, erf_id, rupture_path, output_file, db_file, rsqsim, *args):
    if (site==None or erf_id==None or rupture_path==None or output_file==None):
        print("site, erf_id, rup_path, and output are all required.")
        parser.print_help()
        sys.exit(1)
    server = "moment.usc.edu"
    if db_file is not None:
        server = db_file
    db = "CyberShake"
    user = "cybershk_ro"
    passwd = "CyberShake2007"
    header_rows = 6
    if rsqsim == True:
        header_rows = 5
    #if server[0:9]=="sqlite://":
    #    conn = sqlite3.connect(server[9:])
    #else:
    #    conn = pymysql.connect(host=server, user=user, passwd=passwd, db=db)
    conn = sqlite3.connect(server)
    cur = conn.cursor()
    query = 'select R.Source_ID, R.Rupture_ID from CyberShake_Site_Ruptures R, CyberShake_Sites S where S.CS_Short_Name="%s" and R.CS_Site_ID=S.CS_Site_ID and R.ERF_ID=%d order by R.Source_ID, R.Rupture_ID' % (site, int(erf_id))
    cur.execute(query)
    res = cur.fetchall()
    if len(res)==0:
        print("No ruptures found for site %s." % site)
        sys.exit(2)

    fp_out = open(output_file, "w")
    count = 0
    for r in res:
        sourceID = int(r[0])
        ruptureID = int(r[1])
        fp_out.write("%s/%d/%d/%d_%d.txt nheader=%d latfirst=1\n" % (rupture_path, int(sourceID), int(ruptureID), int(sourceID), int(ruptureID), int(header_rows)))
        count += 1
        if (count%100==0):
            print("Processed %d ruptures." % count)
    fp_out.flush()
    fp_out.close()

    cur.close()
    conn.close()

@task(outputName=FILE_OUT)
def genFaultList(outputName, site, erf_id, rupt_geometries, rsqsim, db_file):
    '''Copies the functionality of gen_faultlist.csh:  it serves as a wrapper to CreateFaultList.java, which queries the database to generate a list of ruptures which are applicable for the given site.'''
    rsqsim_str = ""
    if rsqsim:
        rsqsim_str = "-rsqsim"
    #command = '%s/faultlist_py/CreateFaultList.py %s %s %s %s %s %s' % (sys.path[0], site, erf_id, PATH_TO_RUPTURE_GEOMETRIES, outputName, db_file, rsqsim_str)
    #returnCode = os.system(command)
    returnCode = fault_list_script(site, erf_id, rupt_geometries, outputName, db_file, rsqsim_str=="-rsqsim")
    #if returnCode!=0:
    #    sys.exit((returnCode >> 8) & 0xFF)

@task(radiusfile=FILE_OUT)
def genRadiusFile(radiusfile):
    '''Duplicates the functionality of part of gen_sgtgrid.csh:  it writes the adaptive mesh info to a radius file as part of generating the cordfile.'''
    #magic constants for the adaptive mesh
    RLEV = [10.0, 50.0, 100.0, 1000.0]
    RINC = [10, 15, 25, 50]
    ZLEV = [0.2, 5.0, 24.0, 60.0]
    ZINC = [1, 5, 10, 25]

    output = open(radiusfile, 'w')
    output.write('%d\n' % len(RLEV))
    for r in RLEV:
        output.write('%f ' % r)
    output.write('\n')
    for r in RINC:
        output.write('%d ' % r)
    output.write('\n')
    output.write('%d\n'% len(ZLEV))
    for z in ZLEV:
        output.write('%f ' % z)
    output.write('\n')
    for z in ZINC:
        output.write('%d ' % z)
    output.write('\n')
    output.flush()
    output.close()

@task(sgtcordFileName=FILE_OUT, faultlistFileName=FILE_IN, radiusFileName=FILE_IN, cfg_file=FILE_IN)
def genSgtGrid(sgtcordFileName, cybershake_path, site, ns, src, mlon, mlat, mrot, faultlistFileName, radiusFileName, frequency, spacing, cfg_file):
    '''Copies the functionality of gen_sgtgrid.csh:  it generates a list of grid points that the SGTs should be saved for.  This includes an adaptive mesh as well as the locations of the ruptures.  Essentially this serves as a wrapper for gen_sgtgrid.c.'''
    #magic constants
    HH = 0.1/float(frequency)
    if spacing>0:
        HH = spacing
    IX_MIN = 20
    IX_MAX = ns[0]-20
    IY_MIN = 20
    IY_MAX = ns[1]-20
    IZ_START = 1
    IZ_MAX = ns[2]-20

    #outfile = '../../data/SgtInfo/SgtCords/%s.cordfile' % site
    #faultlist = '../../data/SgtInfo/FaultList/%s.faultlist' % site
    #split into subfiles
    MPI_CMD = getProperty('MPI_CMD', cfg_file)

    if (MPI_CMD == "mpirun"):
        try:
            node_file = os.environ["PBS_NODEFILE"]
            num_nodes = 0
            f = open(node_file)
            lines = f.readlines()
            num_nodes = len(lines)
            MPI_CMD = "%s -np %d -machinefile %s" % \
                    (MPI_CMD, num_nodes, node_file)
        except:
                print("Unable to read nodefile %s" % (node_filw))
                sys.exit(1)
    elif (MPI_CMD == "aprun"):
        np = int(os.environ["PBS_NUM_NODES"])*4
        #No more than 32 cores)
        np = min(np, 32)
        MPI_CMD = "%s -n %d -N 4" % (MPI_CMD, np)

    command = '%s %s/PreSgt/bin/gen_sgtgrid nx=%d ny=%d nz=%d h=%f xsrc=%d ysrc=%d ixmin=%d ixmax=%d iymin=%d iymax=%d izstart=%d izmax=%d radiusfile=%s outfile=%s modellon=%f modellat=%f modelrot=%f faultlist=%s' % (MPI_CMD, cybershake_path, int(ns[0]), int(ns[1]), int(ns[2]), HH, int(src[0]), int(src[1]), IX_MIN, int(IX_MAX), IY_MIN, int(IY_MAX), IZ_START, int(IZ_MAX), radiusFileName, sgtcordFileName, float(mlon), float(mlat), float(mrot), faultlistFileName)
    startTime = time.time()
    returnCode = os.system(command)
    print("Elapsed time: %f\n" % (time.time()-startTime))
    if returnCode!=0:
        sys.exit((returnCode >> 8) & 0xFF)

def steppreSGT(config, id, runID, site, run_path, erf_id, modelbox, gridout, cordfileName, db_file, setup_spacing, frequency, cfg_file_input,
        fdlocFileName, faultlistFileName, radiusFileName, sgtcordFileName):

    spacing = float(setup_spacing)
    rsqsim = True #CHEQUEAR ESTO
    if frequency is None:
        frequency = 0.5
    PATH_TO_RUPTURE_GEOMETRIES = "%s/Ruptures_erf%s" % (getProperty('RUPTURE_ROOT', cfg_file_input), erf_id)
    input = open(modelbox)
    modelboxContents = input.readlines()
    input.close()
    #site = modelboxContents[0].strip()
    
    siteLat, siteLon = getSiteCoords(site, db_file)
    
    modelTokens = modelboxContents[4].split()
    mlon = (float)(modelTokens[1])
    mlat = (float)(modelTokens[3])
    mrot = (float)(modelTokens[5])
    
    src = genFdloc(fdlocFileName, site, siteLon, siteLat, cordfileName)
    genFaultList(faultlistFileName, site, erf_id, PATH_TO_RUPTURE_GEOMETRIES, rsqsim, db_file)
    
    genRadiusFile(radiusFileName)
    
    input = open(gridout)
    gridoutContents = input.readlines();
    input.close()
    ns = []
    ns.append(int((gridoutContents[1].split("="))[1]))
    ns.append(int((gridoutContents[1+ns[0]+2].split("="))[1]))
    ns.append(int((gridoutContents[1+ns[0]+2+ns[1]+2].split("="))[1]))
    #CAFE 12062023
    genSgtGrid(sgtcordFileName, config["input"]["cyberShake"]["path"], site, ns, src, mlon, mlat, mrot, faultlistFileName, radiusFileName, frequency, spacing,
            cfg_file_input)
    return sgtcordFileName, faultlistFileName, radiusFileName, fdlocFileName

@constraint(computing_units="112")#${ComputingUnits}")
@multinode(computing_nodes="1")
@task(returns=1)
def preSGT(config, id, runID, site, run_path, launching_script, ERF_id, model_box, model_Gridout, model_coords, database_path,
           setup_spacing, setup_frequency):
    cmd = (str(launching_script) + " " + str(site) + " " + str(ERF_id) + " " + str(model_box) + " " +
           str(model_Gridout) + " " + str(model_coords) + " " + str(run_path+"/"+site)+".fdloc" + " " + str(run_path+"/"+site)+".faultlist" +
           " " + str(run_path+"/"+site)+".radiusfile" + " " + str(run_path+"/"+site)+".coordfile" + " " + str(database_path) +
           " " + str(setup_spacing) + " " + str(setup_frequency) + " " + str(run_path+"/ruptures"))
    os.system(cmd)
    return 1

@task(gridout=FILE_IN, awp_media=FILE_INOUT)#media_out=FILE_OUT)
def build_media(cybershake_path, site, gridout, rwg_vel_prefix, awp_media):#, media_out):
    fp_in = open(gridout, "r")
    data = fp_in.readlines()
    fp_in.close()
    nx = int((data[1].split("="))[1])
    ny = int((data[1+nx+2].split("="))[1])
    nz = int((data[1+nx+2+ny+2].split("="))[1])
    cmd = "%s/SgtHead/bin/reformat_velocity %d %d %d %s %s" % (cybershake_path, nx, ny, nz, rwg_vel_prefix, awp_media)
    print(cmd)
    exitcode = os.system(cmd)
    fp_in = open(awp_media, "rb")
    data = fp_in.read()
    fp_in.close()
    fp_out = open(awp_media, "wb")
    fp_out.write(data)
    fp_out.flush()
    fp_out.close()

@task(output_file=FILE_OUT, fdloc=FILE_IN, input_file=FILE_IN)
def build_src(output_file, site, fdloc, input_file, awp_comp, frequency, path, filter=None, spacing=None):
    if awp_comp=='x':
        comp = "y"
    elif awp_comp=='y':
        comp = "x"
    elif awp_comp=='z':
        comp = "z"
    else:
        print("Error:  component %s not recognized, aborting." % comp)
        sys.exit(1)
    if spacing is not None:
        #Always simulate 200.0 seconds of time
        #dt = spacing in km/20.0
        dt = spacing/20.0
        nt = int(200.0/dt)
        #Round to nearest 1000
        if (nt % 1000)!=0:
            nt = 1000*((nt/1000)+1)
    else:
        nt = int(frequency*40000.0)
    fp_in = open(fdloc, "r")
    data = fp_in.readline()
    fp_in.close()
    [src_x, src_y] = data.split()
    if filter == None:
        filter = frequency
    fp_in = open(input_file, "r")
    src_data = fp_in.readlines()
    fp_in.close()
    fp_out = open(output_file, "w")
    fp_out.write("%d %d 1\n" % (int(src_y)+1, int(src_x)+1))
    for line in src_data:
        fp_out.write(line)
    fp_out.flush()
    fp_out.close()


@task(output_file=FILE_OUT, input_file=FILE_IN)
def build_IN3D(output_file, input_file, site, gridout, awp_comp, frequency, procs, path_station, velocity_mesh, vel_pref, spacing):
    fp_in = open(input_file)
    data = fp_in.readlines()
    fp_in.close()
    param = dict()
    for line in data:
        pieces = line.split("=")
        param[pieces[0]] = pieces[1].strip()
    #determine igreen
    igreen = 0
    if awp_comp=='x':
        #swap x and y
        comp = 'y'
        igreen = 4
    elif awp_comp=='y':
        #swap x and y
        comp = 'x'
        igreen = 5
    elif awp_comp=='z':
        comp = 'z'
        igreen = 6
    else:
        print("Component %s not recognized, aborting." % comp)
        sys.exit(2)
    param["igreen"] = igreen
    #determine DH, DT, NST, READ_STEP, WRITE_STEP, FP
    if spacing is not None:
        param["DH"] = int(1000.0*spacing)
    else:
        param["DH"] = round(100.0/frequency, 1)
    if spacing is not None:
        param["DT"] = spacing/20.0
    else:
        param["DT"] = 0.005/frequency

    SIMULATED_TIME = float(param["TMAX"])
    #Round up to nearest 1000
    nst = int(SIMULATED_TIME/param["DT"])
    if (nst % 1000)!=0:
        nst = 1000*(nst/1000 + 1)
    param["NST"] = nst
    #Change TMAX based on DT and NST
    param["TMAX"] = param["NST"]*param["DT"]
    #Talked to Kim and Rob, FP should remain 0.5
    param["FP"] = 0.5
    param["READ_STEP"] = param["NST"]
    #Divide by 10 because WRITE_STEP is in units of # of steps being written, not total # of timesteps
    #So if 20000 simulated timesteps and decimation of 10, WRITE_STEP = 2000
    param["WRITE_STEP"] = int(param["READ_STEP"]/int(param["NTISKP_SGT"]))#int(param["READ_STEP"]/int(param["NTISKP_SGT"]))
    param["WRITE_STEP2"] = int(param["WRITE_STEP"])#int(param["WRITE_STEP"])
    #determine NX, NY, NZ
    fp_in = open(gridout, "r")
    data = fp_in.readlines()
    fp_in.close()
    #remember, X and Y are flipped
    ny = int((data[1].split("="))[1])
    nx = int((data[1+ny+2].split("="))[1])
    nz = int((data[1+ny+2+nx+2].split("="))[1])

    param["NX"] = nx
    param["NY"] = ny
    param["NZ"] = nz
    #Check proc values to make sure they are evenly divisible
    if nx % procs[0] != 0:
        print("PX %d must be a factor of NX %d, aborting." % (proc[0], nx))
        sys.exit(2)
        if ny % procs[1] != 0:
            print("PY %d must be a factor of NY %d, aborting." % (proc[1], ny))
            sys.exit(2)
        if nz % procs[2] != 0:
            print("PZ %d must be a factor of NZ %d, aborting." % (proc[2], nz))
            sys.exit(2)
    param["NPX"] = procs[0]
    param["NPY"] = procs[1]
    param["NPZ"] = procs[2]
    #doesn't really matter, but set N<BG|ED><comp> values
    param["NBGX"] = 1
    param["NEDX"] = nx
    param["NBGY"] = 1
    param["NEDY"] = ny
    param["NBGZ"] = 1
    param["NEDZ"] = nz
    param["NBGX2"] = 1
    param["NEDX2"] = nx
    param["NBGY2"] = 1
    param["NEDY2"] = ny
    param["NBGZ2"] = 1
    param["NEDZ2"] = nz
    #paths to INSRC, INVEL, INSGT, SGTGR0
    param["INSRC"] = "%s/comp_%s/input/%s_f%s_src" % (path_station, awp_comp, site, awp_comp)
    param["INVEL"] = "%s/awp.%s.media" % (path_station, site)
    param["INSGT"] = "%s/comp_%s/input/awp.%s.cordfile" % (path_station, awp_comp, site)
    param["SGTGRO"] = "%s/comp_%s/output_sgt/awp-strain-%s-f%s" % (path_station, awp_comp, site, awp_comp)
    fp_out = open(output_file, "w")
    fp_out.write("%9s igreen\n" % (param["igreen"]))
    fp_out.write("%9s TMAX\n\n" % (param["TMAX"]))
    fp_out.write("%9s DH\n" % (param["DH"]))
    fp_out.write("%9s DT\n\n" % (param["DT"]))
    fp_out.write("%9s NPC\n\n" % (param["NPC"]))
    fp_out.write("%9s ND\n" % (param["ND"]))
    fp_out.write("%9s ARBC\n" % (param["ARBC"]))
    fp_out.write("%9s PHT\n\n" % (param["PHT"]))
    fp_out.write("%9s NSRC\n" % (param["NSRC"]))
    fp_out.write("%9s NST\n\n" % (param["NST"]))
    fp_out.write("%9d NX\n" % (param["NX"]))
    fp_out.write("%9d NY\n" % (param["NY"]))
    fp_out.write("%9d NZ\n\n" % (param["NZ"]))
    fp_out.write("%9d NPX\n" % (param["NPX"]))
    fp_out.write("%9d NPY\n" % (param["NPY"]))
    fp_out.write("%9d NPZ\n\n" % (param["NPZ"]))
    fp_out.write("%9s IFAULT\n" % (param["IFAULT"]))
    fp_out.write("%9s CHECKPOINT\n" % (param["CHECKPOINT"]))
    fp_out.write("%9s ISFCVLM\n" % (param["ISFCVLM"]))
    fp_out.write("%9s IMD5\n" % (param["IMD5"]))
    fp_out.write("%9s IVELOCITY\n" % (param["IVELOCITY"]))
    fp_out.write("%9s MEDIARESTART\n" % (param["MEDIARESTART"]))
    fp_out.write("%9s NVAR\n" % (param["NVAR"]))
    fp_out.write("%9s IOST\n" % (param["IOST"]))
    fp_out.write("%9s PARTDEG\n" % (param["PARTDEG"]))
    fp_out.write("%9s IO_OPT\n" % (param["IO_OPT"]))
    fp_out.write("%9s PERF_MEAS\n" % (param["PERF_MEAS"]))
    fp_out.write("%9s IDYNA\n" % (param["IDYNA"]))
    fp_out.write("%9s SOCALQ\n\n" % (param["SOCALQ"]))
    fp_out.write("%9s NVE\n\n" % (param["NVE"]))
    fp_out.write("%9s MU_S\n" % (param["MU_S"]))
    fp_out.write("%9s MU_D\n\n" % (param["MU_D"]))
    fp_out.write("%9s FL\n" % (param["FL"]))
    fp_out.write("%9s FH\n" % (param["FH"]))
    fp_out.write("%9s FP\n\n" % (param["FP"]))
    fp_out.write("%9s READ_STEP\n" % (param["READ_STEP"]))
    fp_out.write("%9s WRITE_STEP\n" % (param["WRITE_STEP"]))
    fp_out.write("%9s WRITE_STEP2\n\n" % (param["WRITE_STEP2"]))
    fp_out.write("%9s NBGX\n" % (param["NBGX"]))
    fp_out.write("%9s NEDX\n" % (param["NEDX"]))
    fp_out.write("%9s NSKPX\n" % (param["NSKPX"]))
    fp_out.write("%9s NBGY\n" % (param["NBGY"]))
    fp_out.write("%9s NEDY\n" % (param["NEDY"]))
    fp_out.write("%9s NSKPY\n" % (param["NSKPY"]))
    fp_out.write("%9s NBGZ\n" % (param["NBGZ"]))
    fp_out.write("%9s NEDZ\n" % (param["NEDZ"]))
    fp_out.write("%9s NSKPZ\n\n" % (param["NSKPZ"]))
    #NTISKP needs to be larger than # of timesteps so that no velocity data is written
    fp_out.write("%9s NTISKP\n" % (2*int(param["NST"])))
    fp_out.write("%9s NBGX2\n" % (param["NBGX2"]))
    fp_out.write("%9s NEDX2\n" % (param["NEDX2"]))
    fp_out.write("%9s NSKPX2\n" % (param["NSKPX2"]))
    fp_out.write("%9s NBGY2\n" % (param["NBGY2"]))
    fp_out.write("%9s NEDY2\n" % (param["NEDY2"]))
    fp_out.write("%9s NSKPY2\n" % (param["NSKPY2"]))
    fp_out.write("%9s NBGZ2\n" % (param["NBGZ2"]))
    fp_out.write("%9s NEDZ2\n" % (param["NEDZ2"]))
    fp_out.write("%9s NSKPZ2\n\n" % (param["NSKPZ2"]))
    fp_out.write("%9s NTISKP2\n\n" % (param["NTISKP2"]))
    fp_out.write("%9s NTISKP_SGT\n\n" % (param["NTISKP_SGT"]))
    fp_out.write("'%s/comp_%s/%s' CHKP\n" % (path_station, awp_comp, param["CHKP"]))
    fp_out.write("'%s/comp_%s/%s' CHKJ\n\n" % (path_station, awp_comp, param["CHKJ"]))
    fp_out.write("'%s' INSRC\n" % (param["INSRC"].split("/")[-1]))
    if vel_pref is not None:
        fp_out.write("'%s' INVEL\n\n" % (param["INVEL"].split("/")[-1]))
    else:
        fp_out.write("'%s' INVEL\n\n" % (velocity_mesh))#(param["INVEL"].split("/")[-1]))
    fp_out.write("'%s' INSGT\n\n" % (param["INSGT"].split("/")[-1]))
    fp_out.write("'%s/comp_%s/%s' SXRGO\n" % (path_station, awp_comp, param["SXRGO"]))
    fp_out.write("'%s/comp_%s/%s' SYRGO\n" % (path_station, awp_comp, param["SYRGO"]))
    fp_out.write("'%s/comp_%s/%s' SZRGO\n\n" % (path_station, awp_comp, param["SZRGO"]))
    fp_out.write("'%s/comp_%s/%s' SXRGO2\n" % (path_station, awp_comp, param["SXRGO2"]))
    fp_out.write("'%s/comp_%s/%s' SYRGO2\n" % (path_station, awp_comp, param["SYRGO2"]))
    fp_out.write("'%s/comp_%s/%s' SZRGO2\n\n" % (path_station, awp_comp, param["SZRGO2"]))
    fp_out.write("'%s' SGTGRO\n\n" % (param["SGTGRO"]))
    fp_out.write("'%s/comp_%s/%s' SGSN\n" % (path_station, awp_comp, param["SGSN"]))

    fp_out.flush()
    fp_out.close()

@task(cordfile=FILE_IN, awp_cordfile=FILE_OUT)
def build_cordfile(cordfile, awp_cordfile, max_depth_index):
    fp_in = open(cordfile, "r")
    fp_out = open(awp_cordfile, "w")
    max_depth_index = int(max_depth_index)

    start_time = time.time()
    data = fp_in.readlines()
    fp_in.close()
    end_time = time.time()
    print("%f sec to read input file." % (end_time-start_time))

    points = dict()

    num_pts_in = int(data[4])

    start_time = time.time()
    for i in range(5, len(data), 1):
        line = data[i]
        pieces = line.split()
        #Encode key like 4th column (x*100000*10000 + y*10000 + z)
        #Add 1 to x, y, and z
        x = int(pieces[0])+1
        y = int(pieces[1])+1
        z = min([int(pieces[2])+1, max_depth_index])
        pt = x*100000*10000 + y*10000 + z
        if pt in points:
            print("Duplicate point entry %s" % point_str)
        else:
            points[pt] = 1
        if (i-4)%500000==0:
            end_time = time.time()
            print("%f sec to read lines %d to %d" % ((end_time-start_time), (i-500004), (i-4)))
            start_time = time.time()

    start_time = time.time()
    p_list = sorted(points)
    end_time = time.time()
    print("%f sec to sort list." % (end_time-start_time))
    fp_out.write("%d\n" % len(p_list))
    start_time = time.time()
    for entry in p_list:
        #Print as (y, x, z)
        x = entry/(100000*10000)
        y = (entry/10000) % 100000
        z = entry % 10000
        #1 already added, don't need it here
        fp_out.write("%d %d %d\n" % (y, x, z))
    fp_out.flush()
    fp_out.close()
    end_time = time.time()
    print("%f sec to write output." % (end_time-start_time))

def build_awp_inputs(run_path, cybershake_path, site, gridout, fdloc, cordfile, px, py, pz, velocity_mesh, frequency, source_frequency, spacing, rwg_vel_prefix=None):
    if site==None or gridout==None or fdloc==None or cordfile==None:
        print("site, gridout, fdloc, and cordfile must be specified.")
        parser.print_help()
        sys.exit(1)
        return
    procs = [px, py, pz]
    if (procs[0]==None or procs[1]==None or procs[2]==None):
        print("px, py, pz must be specified.")
        parser.print_help()
        sys.exit(1)
    source_freq = frequency
    if source_frequency is not None:
        source_freq = source_frequency
    awp_comps = ['x', 'y']
    out_IN3D_x = run_path+"/IN3D."+str(site)+".x"
    out_IN3D_y = run_path+"/IN3D."+str(site)+".y"
    out_src_x = run_path+"/comp_x/input/"+site+"_fx_src"
    out_src_y = run_path+"/comp_y/input/"+site+"_fy_src"
    for c in awp_comps:
        mkdir_p(run_path+"/comp_%s/input" % c)
        mkdir_p(run_path+"/comp_%s/output_ckp" % c)
        mkdir_p(run_path+"/comp_%s/output_sfc" % c)
        mkdir_p(run_path+"/comp_%s/output_vlm" % c)
        mkdir_p(run_path+"/comp_%s/output_sgt" % c)
        #exitcode = os.system("%s setstripe -c 160 -s 5m comp_%s/output_sgt" % (LFS_PATH, c))
        #if exitcode!=0:
        #    print("Error striping with command %s setstripe -c 160 -s 5m comp_%s/output_sgt, exiting." % (LFS_PATH, c))
        print("Building IN3D file for comp %s." % c)
        if c == "x":
            build_IN3D(out_IN3D_x, cybershake_path+"/AWP-ODC-SGT/utils/data/IN3D.ref", site, gridout, c, frequency, procs, run_path, velocity_mesh, rwg_vel_prefix, spacing=spacing)
        elif c == "y":
            build_IN3D(out_IN3D_y, cybershake_path+"/AWP-ODC-SGT/utils/data/IN3D.ref", site, gridout, c, frequency, procs, run_path, velocity_mesh, rwg_vel_prefix, spacing=spacing)
        #out_src = run_path+"/comp_"+c+"/input/"+site+"_f"+c+"_src"
        if spacing is not None:
            #Always simulate 200.0 seconds of time
            #dt = spacing in km/20.0
            dt = spacing/20.0
            nt = int(200.0/dt)
            #Round to nearest 1000
            if (nt % 1000)!=0:
                nt = 1000*((nt/1000)+1)
        else:
            nt = int(frequency*40000.0)
        if c == "x":
            build_src(out_src_x, site, fdloc, cybershake_path+"/AWP-ODC-SGT/utils/data/f"+c+"_src_"+str(nt)+"_"+str(source_frequency)+"hzFilter", c, frequency, run_path, filter=source_frequency, spacing=spacing)
        if c == "y":
            build_src(out_src_y, site, fdloc, cybershake_path+"/AWP-ODC-SGT/utils/data/f"+c+"_src_"+str(nt)+"_"+str(source_frequency)+"hzFilter", c, frequency, run_path, filter=source_frequency, spacing=spacing)
    awp_cordfile_with_path = run_path + "/awp.%s.cordfile" % site
    awp_cordfile = "awp.%s.cordfile" % site
    print("Building cordfile.")
    sys.stdout.flush()
    #Determine max depth index from gridout file
    with open(gridout, "r") as fp_in:
        lines = fp_in.readlines()
        max_depth_index = int(lines[-1].split()[0].strip()) + 1
        fp_in.close()

    rc = build_cordfile(cordfile, awp_cordfile_with_path, max_depth_index)
    for c in awp_comps:
        if os.path.lexists("%s/comp_%s/input/%s" % (run_path, c, awp_cordfile)):
            os.remove("%s/comp_%s/input/%s" % (run_path, c, awp_cordfile))
        os.symlink("%s" % awp_cordfile_with_path, "%s/comp_%s/input/%s" % (run_path, c, awp_cordfile))
    awp_media = run_path+"/awp.%s.media" % (site)
    awp_media_out = "awp.%s.media" % (site)
    if rwg_vel_prefix is not None:
        build_media(cybershake_path, site, gridout, rwg_vel_prefix, awp_media)#, awp_media_out)
        #for c in awp_comps:
        #    if os.path.lexists("%s/comp_%s/input/%s" % (run_path, c, awp_media_out)):
        #        os.remove("%s/comp_%s/input/%s" % (run_path, c, awp_media_out))
        #    #Extra level back because of AWP_SGT_<site> subdirector
        #    os.symlink("%s" % awp_media, "%s/comp_%s/input/%s" % (run_path, c, awp_media_out))
    else:
        if not os.path.exists(awp_media):
            print("Error, since expected velocity file %s does not exist.  Aborting." % awp_media)
            sys.exit(3)
        print("No velocity prefix specified, skipping velocity file reformat.")
        awp_media = velocity_mesh
    for c in awp_comps:
        if os.path.lexists("%s/comp_%s/input/%s" % (run_path, c, awp_media_out)):
            os.remove("%s/comp_%s/input/%s" % (run_path, c, awp_media_out))
        #Extra level back because of AWP_SGT_<site> subdirector
        os.symlink("%s" % awp_media, "%s/comp_%s/input/%s" % (run_path, c, awp_media_out))
    return out_IN3D_x, out_IN3D_y, out_src_x, out_src_y, awp_cordfile_with_path, awp_media

def steppreAWP(run_path, cybershake_path, model_path, site_compound, launching_script, site, grid_out, frequency, px, py, pz, source_frequency, runID, velocity_mesh, cordfile, fdloc):
    os.symlink(model_path, str(run_path) + "/awp." + site + ".media")
    os.chdir(run_path)
    spacing=None
    out_IN3D_x, out_IN3D_y, out_src_x, out_src_y, out_cordfile, out_media = build_awp_inputs(run_path, cybershake_path, site, grid_out, fdloc, cordfile, px, py, pz, velocity_mesh, frequency, source_frequency, spacing)
    return out_IN3D_x, out_IN3D_y, out_src_x, out_src_y, out_cordfile, out_media

@constraint(computing_units="112")#${ComputingUnits}")
@multinode(computing_nodes="1")
@task(returns=1)
def preAWP(run_path, model_path, site_compound, launching_script, site, grid_out, frequency, x, y, z, source_frequency, runID, velocity_mesh, dependency):
    bash_command = "ln -s " + model_path + " " + str(run_path) + "/awp." + site + ".media"
    os.system(bash_command)
    bash_command = "cd " + run_path
    os.system(bash_command)
    cmd = (str(launching_script) + " --site " + str(site) + " --gridout " + str(grid_out) + " --fdloc " + str(run_path+"/"+site) +".fdloc" +
           " --cordfile " + str(run_path+"/"+site)+ ".coordfile" + " --frequency " + str(frequency) + " --px " + str(x) + " --py "
           + str(y) + " --pz " + str(z) + " --source-frequency " + str(source_frequency) + " --run_id " + str(runID) + " --velocity-mesh " + str(velocity_mesh)
           + " --act-dir " + str(run_path))
    os.system(cmd)
    return 1


@mpi(binary='/gpfs/scratch/bsc19/bsc019756/CyberShake/software/AWP-ODC-SGT/bin/pmcl3d', args="{{input_file}}", processes_per_node=112, processes=336, runner='srun')#, computing_nodes="${ComputingNodesAWP}")
@task(input_file=FILE_IN, col_file=COLLECTION_FILE_IN, media_file=FILE_IN, returns=1)
def AWP_task_non_dep(input_file, col_file, media_file):
    pass


@mpi(binary='/gpfs/scratch/bsc19/bsc019756/CyberShake/software/AWP-ODC-SGT/bin/pmcl3d', working_dir="{{work_dir}}", args="{{input_file}}", processes_per_node=112, processes=336, runner='srun')#, computing_nodes="${ComputingNodesAWP}")
@task(returns=1)
def AWP_task(work_dir, input_file, dep_pre):
    pass

@task(returns=1)
def check_AWP(cybershake_path, run_path, variable, site, ined_file):
    if variable == "x":
        op_variable = "y"
    else:
        op_variable = "x"
    cmd = "cd " + run_path
    os.system(cmd)
    cmd = cybershake_path + "/SgtTest/perform_checks.py " + run_path + "/comp_" + variable + "/output_sgt/awp-strain-" + site + "-f" + variable + " " + run_path + "/" + site + ".coordfile " + run_path + "/" + ined_file
    value = os.system(cmd)
    return value

@constraint(computing_units="112")
@multinode(computing_nodes="1")
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

@constraint(computing_units="112")
@multinode(computing_nodes="1")
@task(cordfile=FILE_IN, fdloc=FILE_IN, gridout=FILE_IN, modelbox=FILE_IN, returns=1)# header_file_out=FILE_OUT, returns=1)
def write_head(cybershake_path, modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, source_freq, media, header_name, header_file_out, dep):
    #On Titan, need aprun to execute this on a compute node, not the aprun node
    command = "srun %s/SgtHead/bin/write_head %s %s %s %s %f %d %f %d %s %s %f %s %s" % (cybershake_path, modelbox, cordfile, fdloc, gridout, spacing, nt, dt, decimation, comp, moment, source_freq, media, header_name)
    if not comp=="z":
        #Add c flag to use corrected mu
        command = "%s -c " % command
    print("COMMAND", flush=True)
    print(command, flush=True)
    exitcode = os.system(command)
    return exitcode
    '''fp_in = open(header_name, "rb")
    data = fp_in.read()
    fp_out = open(header_file_out, "wb")
    fp_out.write(data)
    fp_out.close()
    fp_in.close()'''


@task(cordfile=FILE_IN, in3d_file=FILE_IN, returns=6)
def read_previous_files(cordfile, in3d_file):
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

def post_AWP_no_task(cybershake_path, run_path, variable, site, model_box, model_gridout, runID, frequency, source_frequency, gravesPitarka, erf_id, database_path, cordfile, fdloc, media, in3d_file, dependency, skip_md5=False):
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
    total_ts, decimation, spacing, total_ts_decimation, num_sgt_pts, params_dt = read_previous_files(cordfile, in3d_file)
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

@constraint(computing_units="112")
@multinode(computing_nodes="1")
@task(returns=1)
def postAWP(cybershake_path, run_path, variable, site, ined_file, model_box, model_gridout, runID, frequency, source_frequency, gravesPitarka, erf_id, database_path, dependency):
    if variable == "x":
        op_variable = "y"
    else:
        op_variable = "x"
    cmd = "cd " + run_path
    os.system(cmd)
    cmd = cybershake_path + "/SgtTest/perform_checks.py " + run_path + "/comp_" + variable + "/output_sgt/awp-strain-" + site + "-f" + variable + " " + run_path + "/" + site + ".coordfile " + run_path + "/" + ined_file
    os.system(cmd)
    cmd = (cybershake_path + "/AWP-GPU-SGT/utils/prepare_for_pp.py " + site + " " + run_path + "/comp_" + variable + "/output_sgt/awp-strain-" + site + "-f" + variable +
           " " + run_path + "/" + site + "_f" + op_variable + "_1.sgt " + model_box + " "
           + run_path + "/"+site+".coordfile " + run_path + "/"+site+".fdloc " + model_gridout + " " + run_path + "/" + ined_file + " " + run_path + "/" + "awp." + site + ".media " + variable + " " +
           str(runID) + " " + run_path + "/" + site + "_f" + op_variable + "_" + str(runID) + ".sgthead " + str(frequency) + " -s " + str(source_frequency) + " " + run_path)
    os.system(cmd)
    print("COMMAND", flush=True)
    print(cmd, flush=True)
    return 1


@constraint(computing_units="112")
@multinode(computing_nodes="1")
@task(rupture_file=FILE_IN, rupture_out=FILE_OUT, returns=1)
def populate_rvs(cybershake_path, gravesPitarka, erf_id, database_path, run_path, input_region, rupture_file, rupture_out, dependency, dependency2):
    cmd = "module purge"
    os.system(cmd)
    cmd = "module load intel/2023.2.0 impi/2021.10.0 mkl/2023.2.0 ucx/1.15.0 oneapi/2023.2.0 bsc/1.0 fftw/3.3.10"
    os.system(cmd)
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
                            cmd =  cybershake_path + "/RuptureCodes/RupGen-api-%s/utils/num_rvs %s" % (str(gravesPitarka), os.path.join(rup_entry, f))
                            #p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
                            #print("P COMMUNICATE()", flush=True)
                            #print("COMMAND", flush=True)
                            #print(cmd, flush=True)
                            result = subprocess.run(cmd.split(), stdout=subprocess.PIPE)
                            #print("AAAAP", flush=True)
                            #print(result.stdout, flush=True)
                            #aaap = os.popen(cmd).read()
                            #print("AAAAP", flush=True)
                            #print(aaap, flush=True)
                            #print("OS SYSTEM")
                            #print(os.system(cmd), flush=True)
                            #print("P OPEN COMMUNICATE", flush=True)
                            #p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
                            #print("P COMMUNICATE",flush=True)
                            #print(p.communicate(), flush=True)
                            num_rvs = int(result.stdout.split()[-2])
                            #print("NUM RVS")
                            #print(num_rvs, flush=True)
                            for i in range(0, num_rvs):
                                query = "%s (%d, %d, %d, %d, %d, '%s')" % (query_prefix, src, rup, i, erf_id, rup_var_scenario_id, "e%d_rv%d_%d_%d.txt.variation-r%06d" % (erf_id, rup_var_scenario_id, src, rup, i))
                                cur.execute(query)
    conn.commit()
    conn.close()
    dal = DAL.SQLiteHandler(database_path, True)
    dal.generateRuptureFile(rupture_out)
    return 1

def rupVar_no_task(cybershake_path, gravesPitarka, erf_id, database_path, run_path, input_region, rupture_file, dependency, dependency2):
    cmd = "mkdir " + run_path +"/post-processing"
    os.system(cmd)
    rupture_out = run_path + "/post-processing/rupture_file_list_" + str(input_region)
    populate_rvs(cybershake_path, gravesPitarka, erf_id, database_path, run_path, input_region, rupture_file, rupture_out, dependency, dependency2)
    return rupture_out



@task(rupture_file=FILE_IN, rupture_file_list=FILE_IN, returns=1)
def make_lns(rupture_file, ruptures_erf, rupture_file_list, run_path):
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
                os.symlink("%s/%d/%d/%d_%d.txt" % (rupture_root, src, rup, src, rup), run_path + "/post-processing/" + lfn)
            except OSError as e:
                print("Exception")
                print(e)
    return 1


@constraint(computing_units="2")
@multinode(computing_nodes="3", processes_per_node="56")
@task(rupt_file=FILE_IN, returns=1)
def direct_synt(site, cybershake_path, database_path, run_path, gravesPitarka, sgt_handlers, runID, debug, max_buf, ntout, rupt_file, x_sgt_file, y_sgt_file, x_header, y_header, setup_frequency, stoch_freq, psa, rotd, dur, dtouts, out_pointsX, out_pointsY, timesamples, timeskip, highHZ, dep_pre):
    os.environ["SLURM_CPUS_ON_NODE"]="112"
    os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x3)"
    #del os.environ["SLURM_JOB_CPUS_PER_NODE"]
    os.environ["SLURM_NPROCS"]="168"
    os.environ["SLURM_NTASKS"]="168"
    os.environ["SLURM_TASKS_PER_NODE"]="56x(3)"
    # os.environ["SLURM_CPUS_PER_NODE"]="112x(3)"
    os.environ["SLURM_THREADS_PER_CORE"]="1"
    os.environ["SLURM_CPUS_PER_TASK"]="2"
    os.environ["SLURM_MEM_PER_CPU"]="2000"
    os.environ["SLURM_NTASKS_PER_NODE"]="56"
    os.environ["COMPSS_WORKING_DIR"]="/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/"
    os.chdir("/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/")
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
    f_in.write("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + "\n")
    f_in.write("%s\n" % file_python)
    f_in.flush()
    f_in.close()
    dal = DAL.SQLiteHandler(database_path, True)
    lat, lon = dal.getSiteLocation(site)
    f_python = open(file_python, "w")
    f_python.write("#!/usr/bin/env python\n\n")
    f_python.write("import sys\n\n")
    f_python.write("import os\n\n")
    f_python.write("ds_path = \"" + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(x_sgt_file) + " sgt_yfile=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(y_sgt_file) + " x_header=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(x_header) + " y_header=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filter_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_byteswap=no\"\n\n")
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
    return 1


@constraint(computing_units="2")
@multinode(computing_nodes="3", processes_per_node="56")
@task(returns=1)
def direct_synt_este_funciona_y_no_se_toca(site, cybershake_path, database_path, run_path, gravesPitarka, sgt_handlers, runID, debug, max_buf, ntout, rupt_file, x_sgt_file, y_sgt_file, x_header, y_header, setup_frequency, stoch_freq, psa, rotd, dur, dtouts, out_pointsX, out_pointsY, timesamples, timeskip, highHZ, dep_pre):
    os.environ["SLURM_CPUS_ON_NODE"]="112"
    os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x3)"
    #del os.environ["SLURM_JOB_CPUS_PER_NODE"]
    os.environ["SLURM_NPROCS"]="168"
    os.environ["SLURM_NTASKS"]="168"
    os.environ["SLURM_TASKS_PER_NODE"]="56x(3)"
    # os.environ["SLURM_CPUS_PER_NODE"]="112x(3)"
    os.environ["SLURM_THREADS_PER_CORE"]="1"
    os.environ["SLURM_CPUS_PER_TASK"]="2"
    os.environ["SLURM_MEM_PER_CPU"]="2000"
    os.environ["SLURM_NTASKS_PER_NODE"]="56"
    os.environ["COMPSS_WORKING_DIR"]="/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/"
    os.chdir("/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/")
    '''file_sh = run_path + "/post-processing/run_DS_.sh"
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
    #f_in.write("export LIBRARY_PATH=/apps/GPP/ANACONDA/2023.07/envs/libgfortran3/lib:$LIBRARY_PATH\n\n")
    f_in.write("cd " + run_path + "/post-processing/\n\n")
    f_in.write("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + "\n")
    f_in.write("%s\n" % file_python)
    f_in.flush()
    f_in.close()
    dal = DAL.SQLiteHandler(database_path, True)
    lat, lon = dal.getSiteLocation(site)
    f_python = open(file_python, "w")
    f_python.write("#!/usr/bin/env python\n\n")
    f_python.write("import sys\n\n")
    f_python.write("import os\n\n")
    f_python.write("ds_path = \"" + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(x_sgt_file) + " sgt_yfile=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(y_sgt_file) + " x_header=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(x_header) + " y_header=/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filter_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no\"\n\n")
    f_python.write("cmd = \"srun \" + str(ds_path)\n\n")
    f_python.write("os.system(cmd)\n")
    f_python.flush()
    f_python.close()
    os.system("chmod 755 " + run_path + "/post-processing/run_DS_.sh")
    os.system("chmod 755 " + run_path + "/post-processing/script.py")
    os.system("module load python/3.12.1")
    os.system("module purge")
    os.system("module load impi gcc intel fftw")
    os.system("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka))'''
    command = ["/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runDS.sh"]#/runs/SyntheticStation_179_1/runDS.sh"]
    process = subprocess.Popen(command, env=os.environ, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Obtén la salida y los errores
    stdout, stderr = process.communicate()

    # Imprime la salida
    print("Salida:", flush=True)
    print(stdout.decode(), flush=True)

    # Imprime los errores
    print("Errores:")
    print(stderr.decode())
    return 1



@constraint(computing_units="2")
@multinode(computing_nodes="3", processes_per_node="56")
@task(returns=1)
def direct_synt_copy(site, cybershake_path, database_path, run_path, gravesPitarka, sgt_handlers, runID, debug, max_buf, ntout, rupt_file, x_sgt_file, y_sgt_file, x_header, y_header, setup_frequency, stoch_freq, psa, rotd, dur, dtouts, out_pointsX, out_pointsY, timesamples, timeskip, highHZ, dep_pre):
    '''import os
    os.system("module load python/3.12.1")
    os.system("module purge")
    os.system("module load impi gcc intel fftw")
    #cmd = "cd " + run_path + "/post-processing/"
    #os.chdir(run_path + "/post-processing/")
    #print("ENVIRONMNET",flush=True)
    #print(os.environ, flush=True)
    os.environ["SLURM_CPUS_ON_NODE"]="112"
    os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x3)"
    del os.environ["SLURM_JOB_CPUS_PER_NODE"]
    os.environ["SLURM_NPROCS"]="168"
    os.environ["SLURM_NTASKS"]="168"
    os.environ["SLURM_TASKS_PER_NODE"]="56x(3)"
    os.environ["SLURM_CPUS_PER_NODE"]="112x(3)"
    os.environ["SLURM_THREADS_PER_CORE"]="1"
    os.environ["SLURM_CPUS_PER_TASK"]="2"
    os.environ["SLURM_MEM_PER_CPU"]="4000"
    os.environ["COMPSS_WORKING_DIR"]="/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/"
    print("ENVIRONMNET",flush=True)
    print(os.environ, flush=True)'''
    print("LDDDD")
    print(os.system("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + "\n"))
    #return 1
    os.system("module purge")
    os.system("module load impi gcc/13.2.0 intel fftw")
    os.system("cp " + run_path + "/" + str(x_sgt_file) + " " + run_path + "/post-processing/")
    os.system("cp " + run_path + "/" + str(y_sgt_file) + " " + run_path + "/post-processing/")
    os.system("cp " + run_path + "/" + str(x_header) + " " + run_path + "/post-processing/")
    os.system("cp " + run_path + "/" + str(y_header) + " " + run_path + "/post-processing/")
    os.environ["SLURM_CPUS_ON_NODE"]="112"
    os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x3)"
    #del os.environ["SLURM_JOB_CPUS_PER_NODE"]
    os.environ["SLURM_NPROCS"]="168"
    os.environ["SLURM_NTASKS"]="168"
    os.environ["SLURM_TASKS_PER_NODE"]="56x(3)"
    # os.environ["SLURM_CPUS_PER_NODE"]="112x(3)"
    os.environ["SLURM_THREADS_PER_CORE"]="1"
    os.environ["SLURM_CPUS_PER_TASK"]="2"
    os.environ["SLURM_MEM_PER_CPU"]="2000"
    os.environ["SLURM_NTASKS_PER_NODE"]="56"
    os.environ["COMPSS_WORKING_DIR"]="/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/"
    #dal = DAL.SQLiteHandler(database_path, True)
    #lat, lon = dal.getSiteLocation(site)
    cmd = "cd " + run_path + "/post-processing/"
    os.chdir(run_path + "/post-processing/")
    cmd = "/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runDS.sh"
    print("SECOND CMD", flush=True)
    print(cmd, flush=True)
    print("ENVIRONMEEEENT")
    print(os.environ, flush=True)
    print("LDDDD direct synth after loading correct modules")
    print(os.system("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + "\n"))
    os.system("export LD_LIBRARY_PATH=/gpfs/scratch/bsc19/bsc019756/CyberShake/software/AWP-ODC-SGT/bin:$LD_LIBRARY_PATH")
    #os.system("export LD_LIBRARY_PATH=/apps/GPP/ANACONDA/2023.07/envs/libgfortran3/lib:$LD_LIBRARY_PATH")
    print("LDDDD direct synth after loading correct modules and exporting good LD LIBRARYS")
    print(os.system("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + "\n"))
    #os.system(cmd)
    #return 1
    os.system("module purge")
    os.system("module load impi gcc/13.2.0 intel fftw")
    command = ["srun", str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + str(x_sgt_file) + " sgt_yfile=" + str(y_sgt_file) + " x_header=" + str(x_header) + " y_header=" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filer_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no"]
    #output_cmd = os.system("srun " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + str(run_path) + "/" + str(x_sgt_file) + " sgt_yfile=" + str(run_path) + "/" + str(y_sgt_file) + " x_header=" + str(run_path) + "/" + str(x_header) + " y_header=" + str(run_path) + "/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filer_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no")
    #print("OUTPUT CMD", flush=True)
    #print(output_cmd, flush=True)
    #command = os.system("srun " + )
    #process = subprocess.Popen(["srun", str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + str(run_path) + "/" + str(x_sgt_file) + " sgt_yfile=" + str(run_path) + "/" + str(y_sgt_file) + " x_header=" + str(run_path) + "/" + str(x_header) + " y_header=" + str(run_path) + "/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filer_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no"], env=os.environ, stdout=subprocess.PIPE)
    print("ENVIRON", flush=True)
    print(os.environ, flush=True)
    process = subprocess.Popen(["/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/script.py"], env=os.environ, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    #result = subprocess.run(["srun", str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + str(run_path) + "/" + str(x_sgt_file) + " sgt_yfile=" + str(run_path) + "/" + str(y_sgt_file) + " x_header=" + str(run_path) + "/" + str(x_header) + " y_header=" + str(run_path) + "/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filer_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no"], capture_output=True, text=True,shell=True, env = os.environ)
    #subprocess.run()
    '''print("SALIDA", flush=True)
    stdout, stderr = process.communicate()
    print("SALIDA", flush=True)
    print(stdout, flush=True)
    print("ERROR", flush=True)
    print(stderr, flush=True)
    #print(result.stdout, flush=True)
    #process = subprocess.Popen(command, env=os.environ, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Obtén la salida y los errores
    #stdout, stderr = process.communicate()

    # Imprime la salida
    #print("Salida:", flush=True)
    #print(stdout.decode(), flush=True)

    # Imprime los errores
    #print("Errores:")
    #print(stderr.decode())
    return stdout'''
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
    #f_in.write("export LIBRARY_PATH=/apps/GPP/ANACONDA/2023.07/envs/libgfortran3/lib:$LIBRARY_PATH\n\n")
    f_in.write("cd " + run_path + "/post-processing/\n\n")
    f_in.write("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + "\n")
    f_in.write("%s\n" % file_python)
    f_in.flush()
    f_in.close()
    dal = DAL.SQLiteHandler(database_path, True)
    lat, lon = dal.getSiteLocation(site)
    f_python = open(file_python, "w")
    f_python.write("#!/usr/bin/env python\n\n")
    f_python.write("import sys\n\n")
    f_python.write("import os\n\n")
    f_python.write("ds_path = \"" + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + str(x_sgt_file) + " sgt_yfile=" + str(y_sgt_file) + " x_header=" + str(x_header) + " y_header=" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filer_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no\"\n\n")
    #f_python.write("ds_path = str(%s)\n\n" % str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka) + " stat=" + str(site) + " slat=" + str(lat) + " slon="+str(lon) + " sgt_handlers=" + str(sgt_handlers) + " run_id=" + str(runID) + " debug=" + str(debug) + " max_bug_mb=" + str(max_buf) + " rupture_spacing=uniform ntout=" + str(ntout) + " rup_list_file=" + str(rupt_file) + " sgt_xfile=" + str(run_path) + "/" + str(x_sgt_file) + " sgt_yfile=" + str(run_path) + "/" + str(y_sgt_file) + " x_header=" + str(run_path) + "/" + str(x_header) + " y_header=" + str(run_path) + "/" + str(y_header) + " det_max_freq=" + str(setup_frequency) + " stoch_max_freq=" + str(stoch_freq) + " run_psa=" + str(psa) + " run_rotd=" + str(rotd) + " run_durations=" + str(dur) + " dtout=" + str(dtouts) + " simulation_out_pointsX=" + str(out_pointsX) + " simulation_out_pointsY=" + str(out_pointsY) + " simulation_out_timesamples=" + str(timesamples) + " simulation_out_timeskip=" + str(timeskip) + " surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filer_highHZ=" + str(highHZ) + " surfseis_rspectra_apply_bytesawp=no")
    f_python.write("cmd = \"srun \" + str(ds_path)\n\n")
    f_python.write("os.system(cmd)\n")
    f_python.flush()
    f_python.close()
    os.system("chmod 755 " + run_path + "/post-processing/run_DS_.sh")
    os.system("chmod 755 " + run_path + "/post-processing/script.py")
    os.system("module load python/3.12.1")
    os.system("module purge")
    os.system("module load impi gcc intel fftw")
    os.system("ldd " + str(cybershake_path) + "/DirectSynth/bin/direct_synth_v" + str(gravesPitarka))
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
    return 1#stderr.decode()

def runDS_no_task(cybershake_path, database_path, rupture_file, rupture_file_list, run_path, site, runID, input_region, input_erf_ruptures, input_erf_id, gravesPitarka, compute_setup_frequency, config):
    dependency = make_lns(rupture_file, "Ruptures_erf34" + "/", rupture_file_list, run_path)
    direct_synt("SS_179", cybershake_path, database_path, run_path, gravesPitarka, 84, runID, 1, 512, 3000, rupture_file_list,   "SS_179_fx_" + str(runID) + ".sgt", "SS_179_fy_" + str(runID) + ".sgt", "SS_179_fx_" + str(runID) + ".sgthead", "SS_179_fy_" + str(runID) +".sgthead", config["compute"]["setup"]["frequency"], -1.0, 1, 1, 1, 0.1, 2, 1, 3000, 0.1, 5.0, dependency)

@constraint(computing_units="112")
@multinode(computing_nodes="1")
@task(returns=1)
def rupVar(cybershake_path, gravesPitarka, erf_id, database_path, run_path, input_region, dependency, dependency2):
    cmd = (cybershake_path + "/populate_rvs.py " + str(gravesPitarka) + " " + str(erf_id) + " 1 " + str(database_path) + " " + run_path + "/ruptures")
    os.system(cmd)
    cmd = "mkdir " + run_path +"/post-processing"
    os.system(cmd)
    cmd = "cd " + run_path +"/post-processing"
    os.system(cmd)
    dal = DAL.SQLiteHandler(database_path, True)
    dal.generateRuptureFile(run_path+"/post-processing/rupture_file_list_" + str(input_region))
    return 1

@constraint(computing_units="2")
@multinode(computing_nodes="3", processes_per_node="56")
@task(returns=1)
def runDS(cybershake_path, run_path, site, runID, input_region, input_erf_ruptures, input_erf_id, database_path, gravesPitarka, compute_setup_frequency, dependency):
    import os
    #print("ENVIRONMNET",flush=True)
    #print(os.environ, flush=True)
    #os.environ["SLURM_CPUS_ON_NODE"]="112"
    #os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x3)"
    #os.environ["SLURM_NPROCS"]="168"
    #os.environ["SLURM_MEM_PER_CPU"]="2000"
    cmd = "cd " + run_path + "/post-processing/"
    os.chdir(run_path + "/post-processing/")
    #print("LDD AL BINARIO", flush=True)
    #os.system("ldd /gpfs/scratch/bsc19/bsc019756/CyberShake/software/DirectSynth/bin/direct_synth_v5.4.2")
    cmd = "ln -s " + run_path + "/" + site + "_fx_" + str(runID) + ".sgt " +run_path + "/"  + "/post-processing/"+ site+ "_fx_" + str(runID) + ".sgt"
    #cmd = "ln -s "+ run_path + "/post-processing/ " + run_path+"/"+ site + "_fx_" + str(runID) + ".sgt"
    os.system(cmd)
    cmd = "ln -s " + run_path + "/" + site + "_fy_" + str(runID) + ".sgt " +run_path + "/"  + "/post-processing/"+ site+ "_fy_" + str(runID) + ".sgt"
    #cmd = "ln -s "+ run_path + "/post-processing/ " +run_path + "/" + site + "_fy_" + str(runID) + ".sgt"
    os.system(cmd)
    cmd = "ln -s " + run_path + "/" + site + "_fx_" + str(runID) + ".sgthead " + run_path + "/post-processing/"+ site + "_fx_" + str(runID) + ".sgthead"
    #cmd = "ln -s " + run_path + "/post-processing/ "+run_path + "/" + site + "_fx_" + str(runID) + ".sgthead"
    os.system(cmd)
    cmd = "ln -s " + run_path + "/" +site + "_fy_" + str(runID) + ".sgthead " + run_path + "/post-processing/"+ site + "_fy_" + str(runID) + ".sgthead"
    #cmd = "ln -s " + run_path + "/post-processing/ "+run_path + "/" + site + "_fy_" + str(runID) + ".sgthead"
    os.system(cmd)
    cmd = (cybershake_path + "/make_lns.py " + run_path + "/post-processing/rupture_file_list_" + input_region + " " + input_erf_ruptures + "Ruptures_erf"
                   + str(input_erf_id) + "/")
    #print("COMMAND", flush=True)
    #print(cmd, flush=True)
    os.system(cmd)
    dal = DAL.SQLiteHandler(database_path, True)
    lat, lon = dal.getSiteLocation(site)
    cmd = "/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runDS.sh"
    print("SECOND CMD", flush=True)
    print(cmd, flush=True)
    os.system(cmd)
    return 1
    #import os
    #cmd = "cd " + run_path + "/post-processing/"
    #os.chdir(run_path + "/post-processing/")
    #print("ENVIRONMNET",flush=True)
    #print(os.environ, flush=True)
    #os.environ["SLURM_CPUS_ON_NODE"]="112"
    #os.environ["SLURM_JOB_CPUS_PER_NODE"]="112(x3)"
    #del os.environ["SLURM_JOB_CPUS_PER_NODE"]
    #os.environ["SLURM_NPROCS"]="168"
    #os.environ["SLURM_NTASKS"]="168"
    #os.environ["SLURM_TASKS_PER_NODE"]="56x(3)"
    #os.environ["SLURM_CPUS_PER_NODE"]="112x(3)"
    #os.environ["SLURM_THREADS_PER_CORE"]="1"
    #os.environ["SLURM_CPUS_PER_TASK"]="2"
    #os.environ["SLURM_MEM_PER_CPU"]="4000"
    #os.environ["COMPSS_WORKING_DIR"]="/gpfs/scratch/bsc19/bsc019756/Beticas_Dataset/beticas_f627_compss/runs/SyntheticStation_179_1/post-processing/"
    

@constraint(computing_units="2")
@mpi(binary='/gpfs/scratch/bsc19/bsc019756/CyberShake/software/DirectSynth/bin/direct_synth_v5.4.2', working_dir="{{work_dir}}", args="stat={{site}} slat={{lat}} slon={{lon}} sgt_handlers={{sgt_handlers}} run_id={{runID}} debug={{debug}} max_buf_mb={{max_buf}} rupture_spacing=uniform ntout={{ntout}} rup_list_file={{rupt_file}} sgt_xfile={{x_sgt_file}} sgt_yfile={{y_sgt_file}} x_header={{xheader}} y_header={{yheader}} det_max_freq={{setup_frequency}} stoch_max_freq={{stoch_freq}} run_psa={{psa}} run_rotd={{rotd}} run_durations={{dur}} dtout={{dtouts}} simulation_out_pointsX={{out_pointsX}} simulation_out_pointsY={{out_pointsY}} simulation_out_timesamples={{timesamples}} simulation_out_timeskip={{timeskip}} surfseis_rspectra_seismogram_units=cmpersec surfseis_rspectra_output_units=cmpersec2 surfseis_rspectra_output_type=aa surfseis_rspectra_period=all surfseis_rspectra_apply_filter_highHZ={{highHZ}} surfseis_rspectra_apply_byteswap=no", processes_per_node=56, processes=168, runner='srun')#, computing_nodes="${ComputingNodesAWP}")
@task(returns=1)
def runDSMPI(work_dir, site, lat, lon, sgt_handlers, runID, debug, max_buf, ntout, rupt_file, x_sgt_file, y_sgt_file, xheader, yheader, setup_frequency, stoch_freq, psa, rotd, dur, dtouts, out_pointsX, out_pointsY, timesamples, timeskip, highHZ, dep_pre):
    pass

