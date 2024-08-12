import sqlite3
import os
from modules.cybershake_workflows.steps.utils import getProperty
from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *


@constraint(computing_units="${ComputingUnitsSgtGrid}")
@task(sgtcord_file_name=FILE_OUT, fault_list_file_name=FILE_IN, radius_file_name=FILE_IN, cfg_file=FILE_IN)
def gen_sgt_grid(sgtcord_file_name, cybershake_path, site, ns, src, mlon, mlat, mrot, fault_list_file_name, radius_file_name, frequency, spacing, cfg_file):
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

    command = '%s %s/PreSgt/bin/gen_sgtgrid nx=%d ny=%d nz=%d h=%f xsrc=%d ysrc=%d ixmin=%d ixmax=%d iymin=%d iymax=%d izstart=%d izmax=%d radiusfile=%s outfile=%s modellon=%f modellat=%f modelrot=%f faultlist=%s' % (MPI_CMD, cybershake_path, int(ns[0]), int(ns[1]), int(ns[2]), HH, int(src[0]), int(src[1]), IX_MIN, int(IX_MAX), IY_MIN, int(IY_MAX), IZ_START, int(IZ_MAX), radius_file_name, sgtcord_file_name, float(mlon), float(mlat), float(mrot), fault_list_file_name)
    returnCode = os.system(command)
    if returnCode!=0:
        sys.exit((returnCode >> 8) & 0xFF)

@constraint(computing_units="${ComputingUnits}")
@task(radius_file=FILE_OUT)
def gen_radius_file(radius_file):
    '''Duplicates the functionality of part of gen_sgtgrid.csh:  it writes the adaptive mesh info to a radius file as part of generating the cordfile.'''
    #magic constants for the adaptive mesh
    RLEV = [10.0, 50.0, 100.0, 1000.0]
    RINC = [10, 15, 25, 50]
    ZLEV = [0.2, 5.0, 24.0, 60.0]
    ZINC = [1, 5, 10, 25]

    output = open(radius_file, 'w')
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



def fault_list_script(site, erf_id, rupture_path, fault_list_file, db_file, rsqsim, *args):
    if (site==None or erf_id==None or rupture_path==None or fault_list_file==None):
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
    conn = sqlite3.connect(server)
    cur = conn.cursor()
    query = 'select R.Source_ID, R.Rupture_ID from CyberShake_Site_Ruptures R, CyberShake_Sites S where S.CS_Short_Name="%s" and R.CS_Site_ID=S.CS_Site_ID and R.ERF_ID=%d order by R.Source_ID, R.Rupture_ID' % (site, int(erf_id))
    cur.execute(query)
    res = cur.fetchall()
    if len(res)==0:
        print("No ruptures found for site %s." % site)
        sys.exit(2)

    fp_out = open(fault_list_file, "w")
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
    return 0


@constraint(computing_units="${ComputingUnits}")
@task(faultlistfile=FILE_OUT)
def gen_fault_list(faultlistfile, site, erf_id, rupt_geometries, rsqsim, db_file):
    '''Copies the functionality of gen_faultlist.csh:  it serves as a wrapper to CreateFaultList.java, which queries the database to generate a list of ruptures which are applicable for the given site.'''
    rsqsim_str = ""
    if rsqsim:
        rsqsim_str = "-rsqsim"
    returnCode = fault_list_script(site, erf_id, rupt_geometries, faultlistfile, db_file, rsqsim_str=="-rsqsim")
    if returnCode!=0:
        sys.exit((returnCode >> 8) & 0xFF)


@constraint(computing_units="${ComputingUnits}")
@task(fdloc_file=FILE_OUT, cordfileName=FILE_IN, returns=1)
def gen_fdloc(fdloc_file, site, mlon, mlat, cordfileName):
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
    fdlocfile = open(fdloc_file, 'w')
    fdlocfile.write('%s %s\n' % (minLoc[0], minLoc[1]))
    fdlocfile.flush()
    fdlocfile.close()
    return [int(minLoc[0]), int(minLoc[1])]

@constraint(computing_units="${ComputingUnits}")
@task(returns=2)
def get_site_coords(site, db_file):
    connection = sqlite3.connect(db_file)
    cursor = connection.cursor()
    sql_string = 'select CS_Site_Lat, CS_Site_Lon from CyberShake_Sites where CS_Short_Name="%s"' % site
    cursor.execute(sql_string)
    siteCoords = cursor.fetchone()
    return float(siteCoords[0]), float(siteCoords[1])


def step_preSGT(config, id, runID, site, run_path, erf_id, modelbox, gridout, cord_file_name, db_file, setup_spacing, frequency, cfg_file_input, fdloc_file_name, fault_list_file_name, radius_file_name, sgtcord_file_name, rsqsim=True):
    spacing = float(setup_spacing)
    if frequency is None:
        frequency = 0.5
    PATH_TO_RUPTURE_GEOMETRIES = "%s/Ruptures_erf%s" % (getProperty('RUPTURE_ROOT', cfg_file_input), erf_id)
    input = open(modelbox)
    modelboxContents = input.readlines()
    input.close()

    siteLat, siteLon = get_site_coords(site, db_file)

    modelTokens = modelboxContents[4].split()
    mlon = (float)(modelTokens[1])
    mlat = (float)(modelTokens[3])
    mrot = (float)(modelTokens[5])

    src = gen_fdloc(fdloc_file_name, site, siteLon, siteLat, cord_file_name)
    gen_fault_list(fault_list_file_name, site, erf_id, PATH_TO_RUPTURE_GEOMETRIES, rsqsim, db_file)

    gen_radius_file(radius_file_name)

    input = open(gridout)
    gridoutContents = input.readlines();
    input.close()
    ns = []
    ns.append(int((gridoutContents[1].split("="))[1]))
    ns.append(int((gridoutContents[1+ns[0]+2].split("="))[1]))
    ns.append(int((gridoutContents[1+ns[0]+2+ns[1]+2].split("="))[1]))
    gen_sgt_grid(sgtcord_file_name, config["input"]["cyberShake"]["path"], site, ns, src, mlon, mlat, mrot, fault_list_file_name, radius_file_name, frequency, spacing, cfg_file_input)
    return sgtcord_file_name, fault_list_file_name, radius_file_name, fdloc_file_name

