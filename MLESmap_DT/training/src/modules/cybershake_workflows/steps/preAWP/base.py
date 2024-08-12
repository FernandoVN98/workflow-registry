import os
import sys
from pycompss.api.task import task
from pycompss.api.constraint import constraint
from pycompss.api.parameter import *


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as ex:
        if ex.errno==errno.EEXIST and os.path.isdir(path):
            pass
        else: raise


@constraint(computing_units="${ComputingUnits}")
@task(gridout=FILE_IN, awp_media=FILE_INOUT)
def build_media(cybershake_path, site, gridout, rwg_vel_prefix, awp_media):
    fp_in = open(gridout, "r")
    data = fp_in.readlines()
    fp_in.close()
    nx = int((data[1].split("="))[1])
    ny = int((data[1+nx+2].split("="))[1])
    nz = int((data[1+nx+2+ny+2].split("="))[1])
    cmd = "%s/SgtHead/bin/reformat_velocity %d %d %d %s %s" % (cybershake_path, nx, ny, nz, rwg_vel_prefix, awp_media)
    exitcode = os.system(cmd)
    fp_in = open(awp_media, "rb")
    data = fp_in.read()
    fp_in.close()
    fp_out = open(awp_media, "wb")
    fp_out.write(data)
    fp_out.flush()
    fp_out.close()


@constraint(computing_units="${ComputingUnits}")
@task(cordfile=FILE_IN, awp_cordfile=FILE_OUT)
def build_cordfile(cordfile, awp_cordfile, max_depth_index):
    fp_in = open(cordfile, "r")
    fp_out = open(awp_cordfile, "w")
    max_depth_index = int(max_depth_index)

    data = fp_in.readlines()
    fp_in.close()

    points = dict()

    num_pts_in = int(data[4])

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
    p_list = sorted(points)
    fp_out.write("%d\n" % len(p_list))
    for entry in p_list:
        #Print as (y, x, z)
        x = entry/(100000*10000)
        y = (entry/10000) % 100000
        z = entry % 10000
        #1 already added, don't need it here
        fp_out.write("%d %d %d\n" % (y, x, z))
    fp_out.flush()
    fp_out.close()


@constraint(computing_units="${ComputingUnits}")
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


@constraint(computing_units="${ComputingUnits}")
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
    param["WRITE_STEP"] = int(param["READ_STEP"]/int(param["NTISKP_SGT"]))
    param["WRITE_STEP2"] = int(param["WRITE_STEP"])
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
        fp_out.write("'%s' INVEL\n\n" % (velocity_mesh))
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

        if c == "x":
            build_IN3D(out_IN3D_x, cybershake_path+"/AWP-ODC-SGT/utils/data/IN3D.ref", site, gridout, c, frequency, procs, run_path, velocity_mesh, rwg_vel_prefix, spacing=spacing)
        elif c == "y":
            build_IN3D(out_IN3D_y, cybershake_path+"/AWP-ODC-SGT/utils/data/IN3D.ref", site, gridout, c, frequency, procs, run_path, velocity_mesh, rwg_vel_prefix, spacing=spacing)
         
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
        build_media(cybershake_path, site, gridout, rwg_vel_prefix, awp_media)
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

def step_preAWP(run_path, cybershake_path, model_path, site_compound, launching_script, site, grid_out, frequency, px, py, pz, source_frequency, runID, velocity_mesh, cordfile, fdloc):
    os.symlink(model_path, str(run_path) + "/awp." + site + ".media")
    os.chdir(run_path)
    spacing=None
    out_IN3D_x, out_IN3D_y, out_src_x, out_src_y, out_cordfile, out_media = build_awp_inputs(run_path, cybershake_path, site, grid_out, fdloc, cordfile, px, py, pz, velocity_mesh, frequency, source_frequency, spacing)
    return out_IN3D_x, out_IN3D_y, out_src_x, out_src_y, out_cordfile, out_media

