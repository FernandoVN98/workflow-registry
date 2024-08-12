import os
import sys

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
