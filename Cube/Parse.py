import itertools
import os

# Determine filename
def determine_filename(version,geom,N,F,mfp,tol):
    
    N_str   = ('%5u'   % (N  )).replace(' ','_')
    F_str   = ('%2.0f' % (F  )).replace(' ','_')
    mfp_str = ('%6.2f' % (mfp)).replace(' ','_')
    tol_str = ('%3.0f' % (tol)).replace(' ','_')
    
    fname = 'zCube_%s_%s_%s_%s_%s_%s' % \
            (version,geom,N_str,F_str,mfp_str,tol_str)
    return fname

# Parse output files for CTM and NPS
def parse_output_file(fname,version,geom,N,F,mfp):
    
    fname_io = '%s.io'   % (fname)
    fname_out = '%s.out' % (fname)
    
    # Make sure output file is valid
    if (os.path.isfile(direc+fname_io) == False):
        nps = -1
        ctm1 = -1
        ctm2 = -1
    else:
        
        # Grab last 512 bytes
        reader = open(direc+fname_io,'r')
        reader.seek(0,2)
        reader.seek(max(reader.tell()-512,0),0)
        lines = reader.read()
        reader.close()
        
        # Find NPS
        nps_loc = str.find(lines,'nps')
        if nps_loc > 0:
            nps_str = lines[nps_loc+5:nps_loc+20]
            nps = int(nps_str)
        else:
            nps = -1
        
        # Find CTM 1
        ctm1_loc = str.find(lines,'ctm')
        if ctm1_loc > 0:
            ctm1_str = lines[ctm1_loc+5:ctm1_loc+18]
            ctm1 = float(ctm1_str)
        else:
            ctm1 = -1
        
        # Find CTM 2
        ctm2_loc = str.find(lines,'time =')
        if ctm2_loc > 0:
            ctm2_str = lines[ctm2_loc+6:ctm2_loc+15]
            ctm2 = float(ctm2_str)
        else:
            ctm2 = -1
    
    # Write NPS to text file
    writer = open(resdir+'Cube_results.txt','a')
    
    if writer.tell() == 0:
        print >> writer, '|---------|------|-------|---------|-------|---------|---------|------------|'
        print >> writer, '| version | geom |   N   |    F    |  mfp  |  ctm1   |  ctm2   |    nps     |'
        print >> writer, '|---------|------|-------|---------|-------|---------|---------|------------|'
        
    print >> writer, '|   %s   |  %s  | %5u | %7.3f | %5.0e | %7.2f | %7.2f | %10u |' % (version,geom,N,F,mfp,ctm1,ctm2,nps)
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params_str = reader.readlines()

global direc
direc    =                    params_str[0].split()[0 ]                         # directory to place files
resdir   =                    params_str[1].split()[0 ]                         # directory to place results
versions =                    params_str[2].split()[1:]                         # list of MCNP versions
geoms    =                    params_str[3].split()[1:]                         # list of geometry configurations
N_vals   = [  int(i) for i in params_str[4].split()[1:]]                        # list of values of N
F_vals   = [float(i) for i in params_str[5].split()[1:]]                        # list of values of F
mfps     = [float(i) for i in params_str[6].split()[1:]]                        # list of number of mean free paths in 1 m
ctme     =  float(            params_str[7].split()[1 ])                        # computer time
mfp_in   =  float(            params_str[8].split()[1 ])                        # mfp for pre-existing LCAD file
toltype  =                    params_str[9].split()[1 ]  # type of faceting tolerance to use for DAG
tols     = [float(i) for i in params_str[9].split()[2:]] # list of faceting tolerances to use for DAG

reader.close()

if os.path.isfile(resdir+'Cube_results.txt') == True:
    os.remove(resdir+'Cube_results.txt')

# Parse output
for params in list(itertools.product(versions,geoms,N_vals,F_vals,mfps,tols)):
    
    version = params[0]                                                         # loop for all MCNP versions
    geom    = params[1]                                                         # loop for all geometry configurations
    N       = params[2]                                                         # loop for all values of N
    F       = params[3]                                                         # loop for all values of F
    mfp     = params[4]                                                         # loop for all values of mfp
    tol     = params[5]                                                         # loop for all values of tol
    
    fname = determine_filename(version,geom,N,F,mfp,tol)                            # base filename (no extension)
    
    parse_output_file(fname,version,geom,N,F,mfp)                               # parse output file for CTM and NPS
    print fname

# Display output in command window
f = open(resdir+'Cube_results.txt','r')
print f.read()
f.close()
