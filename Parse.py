import itertools
import os

# Parse output files for CTM and NPS
def parse_output_file(fname,version,geom,N,F,mfp,fname_res):
    
    fname_io = '%s.io' % (fname)
    fname_out = '%s.out' % (fname)
    if os.path.isfile(direc+fname_io) == False:
        return
    # if os.path.getsize(direc+fname_out) == 0:
    #     return
    
    reader = open(direc+fname_io,'r')
    
    # Grab last 512 bytes
    reader.seek(0,2)
    reader.seek(max(reader.tell()-512,0),0)
    lines = reader.read()
    
    reader.close()
    
    # Find NPS
    nps_loc = str.find(lines,'nps')
    nps_str = lines[nps_loc+5:nps_loc+20]
    nps = int(nps_str)
    
    # Find CTM 1
    ctm1_loc = str.find(lines,'ctm')
    ctm1_str = lines[ctm1_loc+5:ctm1_loc+18]
    ctm1 = float(ctm1_str)
    
    # Find CTM 2
    ctm2_loc = str.find(lines,'time =')
    ctm2_str = lines[ctm2_loc+6:ctm2_loc+15]
    ctm2 = float(ctm2_str)

    writer = open(fname_res,'a')
    
    if writer.tell() == 0:
        print >> writer, '|---------|------|-------|---------|----------|--------|--------|-----------|'
        print >> writer, '| version | geom |   N   |    F    |   mfp    |  ctm1  |  ctm2  |    nps    |'
        print >> writer, '|---------|------|-------|---------|----------|--------|--------|-----------|'
        
    # Write NPS to text file
    print >> writer, '|   %s   |  %s  | %5u | %7.3f | %8.3f | %6.2f | %6.2f | %9u |' % (version,geom,N,F,mfp,ctm1,ctm2,nps)
    
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params_str = reader.readlines()

global direc
direc      =                    params_str[0].split()[0 ]                        # directory to place files
versions   =                    params_str[1].split()[1:]                        # list of MCNP versions
geoms      =                    params_str[2].split()[1:]                        # list of geometry configurations
N_vals     = [  int(i) for i in params_str[3].split()[1:]]                       # list of values of N
F_vals     = [float(i) for i in params_str[4].split()[1:]]                       # list of values of F
mfps       = [float(i) for i in params_str[5].split()[1:]]                       # list of number of mean free paths in 1 m
ctme       =  float(            params_str[6].split()[1 ])                       # computer time
mfp_in     =  float(            params_str[7].split()[1 ])                       # mfp for pre-existing LCAD file

reader.close()

if os.path.isfile('Cube_results.txt') == True:
    os.remove('Cube_results.txt')

# Parse output
for params in list(itertools.product(versions,geoms,N_vals,F_vals,mfps)):
    version = params[0]                                                         # loop for all MCNP versions
    geom    = params[1]                                                         # loop for all geometry configurations
    N       = params[2]                                                         # loop for all values of N
    F       = params[3]                                                         # loop for all values of F
    mfp     = params[4]                                                         # loop for all values of mfp
    
    fname = 'zCube_%s_%s_%u_%.0f_%.2f' % (version,geom,N,F,mfp)                 # base filename (no extension)
    
    parse_output_file(fname,version,geom,N,F,mfp,'Cube_results.txt')            # parse output file for CTM and NPS
    print fname

# Display output in command window
f = open('Cube_results.txt','r')
print f.read()
f.close()
