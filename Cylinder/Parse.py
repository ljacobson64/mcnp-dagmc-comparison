import itertools
import os

# Determine filename
def determine_filename(version,H,D,mfp,tol):
    
    H_str   = ('%4.0f' % (H  )).replace(' ','_')
    D_str   = ('%4.0f' % (D  )).replace(' ','_')
    mfp_str = ('%6.2f' % (mfp)).replace(' ','_')
    tol_str = ('%6.2f' % (tol)).replace(' ','_')
    
    fname = 'zCyl_%s_%s_%s_%s_%s' % (version,H_str,D_str,mfp_str,tol_str)
    return fname

# Parse output files for CTM and NPS
def parse_output_file(fname,version,H,D,mfp,tol):
    
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
        print >> writer, '|---------|-------|-------|-------|---------|---------|---------|-------------|'
        print >> writer, '| version |   H   |   D   |  mfp  |   tol   |  ctm1   |  ctm2   |     nps     |'
        print >> writer, '|---------|-------|-------|-------|---------|---------|---------|-------------|'
        
    print >> writer, '|   %s   | %5u | %5u | %5.0e | %7.2f | %7.2f | %7.2f | %11u |' % (version,H,D,mfp,tol,ctm1,ctm2,nps)
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params_str = reader.readlines()

global direc
direc    =                    params_str[0].split()[0 ]                         # directory to place files
resdir   =                    params_str[1].split()[0 ]                         # directory to place results
versions =                    params_str[2].split()[1:]                         # list of MCNP versions
H_vals   = [  int(i) for i in params_str[3].split()[1:]]                        # list of values of N
D_vals   = [float(i) for i in params_str[4].split()[1:]]                        # list of values of F
mfps     = [float(i) for i in params_str[5].split()[1:]]                        # list of number of mean free paths in 1 m
tols     = [float(i) for i in params_str[6].split()[1:]]                        # list of faceting tolerances
ctme     =  float(            params_str[7].split()[1 ])                        # computer time

reader.close()

if os.path.isfile(resdir+'Cube_results.txt') == True:
    os.remove(resdir+'Cube_results.txt')

# Parse output
for params in list(itertools.product(versions,H_vals,D_vals,mfps,tols)):
    
    version = params[0]                                                         # loop for all MCNP versions
    H       = params[1]                                                         # loop for all values of H
    D       = params[2]                                                         # loop for all values of D
    mfp     = params[3]                                                         # loop for all values of mfp
    tol     = params[4]                                                         # loop for all values of tol
    
    fname = determine_filename(version,H,D,mfp,tol)                             # base filename (no extension)
    
    parse_output_file(fname,version,H,D,mfp,tol)                                # parse output file for CTM and NPS
    print fname

# Display output in command window
f = open(resdir+'Cube_results.txt','r')
print f.read()
f.close()
