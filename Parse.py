import os

# Parse output files for CTM and NPS
def parse_output_file(fname,N,F,geom,fname_res):
    
    fname_io = '%s.io' % (fname)
    fname_out = '%s.out' % (fname)
    if os.path.isfile(direc+fname_io) == False:
        return
    if os.path.getsize(direc+fname_out) == 0:
    	return
        
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
        print >> writer, '|------|-------|-----|--------|--------|-----------|'
        print >> writer, '| type |   N   |  F  |  ctm1  |  ctm2  |    nps    |'
        print >> writer, '|------|-------|-----|--------|--------|-----------|'
        
    # Write NPS to text file
    print >> writer, '|  %s  | %5u | %3.0f | %6.2f | %6.2f | %9u |' % (geom,N,F,ctm1,ctm2,nps)
    
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params = reader.readlines()

global direc
direc      =                    params[0].split()[0 ]                            # directory to place files
N_vals     = [  int(i) for i in params[1].split()[1:]]                           # list of values of N
F_vals     = [float(i) for i in params[2].split()[1:]]                           # list of values of F
geoms      =                    params[3].split()[1:]                            # list of geometry configurations
ctme       =  float(            params[4].split()[1 ])                           # computer time
versions   =                    params[5].split()[1:]                            # list of MCNP versions

reader.close()

max_N_2s = 40000
max_N_2j =  1000

if os.path.isfile('Cube_results.txt') == True:
    os.remove('Cube_results.txt')

# Parse output
for version in versions:                                                        # loop for all MCNP versions
    for geom in geoms:                                                          # loop for all geometry configurations
        for N in N_vals:                                                        # loop for all values of N
            for F in F_vals:                                                    # loop for all values of F
                
                fname = 'zCube_%s_%s_%u_%.0f' % (version,geom,N,F)              # base filename (no extension)
                
                # Parse MCNP output file
                parse_output_file(fname,N,F,geom,'Cube_results.txt')            # parse output file for CTM and NPS
                
                print fname

# Display output in command window
f = open('Cube_results.txt','r')
print f.read()
f.close()
