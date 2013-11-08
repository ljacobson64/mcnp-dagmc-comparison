# Parse output files for CTM and NPS
def parse_output_file(filename,N,F,geom_type,filename_res):
    
    filename_o = '%s.io' % (filename)
    reader = open(filename_o,'r')
    
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

    writer = open(filename_res,'a')
    
    # Write NPS to text file
    if writer.tell() == 0:
        print >> writer, '|------|-------|-----|-------|-------|----------|'
        print >> writer, '| type |   N   |  F  | ctm1  | ctm2  |   nps    |'
        print >> writer, '|------|-------|-----|-------|-------|----------|'
    print >> writer, '|  %s  | %5u | %3.0f | %5.2f | %5.2f | %8u |' % (geom_type,N,F,ctm1,ctm2,nps)
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params = reader.readlines()

N_vals     = [  int(i) for i in params[0].split()[1:]]                         # list of values of N
F_vals     = [float(i) for i in params[1].split()[1:]]                         # list of values of F
geom_types =                    params[2].split()[1:]                          # list of geometry configurations
ctme       =  float(            params[3].split()[1 ])                         # computer time

reader.close()

# Parse output
for geom_type in geom_types:                                                   # loop for all geometry configurations
    for N in N_vals:                                                           # loop for all values of N
        for F in F_vals:                                                       # loop for all values of F
            
            if geom_type == '2s':
                valid_N = N <= max_N_2s
            if geom_type == '2j':
                valid_N = N <= max_N_2j
            
            if valid_N:
                
                filename = 'Cube_%s_%u_%.0f' % (geom_type,N,F)                     # Base filename (no extension)
                
                # Parse MCNP output file
                parse_output_file(filename,N,F,geom_type,'Cube_results.txt')       # parse output file for CTM and NPS

# Display output in command window
f = open('Cube_results.txt','r')
print f.read()
f.close()
