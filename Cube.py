from subprocess import call
import os

# Write input file; 2 dimensions with separate cells
def write_mcnp_input_2s(filename,N,F,ctme):
	
    writer = open(filename,'w')
    
    # Write title cards
    print >> writer, 'Small geometric features test'
    print >> writer, 'c Cube geometry with separate cells in 2 dimensions'
    print >> writer, 'c N = %u; F = %4.1f' % (N,F)
    
    # Write cell cards
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0        -09000'
    print >> writer, '   12 0  09999'
    print >> writer, '   13 0  09000 -09999        -01000'
    print >> writer, '   14 0  09000 -09999  04999'
    print >> writer, '   15 0  09000 -09999  01000 -04999        -05000'
    print >> writer, '   16 0  09000 -09999  01000 -04999  08999'
    print >> writer, '    1 0  09000 -09999  %u -04999  %u -08999' % (10000+N,50000+N)
    print >> writer, '    2 0  09000 -09999  01000 -%u  %u -08999' % (10000+N,50000+N)
    print >> writer, '    3 0  09000 -09999  %u -04999  05000 -%u' % (10000+N,50000+N)
    for j in range(1,N):
        print >> writer, '%u 0  09000 -09999  %u -%u  %u -%u' % (10000+j,10000+j,10000+N,50000+N-j,50000+N-j+1)
        print >> writer, '%u 0  09000 -09999  01000 -%u  %u -%u' % (50000+j,10000+j,50000+N-j,50000+N-j+1)
    print >> writer, '%u 0  09000 -09999  01000 -%u  05000 -50001' % (50000+N,10000+N)
    print >> writer, ''
    
    # Write surface cards
    print >> writer, 'c SURFACE CARDS'
    print >> writer, '01000 PX 0'
    print >> writer, '05000 PY 0'
    print >> writer, '09000 PZ 0'
    print >> writer, '04999 PX 100'
    print >> writer, '08999 PY 100'
    print >> writer, '09999 PZ 100'
    for j in range(1,N+1):
        print >> writer, '%u PX %11.8f' % (10000+j,F*j/N)
        print >> writer, '%u PY %11.8f' % (50000+j,F*j/N)
    print >> writer, ''
    
    # Write data cards
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'IMP:N 0 5R 1 %uR' % (2*N+1)
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % (ctme)
    
    writer.close()

# Write input file; 2 dimensions with joined cells
def write_mcnp_input_2j(filename,N,F,ctme):
	
    writer = open(filename,'w')
    
    # Write title cards
    print >> writer, 'Small geometric features test'
    print >> writer, 'c Cube geometry with joined cells in 2 dimensions'
    print >> writer, 'c N = %u; F = %4.1f' % (N,F)
    
    # Write cell cards
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0  09999:-09000: 04999:-01000: 08999:-05000'
    print >> writer, '    1 0  09000 -09999 -04999 -08999 ('
    print >> writer, '         01000  %u:' % (50000+N)
    for j in range(1,N):
        print >> writer, '         %u  %u:' % (10000+j,50000+N-j)
    print >> writer, '         %u  05000)' % (10000+N)
    print >> writer, '    2 0  09000 -09999  01000  05000 ('
    for j in range(1,N):
        print >> writer, '        -%u -%u:' % (10000+j,50000+N-j+1)
    print >> writer, '        -%u -%u)' % (10000+N,50000+1)
    print >> writer, ''
    
    # Write surface cards
    print >> writer, 'c SURFACE CARDS'
    print >> writer, '01000 PX 0'
    print >> writer, '05000 PY 0'
    print >> writer, '09000 PZ 0'
    print >> writer, '04999 PX 100'
    print >> writer, '08999 PY 100'
    print >> writer, '09999 PZ 100'
    for j in range(1,N+1):
        print >> writer, '%u PX %11.8f' % (10000+j,F*j/N)
        print >> writer, '%u PY %11.8f' % (50000+j,F*j/N)
    print >> writer, ''

    
    # Write data cards
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'IMP:N 0 1 1'
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % ctme
    
    writer.close()

# Parse output files for CTM and NPS
def parse_output_file(filename,N,F,geom_type,filename_res):
    
    reader = open(filename+'o','r')
    
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

# Write input files, run MCNP, parse output
for geom_type in geom_types:                                                   # loop for all geometry configurations
    for N in N_vals:                                                           # loop for all values of N
        for F in F_vals:                                                       # loop for all values of F
            
            # Write MCNP input file
            filename = 'Cube_%s_%u_%.0f.i' % (geom_type,N,F)                   # input filename
            if geom_type=='2s':                                                # 2 dimensions with separate cells
                write_mcnp_input_2s(filename,N,F,ctme)                         # write input file
            if geom_type=='2j':                                                # 2 dimensions with joined cells
                write_mcnp_input_2j(filename,N,F,ctme)                         # write input file
            
            # Run MCNP
            mcnp_exec_str = 'mcnp5 n=%s' % (filename)                          # MCNP5 execution command
            print mcnp_exec_str                                                # print MCNP5 execution command
            call(mcnp_exec_str,shell=True)                                     # run MCNP5
            
            # Parse MCNP output file
            parse_output_file(filename,N,F,geom_type,'Cube_results.txt')       # parse output file for CTM and NPS

# Display output in command window
f = open('Cube_results.txt','r')
print f.read()
f.close()
