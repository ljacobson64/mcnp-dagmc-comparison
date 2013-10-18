from subprocess import call
import os

def write_mcnp_input_2s(filename,N,F,ctme):
	
    writer = open(filename,'w')
    
    # Write title cards
    print >> writer, 'Small geometric features test'
    print >> writer, 'c Cube geometry with separate cells in 2 dimensions'
    print >> writer, 'c N = %u; F = %4.1f' % (N,F)
    
    # Write cell cards
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0 -09000'                             # Outside cube
    print >> writer, '   12 0  09999'                             # Outside cube
    print >> writer, '   13 0  09000 -09999 -01000'               # Outside cube
    print >> writer, '   14 0  09000 -09999  04999'               # Outside cube
    print >> writer, '   15 0  09000 -09999  01000 -04999 -05000' # Outside cube
    print >> writer, '   16 0  09000 -09999  01000 -04999  08999' # Outside cube
    print >> writer, '    1 0  09000 -09999  %u -04999  %u -08999' % (10000+N,50000+N)
    print >> writer, '    2 0  09000 -09999  01000 -%u  %u -08999' % (10000+N,50000+N)
    print >> writer, '    3 0  09000 -09999  %u -04999  05000 -%u' % (10000+N,50000+N)
    for j in range(1,N):
        print >> writer, '%u 0  09000 -09999  %u -%u  %u -%u' % (10000+j,10000+j,10000+N,50000+N-j,50000+N-j+1)
        print >> writer, '%u 0  09000 -09999  01000 -%u  %u -%u' % (50000+j,10000+j,50000+N-j,50000+N-j+1)
    print >> writer, '%u 0  09000 -09999  01000 -%u  05000 -%u' % (50000+N,10000+N,50000+N-N+1)
    print >> writer, ''
    
    # Write surface cards
    print >> writer, 'c SURFACE CARDS'
    print >> writer, '01000 PX 0'
    print >> writer, '05000 PY 0'
    print >> writer, '09000 PZ 0'   # bottom side
    print >> writer, '04999 PX 100'
    print >> writer, '08999 PY 100'
    print >> writer, '09999 PZ 100' # top side
    for j in range(1,N+1):
        print >> writer, '%u PX %11.8f' % (10000+j,F*j/N)
        print >> writer, '%u PY %11.8f' % (50000+j,F*j/N)
    print >> writer, ''
    
    # Write data cards
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'IMP:N 0 5R 1 %uR' % (2*N+1)
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % ctme
    
    writer.close()

def write_mcnp_input_2j(filename,N,F,ctme):
	
    writer = open(filename,'w')
    
    # Write title cards
    print >> writer, 'Small geometric features test'
    print >> writer, 'c Cube geometry with joined cells in 2 dimensions'
    print >> writer, 'c N = %u; F = %4.1f' % (N,F)
    
    # Write cell cards
    print >> writer, 'c CELL CARDS'
    print >> writer, '    1 0 -1 (-%u:' % (10000+1)
    
    for j in range(2,N):
        print >> writer, '            -%u:' % (10000+j)  
    print >> writer, '            -%u)' % (10000+N)
    print >> writer, '    2 0 -1   %u' % (10000+1)
    for j in range(2,N+1):
    	print >> writer, '             %u' % (10000+j)
    print >> writer, '    3 0  1'
    print >> writer, ''
    
    # Write surface cards
    print >> writer, 'c SURFACE CARDS'
    print >> writer, '    1 RPP 0 100 0 100 0 100'
    for j in range(1,N+1):
        print >> writer, '%u RPP 0 %11.8f 0 %11.8f 0 100' % (10000+j, F*(1-(j-1)/float(N)), F*j/N)
    print >> writer, ''
    
    # Write data cards
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'IMP:N 1 1 0'
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % ctme
    
    writer.close()

def parse_output_file(filename,N,F,filename_results):
    
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
    
    # Find CTM
    ctm_loc = str.find(lines,'ctm')
    ctm_str = lines[ctm_loc+5:ctm_loc+18]
    ctm = float(ctm_str)
    
    writer = open(filename_results,'a')
    
    # Write NPS to text file
    if writer.tell() == 0:
        print >> writer, '|-------|-----|-------|----------|'
        print >> writer, '|   N   |  F  |  ctm  |   nps    |'
        print >> writer, '|-------|-----|-------|----------|'
    print >> writer, '| %5u | %3.0f | %5.2f | %8u |' % (N,F,ctm,nps)
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params = reader.readlines()

N_vals     = [  int(i) for i in params[0].split()[1:]]
F_vals     = [float(i) for i in params[1].split()[1:]]
geom_types =                    params[2].split()[1:]
ctme       =  float(            params[3].split()[1 ])

reader.close()

# Delete all old files
[os.remove(f) for f in os.listdir('.') if f.startswith('Cube_')]

for geom_type in geom_types:
    for N in N_vals:
        for F in F_vals:
            
            # Write MCNP input file
            filename = 'Cube_%s_%u_%.0f.i' % (geom_type,N,F)
            if geom_type=='2s': write_mcnp_input_2s(filename,N,F,ctme)
            if geom_type=='2j': write_mcnp_input_2j(filename,N,F,ctme)
            
            # Run MCNP
            mcnp_exec_str = 'mcnp5 n=%s' % (filename)
            print mcnp_exec_str
            call(mcnp_exec_str,shell=True)
            os.remove(filename+'r')
            
            # Parse MCNP output file
            filename_results = 'Cube_%s_results.txt' % geom_type
            parse_output_file(filename,N,F,filename_results)

for geom_type in geom_types:
    filename_results = 'Cube_%s_results.txt' % geom_type
    f = open(filename_results,'r')
    print f.read()
    f.close()




