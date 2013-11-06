# Write input file; 2 dimensions with separate cells
def write_mcnp_input_2s(filename,N,F,ctme):
    
    filename_i = '%s.i' % (filename)
    writer = open(filename_i,'w')
    
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
	
    filename_i = '%s.i' % (filename)
    writer = open(filename_i,'w')
    
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

# Write job script
def write_job_script(filename):
    
    filename_j = '%s.sh' % (filename)
    writer = open(filename_j,'w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, 'get_until_got(){' 
    print >> writer, 'wget  -c -t 5 --waitretry=20 --read-timeout=10 $1 '
    print >> writer, 'while [[ $? != 0 ]]'
    print >> writer, 'do'
    print >> writer, 'wget $1'
    print >> writer, 'done'
    print >> writer, '}'
    print >> writer, ''
    print >> writer, '# Determine absolute path'
    print >> writer, 'cwd=$PWD'
    print >> writer, ''
    print >> writer, '# Get and set the gcc compiler suite and set ld and paths'
    print >> writer, 'get_until_got http://proxy.chtc.wisc.edu/SQUID/lucasj/compiler_tools.tar.gz' 
    print >> writer, ''
    print >> writer, '# Unpack compiler tools'
    print >> writer, 'tar -zxf compiler_tools.tar.gz '
    print >> writer, ''
    print >> writer, '# Set library paths'
    print >> writer, 'export LD_LIBRARY_PATH=$cwd/compiler/gcc-4.8.1/lib:$cwd/compiler/gcc-4.8.1/lib64:$cwd/compiler/gmp-5.1.2/lib:$cwd/compiler/mpc-1.0.1/lib:$cwd/compiler/mpfr-3.1.2/lib  '
    print >> writer, ''
    print >> writer, '# Get and set the MOAB and HDF5 libs' 
    print >> writer, 'get_until_got http://proxy.chtc.wisc.edu/SQUID/ljjacobson/moab_tools.tar.gz '
    print >> writer, ''
    print >> writer, '# set the MOAB path'
    print >> writer, 'tar -zxf moab_tools.tar.gz '
    print >> writer, 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$cwd/moab_tools/hdf5-1.8.4/lib'
    print >> writer, 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$cwd/moab_tools/moab-4.6.0/lib' 
    print >> writer, ''
    print >> writer, '# Get and set the required MCNP5 paths '
    print >> writer, 'get_until_got http://proxy.chtc.wisc.edu/SQUID/ljjacobson/mcnp5.tar.gz '
    print >> writer, 'mkdir mcnp5 '
    print >> writer, 'cp mcnp5.tar.gz mcnp5/.' 
    print >> writer, 'cd mcnp5 '
    print >> writer, 'tar -zxf mcnp5.tar.gz' 
    print >> writer, 'cd .. '
    print >> writer, 'export PATH=$cwd/mcnp5:$PATH '
    print >> writer, ''
    print >> writer, 'get_until_got %s.i' % (filename)
    print >> writer, 'mcnp5 i=%s.i o=%s.io' % (filename,filename)
    print >> writer, 'ls | grep -v %s.io | xargs rm -rf' % (filename)
    
    writer.close()
    
# Write command file
def write_command_file(filename):
    
    filename_c = '%s.cmd' % (filename)
    writer = open(filename_c,'w')
    
    print >> writer, 'executable = %s.sh' % (filename)
    print >> writer, ''
    print >> writer, 'copy_to_spool = false' 
    print >> writer, 'should_transfer_files = yes' 
    print >> writer, 'when_to_transfer_output = on_exit' 
    print >> writer, 'output = %s.out' % (filename)
    print >> writer, 'log = %s.log' % (filename)
    print >> writer, 'error = %s.err' % (filename)
    print >> writer, 'transfer_input_files = %s.sh' % (filename)
    print >> writer, '+AccountingGroup = EngrPhysics_Wilson' 
    print >> writer, 'Queue'
    
    writer. close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params = reader.readlines()

N_vals     = [  int(i) for i in params[0].split()[1:]]                         # list of values of N
F_vals     = [float(i) for i in params[1].split()[1:]]                         # list of values of F
geom_types =                    params[2].split()[1:]                          # list of geometry configurations
ctme       =  float(            params[3].split()[1 ])                         # computer time

reader.close()

# Write input files and job scripts
for geom_type in geom_types:                                                   # loop for all geometry configurations
    for N in N_vals:                                                           # loop for all values of N
        for F in F_vals:                                                       # loop for all values of F
            
            filename = 'Cube_%s_%u_%.0f' % (geom_type,N,F)                     # Base filename (no extension) 
            
            # Write MCNP input file
            if geom_type=='2s':                                                # 2 dimensions with separate cells
                write_mcnp_input_2s(filename,N,F,ctme)                         # write input file
            if geom_type=='2j':                                                # 2 dimensions with joined cells
                write_mcnp_input_2j(filename,N,F,ctme)                         # write input file
            
            # Write job script
            write_job_script(filename)                                         # write job script
            
            # Write command file
            write_command_file(filename)                                       # write command file
            
