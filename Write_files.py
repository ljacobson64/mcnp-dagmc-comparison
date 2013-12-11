from subprocess import call

# Write title cards
def write_title(fname,N,F,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'w')
    
    print >> writer, 'Small geometric features test'
    print >> writer, 'c Cube geometry with separate cells in 2 dimensions'
    print >> writer, 'c N = %u; F = %4.1f' % (N,F)
    
    writer.close()

# Write cell cards; 2 dimensions, separate cells
def write_cells_2s(fname,N,F,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
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
    
    writer.close()

# Write cell cards; 2 dimensions, joined cells
def write_cells_2j(fname,N,F,ctme):
	
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
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
    
    writer.close()

# Write surface cards; 2 dimensions, separate or joined cells
def write_surfaces_2(fname,N,F,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
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
    
    writer.close()

# Write data cards; native MCNP, 2 dimensions, separate cells
def write_data_2s(fname,N,F,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'IMP:N 0 5R 1 %uR' % (2*N+1)
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % (ctme)
    
    writer.close()

# Write data cards; native MCNP, 2 dimensions, joined cells
def write_data_2j(fname,N,F,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'IMP:N 0 1 1'
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % ctme
    
    writer.close()

# Write data cards; DAG-MCNP
def write_data_2d(fname,N,F,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'SDEF X=99.9999 Y=99.9999 Z=99.9999'
    print >> writer, 'CTME %f' % ctme
    
    writer.close()

# Write job script; native MCNP
def write_jobs_n(fname,N,F,ctme,geom):
    
    fname_sh = '%s.sh' % (fname)
    writer = open(direc+fname_sh,'w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, '#SBATCH --partition=univ'                                 # default "univ" if not specified
    
    if   geom == '2s' and N > 4000:
        print >> writer, '#SBATCH --time=0-01:00:00'                            # run time in days-hh:mm:ss
    elif geom == '2j' and N >  200:
        print >> writer, '#SBATCH --time=0-01:00:00'                            # run time in days-hh:mm:ss
    else:
        print >> writer, '#SBATCH --time=0-00:06:00'                            # run time in days-hh:mm:ss
    
    print >> writer, '#SBATCH --ntasks=1'                                       # number of CPUs
    print >> writer, '#SBATCH --mem-per-cpu=2000'                               # RAM in MB (default 4GB, max 8GB)
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err' % (fname)
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out' % (fname)
    print >> writer, ''
    print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i' % (fname)
    
    writer.close()
    # Write job script; MCNP

# Write job script; DAG-MCNP
def write_jobs_d(fname,N,F,ctme,geom):
    
    fname_sh = '%s.sh' % (fname)
    writer = open(direc+fname_sh,'w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, '#SBATCH --partition=univ'                                 # default "univ" if not specified
    
    if   geom == '2s' and N > 4000:
        print >> writer, '#SBATCH --time=0-01:00:00'                            # run time in days-hh:mm:ss
    elif geom == '2j' and N  > 200:
        print >> writer, '#SBATCH --time=0-01:00:00'                            # run time in days-hh:mm:ss
    else:
        print >> writer, '#SBATCH --time=0-00:06:00'                            # run time in days-hh:mm:ss
    
    print >> writer, '#SBATCH --ntasks=1'                                       # number of CPUs
    print >> writer, '#SBATCH --mem-per-cpu=2000'                               # RAM in MB (default 4GB, max 8GB)
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err' % (fname)
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out' % (fname)
    print >> writer, ''
    print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i g=%s.h5m' % (fname,fname)
    
    writer.close()
    
# Write master job script
def write_master_jobs(fname):
    
    writer = open(direc+'submit_jobs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s.sh' % (fname)
    print >> writer, 'sbatch %s.sh' % (fname)
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params = reader.readlines()

global direc
direc      =                    params[0].split()[0 ]                           # directory to place files
N_vals     = [  int(i) for i in params[1].split()[1:]]                          # list of values of N
F_vals     = [float(i) for i in params[2].split()[1:]]                          # list of values of F
geoms      =                    params[3].split()[1:]                           # list of geometry configurations
ctme       =  float(            params[4].split()[1 ])                          # computer time
versions   =                    params[5].split()[1:]                           # list of MCNP versions

reader.close()

max_N_2s = 40000
max_N_2j =  1000

use_local_hd = True

# Write input files and job scripts
for version in versions:                                                        # loop for all MCNP versions
    for geom in geoms:                                                          # loop for all geometry configurations
        for N in N_vals:                                                        # loop for all values of N
            for F in F_vals:                                                    # loop for all values of F
                
                # Determine whether parameters are valid
                if   version == 'nat' and geom == '2s':
                    valid_N = N <= max_N_2s
                elif version == 'nat' and geom == '2j':
                    valid_N = N <= max_N_2j
                elif version == 'dag' and (geom == '2s' or geom == '2j'):
                    valid_N = True
                else:
                    valid_N = False
                valid_F = F >= 0 and F <= 100
                
                fname = 'zCube_%s_%s_%u_%.0f' % (version,geom,N,F)              # base filename (no extension)
                
                # Write input files
                if valid_N and valid_F:
                    
                    if version == 'nat' and (geom == '2s' or geom == '2j'):     # Native MCNP, 2 dimensions
                        
                        write_title(fname,N,F,ctme)                             # Write title cards
                        
                        if   geom == '2s':
                            write_cells_2s(fname,N,F,ctme)                      # Append cell cards; separate cells
                        elif geom == '2j':
                            write_cells_2j(fname,N,F,ctme)                      # Append cell cards; joined cells
                        
                        write_surfaces_2(fname,N,F,ctme)                        # Append surface cards
                        
                        if   geom == '2s':
                            write_data_2s(fname,N,F,ctme)                       # Append data cards; separate cells
                        elif geom == '2j':
                            write_data_2j(fname,N,F,ctme)                       # Append data cards; joined cells
                        
                        write_jobs_n(fname,N,F,ctme,geom)                       # Write job file
                        write_master_jobs(fname)                                # Append instructions to master job file
                        
                    if version == 'dag' and (geom == '2s' or geom == '2j'):     # DAG-MCNP, 2 dimensions
                        
                        write_title(fname,N,F,ctme)                             # Write title cards
                        write_cells_2s(fname,N,F,ctme)                          # Append cell cards; separate cells
                        write_surfaces_2(fname,N,F,ctme)                        # Append surface cards; separate cells
                        write_data_2s(fname,N,F,ctme)                           # Append data cards; separate cells
                        call('mcnp2cad '+direc+fname+'.i -o '+direc+fname+'.sat',shell=True)
                                                                                # Convert MCNP to CAD
                      # call('mv out.sat '+direc+fname+'.sat',shell=True)       # Rename CAD file
                        call('rm -f '+direc+fname+'.i',shell=True)              # Delete original MCNP input file
                        write_data_2d(fname,N,F,ctme)                           # Write MCNP data cards for DAG
                      # if geom == '2j':
                            # CUBIT tinkering                                   # Join cells in CUBIT
                        call('dagmc_preproc -f 1.0e-4'+direc+fname+'.sat -o '+direc+fname+'.h5m',shell=True)
                                                                                # DAGMC pre-processing
                        write_jobs_d(fname,N,F,ctme,geom)                       # Write job file
                        write_master_jobs(fname)                                # Append instructions to master job file
                    
                    print direc+fname

# Make master job script executable
call('chmod 770 '+direc+'submit_jobs.sh',shell=True)
