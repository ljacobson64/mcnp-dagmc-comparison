"""

Automated MCNP5 Input File Generator

Author: Lucas Jacobson

Advisers: P.P.H. Wilson, Andrew Davis

The file "params.txt" should contain the following entries:
Line  1: Directory to place input files. Must include final slash.
Line  2: Directory to place results files. Must include final slash.
Line  3: List of versions. "nat" represents native MCNP5 and "dag"
         represents DAG-MCNP5.
Line  4: List of geometry configurations. "2s" represents two-
         dimensional steps with separate cell definitions; "2j"
         represents two-dimensional steps with joined cell definitions.
         "3s" and "3j" may become options in the future.
Line  5: List of values of N (integer from 1 to 40000). N is the number
         of steps in the geometry.
Line  6: List of values of F (0.001 to 100). F is the portion of the
         geometry containing the steps. For example, an F of 0.001 would
         mean that all the steps would be condensed into a very small
         space in a corner, whereas an F of 100 would mean that the
         steps would comprise an entire square diagonal.
Line  7: List of values of mfp (0 to 100). mfp is the number of mean-
         free paths of pure deuterium in one meter. The required density
         is obtained simply by multiplying mfp by 0.0078958.
Line  8: Computer time in minutes.
Line  9: If running in DAG-MCNP mode and mfp_in > 0, CUBIT and
         dagmc_preproc are not run. Rather, an existing .h5m file with
         identical geometry but different density is used, and the only
         things changed are the density entries on the LCAD card.
Line 10: Faceting tolerance type (F, L, or A) followed by a list of
         values of the faceting tolerance. (DAG-MCNP5 only)

"""


from subprocess import call
import itertools
import math
import time


# Determine filename
def determine_filename(version,geom,N,F,mfp,tol):
    
    N_str   = ('%5u'   % (N  )).replace(' ','_')
    F_str   = ('%2.0f' % (F  )).replace(' ','_')
    mfp_str = ('%6.2f' % (mfp)).replace(' ','_')
    tol_str = ('%3.0f' % (tol)).replace(' ','_')
    
    fname = 'zCube_%s_%s_%s_%s_%s_%s' % \
            (version,geom,N_str,F_str,mfp_str,tol_str)
    return fname
    
    
# Write title cards
def write_title(fname,geom,N,F,rho):
    
    # Write title cards to the specified file
    writer = open(direc+fname+'.i','w')
    
    # Print the number of steps (N), the portion of the cube taken up by
    # the steps (F), the mass density of deuterium in the cube (rho) and
    # the geometry type (geom)
    print >> writer, 'Small geometric features test'
    print >> writer, 'c    N = %u'   % (N)
    print >> writer, 'c    F = %.3f' % (F)
    print >> writer, 'c  rho = %.8f' % (rho)
    print >> writer, 'c geom = %s'   % (geom)
    print >> writer, 'c'
    
    writer.close()
    
    
# Write cell cards; 2 dimensions; separate cells
def write_cells_2s(fname,geom,N,F,rho):
    
    # Determine material specification syntax
    # If rho=0, set material to void
    if rho == 0:
        mat_str = '0            '
    # Otherwise, set material number to 1 with specified rho
    else:
        mat_str = '1 -%10.8f' % (rho)
    
    # Write cell cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c CELL CARDS'
    
    # Rest of universe
    print >> writer, '   11 0                    -09000'
    print >> writer, '   12 0              09999'
    print >> writer, '   13 0              09000 -09999        -01000'
    print >> writer, '   14 0              09000 -09999  04999'
    print >> writer, '   15 0              09000 -09999  01000 '+ \
                                         '-04999        -05000'
    print >> writer, '   16 0              09000 -09999  01000 '+ \
                                         '-04999  08999'
    
    # Top-right portion of cube (source is always in this cell)
    print >> writer, '    1 %s  09000 -09999  %u -04999  %u -08999' \
                     % (mat_str,10000+N,50000+N)
    
    # Top-left portion of cube
    print >> writer, '    2 %s  09000 -09999  01000 -%u  %u -08999' \
                     % (mat_str,10000+N,50000+N)
    
    # Bottom-right portion of cube
    print >> writer, '    3 %s  09000 -09999  %u -04999  05000 -%u' \
                     % (mat_str,10000+N,50000+N)
    
    for j in range(1,N):
        
        # Left step cells
        print >> writer, '%u %s  09000 -09999  %u -%u  %u -%u' \
                         % (10000+j,mat_str,10000+j,10000+N,50000+N-j,
                            50000+N-j+1)
        
        # Right step cells
        print >> writer, '%u %s  09000 -09999  01000 -%u  %u -%u' \
                         % (50000+j,mat_str,10000+j,50000+N-j,
                            50000+N-j+1)
        
    # Bottom step cell
    print >> writer, '%u %s  09000 -09999  01000 -%u  05000 -50001' \
                     % (50000+N,mat_str,10000+N)
    print >> writer, ''
    
    writer.close()
    
    
# Write cell cards; 2 dimensions; joined cells
def write_cells_2j(fname,geom,N,F,rho):
    
    # Determine material specification syntax
    # If rho=0, set material to void
    if rho == 0:
        mat_str = '0            '
    # Otherwise, set material number to 1 with specified rho
    else:
        mat_str = '1 -%10.8f' % (rho)
    
    # Write cell cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c CELL CARDS'
    
    # Rest of universe
    print >> writer, '   11 0              09999:-09000: 04999:'+ \
                                         '-01000: 08999:-05000'
    
    # Bottom-left cell
    print >> writer, '    1 %s  09000 -09999 -04999 -08999 (' \
                     % (mat_str)
    print >> writer, '                     01000  %u:' % (50000+N)
    for j in range(1,N):
        print >> writer, '                     %u  %u:' \
                         % (10000+j,50000+N-j)
    print >> writer, '                     %u  05000)' % (10000+N)
    
    # Top-right cell (source is always in this cell)
    print >> writer, '    2 %s  09000 -09999  01000  05000 (' \
                     % (mat_str)
    for j in range(1,N):
        print >> writer, '                    -%u -%u:' \
                         % (10000+j,50000+N-j+1)
    print >> writer, '                    -%u -%u)' \
                     % (10000+N,50000+1)
    print >> writer, ''
    
    writer.close()
    
    
# Write surface cards; 2 dimensions
def write_surfaces_2(fname,geom,N,F,rho):
    
    # Write surface cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c SURFACE CARDS'
    
    # Outer cube boundaries
    print >> writer, ' 01000 PX 0'   # Left boundary
    print >> writer, ' 05000 PY 0'   # Front boundary
    print >> writer, ' 09000 PZ 0'   # Bottom boundary
    print >> writer, '*04999 PX 100' # Right boundary (reflecting)
    print >> writer, '*08999 PY 100' # Back boundary (reflecting)
    print >> writer, ' 09999 PZ 200' # Top boundary
    
    for j in range(1,N+1):
        
        # Planes normal to the x-axis for the N steps
        print >> writer, ' %u PX %11.8f' % (10000+j,F*j/N)
        
        # Planes normal to the y-axis for the N steps
        print >> writer, ' %u PY %11.8f' % (50000+j,F*j/N)
        
    print >> writer, ''
    
    writer.close()
    
    
# Write data cards; 2 dimensions
def write_data_2(fname,version,geom,N,F,rho,ctme,mfp_in):
    
    # Write data cards to the specified file
    if version == 'nat':
        writer = open(direc+fname+'.i','a')
    elif version == 'dag':
        writer = open(direc+fname+'.i','w')
    
    print >> writer, 'c DATA CARDS'
    
    # Neutron transport problem
    print >> writer, 'MODE N'
    
    if version == 'nat':
        
        # Neutron importance of 1 in the cube, 0 in the rest of universe
        if geom == '2s':
            print >> writer, 'IMP:N 0 5R 1 %uR' % (2*N+1)
        elif geom == '2j':
            print >> writer, 'IMP:N 0 1 1'
        
        # Source located slightly offset from the center of the cube
        # (when extrapolating based on the reflecting boundary
        # conditions)
        print >> writer, 'SDEF ERG=0.0000000253 X=99.999 Y=99.999 '+ \
                         'Z=99.999'
        
    elif version == 'dag':
        
        # Source located slightly offset from the center of the cube;
        # in cell 1
        print >> writer, 'SDEF ERG=0.0000000253 X=99.999 Y=99.999 '+ \
                         'Z=99.999 CELL=1'
        
        # Do not check source cell
        print >> writer, 'DAGMC CHECK_SRC_CELL=OFF'
        
    # Material 1 is pure deuterium
    print >> writer, 'M1 1002.70C 1'
    
    # Computer time (currently unused since runs are currently done via
    # continue run)
    print >> writer, 'CTME %.8f' % (ctme)
    
    print >> writer, ''
    
    writer.close()
    
    
# Write continue run script
def write_continue(ctme):
    
    writer = open(direc+'cont.i','w')
    
    print >> writer, 'CONTINUE'
    
    # Run for 16 times the specified computer time since 16 cores are
    # used on the ACI cluster
    print >> writer, 'CTME %.8f' % (ctme*16)
    
    # Decrease overhead and I/O
    print >> writer, 'PRDMP 1E12 1E12 0 100 1E12'
    
    writer.close()
    
    
# Write CUBIT commands; 2 dimensions
def write_cubit_2(fname,geom,N,F,rho):
    
    # Write CUBIT commands to the specified file
    writer = open(direc+fname+'.jou','w')
    
    # Volume 1: top-right (source is always in this brick)
    print >> writer, 'brick x %11.8f y %11.8f z 200' \
                     % (100-F,100-F)
    print >> writer, 'move volume %u location x %11.8f y %11.8f z 100' \
                     % (1,(100+F)/2,(100+F)/2)
    
    # Volume 2: top-left
    print >> writer, 'brick x %11.8f y %11.8f z 200' \
                     % (F,100-F)
    print >> writer, 'move volume %u location x %11.8f y %11.8f z 100' \
                     % (2,F/2,(100+F)/2)
    
    # Volume 3: bottom-right
    print >> writer, 'brick x %11.8f y %11.8f z 200' \
                     % (100-F,F)
    print >> writer, 'move volume %u location x %11.8f y %11.8f z 100' \
                     % (3,(100+F)/2,F/2)
    
    # Make top and right boundaries reflecting
    print >> writer, 'group "spec.reflect" add surf 5 6'
    
    # Volumes 4 to N+3: bricks left of the steps
    for j in range(1,N+1):
        print >> writer, 'brick x %11.8f y %11.8f z 200' \
                         % (j*F/N,F/N)
        print >> writer, ('move volume %u location x %11.8f y %11.8f'+ \
                          ' z 100') % (3+j,F*j/(2*N),(N-j+0.5)*F/N)
        
    # Volumes N+4 to 2N+2: bricks right of the steps
    for j in range(1,N):
        print >> writer, 'brick x %11.8f y %11.8f z 200' \
                         % ((N-j)*F/N,F/N)
        print >> writer, ('move volume %u location x %11.8f y %11.8f'+ \
                          ' z 100') \
                          % (3+N+j,(N+j)*F/(2*N),(N-j+0.5)*F/N)
        
    # Unite the volumes so that there is one volume to the left of the
    # steps and one volume to the right
    if geom == '2j':
        if N != 1:
            print >> writer, 'unite body       %s' \
                             % (' '.join(map(str,range(4,N+4))))
            print >> writer, 'unite body 1 2 3 %s' \
                             % (' '.join(map(str,range(N+4,2*N+3))))
        else:
            print >> writer, 'unite body 1 2 3'
            
    # Volume 2N+3: inner graveyard brick
    print >> writer, 'brick x 110 y 110 z 210'
    print >> writer, 'move volume %u location 50 50 100' % (2*N+3)
    
    # Volume 2N+4: outer graveyard brick
    print >> writer, 'brick x 120 y 120 z 220'
    print >> writer, 'move volume %u location 50 50 100' % (2*N+4)
    
    # Volume 2N+5: subtract inner from outer
    print >> writer, 'subtract %u from %u' % (2*N+3,2*N+4)
    
    # Add the graveyard volume to the graveyard group
    print >> writer, 'group "graveyard" add volume %u' % (2*N+5)
    
    # If the density is nonzero, specify the required materials.
    if rho != 0:
        
        # Add all volumes to the group representing material 1 and the
        # specified rho
        print >> writer, 'group "mat_1_rho_-%.8f" add volume all' \
                         % (rho)
        
        # Remove the graveyard volume from this group
        print >> writer, 'group "mat_1_rho_-%.8f" remove volume %u' \
                         % (rho,2*N+5)
        
    # Imprint and merge
    print >> writer, 'merge all'
    print >> writer, 'imprint body all'
    print >> writer, 'set geometry version 1902'
    print >> writer, 'set attribute on'
    
    # Export to .sat file
    print >> writer, 'export acis "%s%s.sat"' % (direc,fname)
    
    writer.close()
    
    
# Modify LCAD files to change density
def modify_lcad(fname_in,fname_out,rho):
    
    # Read from the input LCAD file and write to the output file
    reader = open(direc+fname_in+'.lcad','r')
    writer = open(direc+fname_out+'.lcad','w')
    
    lines = reader.readlines()
    
    # New density is zero
    if rho == 0:
        
        # Input rho string contains material ID and density
        in_str = ' '.join(lines[0].split()[1:3])
        
        # Output rho string is zero
        out_str = str(0)
        
    # New density is nonzero
    else:
        
        # Input rho string contains material density
        in_str = lines[0].split()[2]
        
        # Output rho string is the new material density
        out_str = '-%10.8f' % (rho)
        
    # Replace the original rho string with new rho string; leave
    # everything else as is
    line_ind = 0
    for line in lines:
        line_ind += 1
        if not line.strip():
            writer.write(line)
            break
        elif in_str in line:
            writer.write(line.replace(in_str,out_str))
        else:
            writer.write(line)
    for line in lines[line_ind:]:
        writer.write(line)
    
    reader.close()
    writer.close()
    
    
# Write script to create runtpe files for native runs and do CUBIT and
# DAGMC-preprocessing and create runtpe files for DAG runs
def write_setup(fname,fname_in,version,toltype,tol):
    
    writer = open(direc+fname+'.setup.sh','w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    
    # Native MCNP
    if version == 'nat':
        
        # Create runtpe file; do not simulate any particles
        print >> writer, 'mcnp5.mpi ix n=%s.i o=%s.setup.io' \
                         % (fname,fname)
    
    # DAG-MCNP
    elif version == 'dag':
        
        # Not using pre-existing FCAD file
        if fname == fname_in:
            
            # Play CUBIT journal file
            print >> writer, ('cubit -batch -nographics -nojournal '+ \
                              '-information=off %s.jou') % (fname)
            
            # Run DAGMC-preprocessing with the specified faceting
            # tolerance
            print >> writer, 'dagmc_preproc -%s %.2e %s.sat -o %s.h5m' \
                             % (toltype,tol,fname,fname)
            
            # Create runtpe file; do not simulate any particles
            print >> writer, ('mcnp5.mpi ix n=%s.i g=%s.h5m '+ \
                              'o=%s.setup.io l=%s.lcad '+ \
                              'f=%s.setup.fcad') \
                              % (fname,fname,fname,fname,fname)
            
            # Create .stl file for viewing in VisIt
           #print >> writer, 'mbconvert %s.h5m %s.stl' % (fname,fname)
            
        # Using pre-existing FCAD file
        else:
            
            # Create runtpe file; do not simulate any particles
            print >> writer, ('mcnp5.mpi ix n=%s.i g=%s.h5m '+ \
                              'o=%s.setup.io l=%s.lcad '+ \
                              'f=%s.setup.fcad') \
                              % (fname,fname_in,fname,fname,fname)
            
            # Remove the duplicate FCAD file
            print >> writer, 'rm %s.setup.fcad' % (fname)
            
    # Save a copy of the original runtpe file
    print >> writer, 'cp %s.ir %s.setup.ir' % (fname,fname)
    
    writer.close()
    
    
# Write script to run all setup scripts
def write_setup_master(fname,fname_in,fid,max_concurrent):
    
    if fname == fname_in:
        writer = open(direc+'setup_files.sh','a')
    else:
        writer = open(direc+'setup_files_2.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s' % (fname)
    
    # Create runtpe file, wait for completion
    if (fid % max_concurrent) == 0:
        print >> writer, 'sh %s.setup.sh > %s.setup.out' \
                         % (fname,fname)
        
    # Create runtpe file, do not wait for completion
    else:
        print >> writer, 'sh %s.setup.sh > %s.setup.out &' \
                         % (fname,fname)
        
    writer.close()
    
# Write local run script
def write_local_run(fname,fname_in,version):
    
    writer = open(direc+'local_runs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
        
    print >> writer, 'echo %s' % (fname)
    
    # Native MCNP
    if version == 'nat':
        print >> writer, ('mpirun -np 2 mcnp5.mpi c i=cont.i'+ \
                          ' o=%s.io r=%s.ir > %s.out') \
                          % (fname,fname,fname)
        
    # DAG-MCNP
    elif version == 'dag':
        
        # Do not use pre-existing FCAD file
        if fname == fname_in:
            print >> writer, ('mpirun -np 2 mcnp5.mpi c i=cont.i '+ \
                              'g=%s.setup.fcad o=%s.io r=%s.ir '+ \
                              'l=%s.lcad f=%s.fcad > %s.out') \
                              % (fname,fname,fname,fname,fname,fname)
        
        # Use pre-existing FCAD file
        else:
            print >> writer, ('mpirun -np 2 mcnp5.mpi c i=cont.i '+ \
                              'g=%s.setup.fcad o=%s.io r=%s.ir '+ \
                              'l=%s.lcad f=%s.fcad > %s.out') \
                              % (fname_in,fname,fname,fname,fname,fname)
            
    writer.close()
    
    
# Write ACI job script
def write_job(fname,fname_in,version):
    
    # Extra time to add to the ACI runs (in minutes)
    time_extra = 10
    
    writer = open(direc+fname+'.aci.sh','w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    
    # If the run time is less than one hour, use the pre-emptable
    # partition. Otherwise, use the university partition.
    if time_extra+ctme <= 60:
        print >> writer, '#SBATCH --partition=pre'
    else:
        print >> writer, '#SBATCH --partition=univ'
    
    # Run time in days-hh:mm:ss
    print >> writer, '#SBATCH --time=0-00:%u:00' % (time_extra+ctme)
    
    # Use 16 tasks
    print >> writer, '#SBATCH --ntasks=16'
    
    # Request 1 GB of RAM (max 8GB)
    print >> writer, '#SBATCH --mem-per-cpu=1024'
    
    # Specify location of command error and output files
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err' % (fname)
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out' \
                     % (fname)
    
    print >> writer, ''
    
    # Native MCNP
    if version == 'nat':
        print >> writer, ('srun /home/adavis23/dagmc/mcnp5v16_dag/'+ \
                          'bin/mcnp5.mpi c i=cont.i o=%s.io r=%s.ir') \
                          % (fname,fname)
    
    # DAG-MCNP
    elif version == 'dag':
        
        # Do not use pre-existing FCAD file
        if fname == fname_in:
            print >> writer, ('srun /home/adavis23/dagmc/'+ \
                              'mcnp5v16_dag/bin/mcnp5.mpi c '+ \
                              'i=cont.i g=%s.setup.fcad o=%s.io '+ \
                              'r=%s.ir l=%s.lcad f=%s.fcad') \
                              % (fname,fname,fname,fname,fname)
        
        # Use pre-existing FCAD file
        else:
            print >> writer, ('srun /home/adavis23/dagmc/'+ \
                              'mcnp5v16_dag/bin/mcnp5.mpi c '+ \
                              'i=cont.i g=%s.setup.fcad o=%s.io '+ \
                              'r=%s.ir l=%s.lcad f=%s.fcad') \
                              % (fname_in,fname,fname,fname,fname)
    
    print >> writer, ''
    
    writer.close()
    
    
# Write master job script
def write_job_master(fname):
    
    writer = open(direc+'submit_jobs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s.sh' % (fname)
    
    # Execute the job script
    print >> writer, 'sbatch %s.aci.sh' % (fname)
    
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
toltype  =                    params_str[9].split()[1 ]                         # type of faceting tolerance to use for DAG
tols     = [float(i) for i in params_str[9].split()[2:]]                        # faceting tolerance to use for DAG

reader.close()

min_N        =      1
max_N_nat_2s =  40000
max_N_nat_2j =   2000
max_N_dag    =  40000
min_F        =      0.1
max_F        =     99.9

fid = 0

# Write input files and job scripts
for params in list(itertools.product(versions,geoms,N_vals,F_vals,mfps,tols)):
    
    fid = fid+1
    
    version = params[0]                                                         # loop for all MCNP versions
    geom    = params[1]                                                         # loop for all geometry configurations
    N       = params[2]                                                         # loop for all values of N
    F       = params[3]                                                         # loop for all values of F
    mfp     = params[4]                                                         # loop for all values of mfp
    tol     = params[5]                                                         # loop for all values of tol
    
    # Calculate mass density of deuterium from input number of mean free paths in 1 meter
    # Setting the mean free path to 1 meter results in a density of 0.0078958 g/cm^3
    rho = mfp*0.0078958
    
    # Determine whether parameters are valid
    valid_version = version == 'nat' or version == 'dag'                        # version must be nat or dag
    valid_geom = geom == '2s' or geom == '2j'                                   # geom must be 2s or 2j
    
    if version == 'nat' and geom == '2s':                                       # N must be between the minimum value and maximum value
        valid_N = N >= min_N and N <= max_N_nat_2s
    elif version == 'nat' and geom == '2j':
        valid_N = N >= min_N and N <= max_N_nat_2j
    elif version == 'dag' and (geom == '2s' or geom == '2j'):
        valid_N = N >= min_N and N <= max_N_dag
    else:
        valid_N = False
    
    valid_F = F >= min_F                                                        # F must be between the minimum value and maximum value
    valid_mfp = mfp >= 0                                                        # mfp must be greater than or equal to zero
    
    if not valid_version or not valid_geom or not valid_N or not valid_F or not valid_mfp:
        print 'Invalid parameters'
        continue

    if version == 'nat':                                                        # maximum number of setup jobs running at once
        max_concurrent = len(mfps)*2
    elif (version == 'dag') and (mfp_in == 0):
        max_concurrent = len(F_vals)
    elif (version == 'dag') and (mfp_in > 0):
        max_concurrent = len(mfps)*2
        
    fname = determine_filename(version,geom,N,F,mfp,tol)                        # base filename (no extension)
    print fname
    fname_in = determine_filename(version,geom,N,F,mfp,tol)
    
    # Write input files
    if version == 'nat' and (geom == '2s' or geom == '2j'):                     # Native MCNP, 2 dimensions
        
        write_title(fname,geom,N,F,rho)                                         # Write title cards
        if geom == '2s':
            write_cells_2s(fname,geom,N,F,rho)                                  # Append cell cards; separate cells
        elif geom == '2j':
            write_cells_2j(fname,geom,N,F,rho)                                  # Append cell cards; joined cells
        write_surfaces_2(fname,geom,N,F,rho)                                    # Append surface cards
        
    elif version == 'dag' and (geom == '2s' or geom == '2j'):                   # DAG-MCNP, 2 dimensions
        
        if mfp_in > 0:                                                          # Use an existing LCAD file with a different rho
            fname_in = determine_filename(version,geom,N,F,mfp_in,tol)
            modify_lcad(fname_in,fname,rho)
        else:
            write_cubit_2(fname,geom,N,F,rho)                                   # Write CUBIT instructions
        
    write_data_2(fname,version,geom,N,F,rho,ctme,mfp_in)                        # Write MCNP data cards
    
    write_setup(fname,fname_in,version,toltype,tol)                             # Write script to run CUBIT and DAGMC pre-processing and/or create runtpe file
    write_setup_master(fname,fname_in,fid,max_concurrent)                       # Append instructions to setup file
    write_local_run(fname,fname_in,version)                                     # Append instructions to local run file
    write_job(fname,fname_in,version)                                           # Write ACI job file
    write_job_master(fname)                                                     # Append instructions to master job file
    
# Write continue run script
write_continue(ctme)

# Make all scripts executable
call('chmod 740 '+direc+'*.sh',shell=True)
