"""

Automated MCNP5 Input File Generator

Author: Lucas Jacobson

Advisors: P.P.H. Wilson, Andrew Davis

The file "params.txt" should contain the following entries:
Line 1: Directory to place input files. Must include final slash.
Line 2: Directory to place results files. Must include final slash.
Line 3: List of versions. "nat" represents native MCNP5 and "dag" represents DAG-MCNP5.
Line 4: List of geometry configurations. "2s" represents two-dimensional steps with separate cell
        definitions; "2j" represents two-dimensional steps with joined cell definitions. "3s" and
        "3j" may become options in the future.
Line 5: List of values of N (integer from 1 to 40000). N is the number of steps in the geometry.
Line 6: List of values of F (0.001 to 100). F is the portion of the geometry containing the steps.
        For example, an F of 0.001 would mean that all the steps would be condensed into a very
        small space in a corner, whereas an F of 100 would mean that the steps would comprise an
        entire square diagonal.
Line 7: List of values of mfp (0 to 100). mfp is the number of mean-free paths of pure deuterium in
        one meter. The required density is obtained simply by multiplying mfp by 0.0078958.
Line 8: Computer time in minutes.
Line 9: If running in DAG-MCNP mode and mfp_in > 0, CUBIT and dagmc_preproc are not run. Rather, an
        existing .h5m file with identical geometry but different density is used, and the only
        things changed are the density entries on the LCAD card.

"""

from subprocess import call
import itertools
import math
import time

# Determine filename
def determine_filename(version,geom,N,F,mfp):
    
    N_str   = ('%5u'   % (N  )).replace(' ','_')
    F_str   = ('%2.0f' % (F  )).replace(' ','_')
    mfp_str = ('%6.2f' % (mfp)).replace(' ','_')
    
    fname = 'zCube_%s_%s_%s_%s_%s' % (version,geom,N_str,F_str,mfp_str)         # base filename (no extension)
    return fname

# Write title cards
def write_title(fname,geom,N,F,rho):
    
    # Write title cards to the specified file
    writer = open(direc+fname+'.i','w')
    
    print >> writer, 'Small geometric features test'
    print >> writer, 'c    N = %u'   % (N)                                      # number of steps
    print >> writer, 'c    F = %.3f' % (F)                                      # portion of cube taken up by the steps
    print >> writer, 'c  rho = %.8f' % (rho)                                    # mass density of deuterium in the cube
    print >> writer, 'c geom = %s'   % (geom)                                   # type of geometry (2D, 3D, separate cells, joined cells)
    print >> writer, 'c'
    
    writer.close()

# Write cell cards; 2 dimensions; separate cells
def write_cells_2s(fname,geom,N,F,rho):
    
    # Determine material specification syntax
    if rho == 0:
        mat_str = '0            '                                               # if rho=0, set material to void
    else:
        mat_str = '1 -%10.8f' % (rho)                                           # otherwise, set material number to 1 with specified rho
    
    # Write cell cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0                    -09000'                        # rest of universe
    print >> writer, '   12 0              09999'
    print >> writer, '   13 0              09000 -09999        -01000'
    print >> writer, '   14 0              09000 -09999  04999'
    print >> writer, '   15 0              09000 -09999  01000 -04999        -05000'
    print >> writer, '   16 0              09000 -09999  01000 -04999  08999'
    print >> writer, '    1 %s  09000 -09999  %u -04999  %u -08999' % (mat_str,10000+N,50000+N) # top-right portion of cube (source is always in this cell)
    print >> writer, '    2 %s  09000 -09999  01000 -%u  %u -08999' % (mat_str,10000+N,50000+N) # top-left portion of cube
    print >> writer, '    3 %s  09000 -09999  %u -04999  05000 -%u' % (mat_str,10000+N,50000+N) # bottom-right portion of cube
    for j in range(1,N):
        print >> writer, '%u %s  09000 -09999  %u -%u  %u -%u'      % (10000+j,mat_str,10000+j,10000+N,50000+N-j,50000+N-j+1) # left step cells
        print >> writer, '%u %s  09000 -09999  01000 -%u  %u -%u'   % (50000+j,mat_str,10000+j,50000+N-j,50000+N-j+1)         # right step cells
    print >> writer, '%u %s  09000 -09999  01000 -%u  05000 -50001' % (50000+N,mat_str,10000+N)                               # bottom step cell
    print >> writer, ''
    
    writer.close()

# Write cell cards; 2 dimensions; joined cells
def write_cells_2j(fname,geom,N,F,rho):
	
    # Determine material specification syntax
    if rho == 0:
        mat_str = '0            '                                               # if rho=0, set material to void
    else:
        mat_str = '1 -%10.8f' % (rho)                                           # otherwise, set material number to 1 with specified rho
    
    # Write cell cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0              09999:-09000: 04999:-01000: 08999:-05000' # rest of universe
    print >> writer, '    1 %s  09000 -09999 -04999 -08999 ('      % (mat_str)       # bottom-left cell
    print >> writer, '                     01000  %u:'             % (50000+N)
    for j in range(1,N):
        print >> writer, '                     %u  %u:'            % (10000+j,50000+N-j)
    print >> writer, '                     %u  05000)'             % (10000+N)
    print >> writer, '    2 %s  09000 -09999  01000  05000 ('      % (mat_str)       # top-right cell (source is always in this cell)
    for j in range(1,N):
        print >> writer, '                    -%u -%u:'            % (10000+j,50000+N-j+1)
    print >> writer, '                    -%u -%u)'                % (10000+N,50000+1)
    print >> writer, ''
    
    writer.close()

# Write surface cards; 2 dimensions
def write_surfaces_2(fname,geom,N,F,rho):
    
    # Write surface cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c SURFACE CARDS'
    print >> writer, ' 01000 PX 0'                                              # left boundary of cube
    print >> writer, ' 05000 PY 0'                                              # front boundary of cube
    print >> writer, ' 09000 PZ 0'                                              # bottom boundary of cube
    print >> writer, '*04999 PX 100'                                            # right boundary of cube (reflecting)
    print >> writer, '*08999 PY 100'                                            # back boundary of cube (reflecting)
    print >> writer, ' 09999 PZ 200'                                            # top boundary of cube
    for j in range(1,N+1):
        print >> writer, ' %u PX %11.8f' % (10000+j,F*j/N)                      # planes normal to the x-axis for the N steps
        print >> writer, ' %u PY %11.8f' % (50000+j,F*j/N)                      # planes normal to the y-axis for the N steps
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
    print >> writer, 'MODE N'                                                   # neutron transport problem
    if version == 'nat':
        if geom == '2s':                                                        # neutron importance of 1 in the cube, 0 in the rest of universe
            print >> writer, 'IMP:N 0 5R 1 %uR' % (2*N+1)
        elif geom == '2j':
            print >> writer, 'IMP:N 0 1 1'
        print >> writer, 'SDEF ERG=0.0000000253 X=99.999 Y=99.999 Z=99.999'
    elif version == 'dag':
        print >> writer, 'SDEF ERG=0.0000000253 X=99.999 Y=99.999 Z=99.999 CELL=1'
        print >> writer, 'DAGMC CHECK_SRC_CELL=OFF'                             # do not check source cell
    print >> writer, 'M1 1002.70C 1'                                            # material 1 is pure deuterium
    print >> writer, 'CTME %.8f' % (ctme)                                       # do not run simulation; only create runtpe file
    print >> writer, ''
    
    writer.close()

# Write continue run script
def write_continue(ctme):
    
    writer = open(direc+'cont.i','w')
    
    print >> writer, 'CONTINUE'
    print >> writer, 'CTME %.8f' % (ctme)
    print >> writer, 'PRDMP 1E12 1E12 0 100 1E12'
    
    writer.close()

# Write CUBIT commands; 2 dimensions
def write_cubit_2(fname,geom,N,F,rho):
    
    # Write CUBIT commands to the specified file
    writer = open(direc+fname+'.jou','w')
    
    print >> writer, 'brick x %11.8f y %11.8f z 200'                       % (100-F,100-F)   # volume 1: top-right (source is always in this brick)
    print >> writer, 'move volume %u location x %11.8f y %11.8f z 100'     % (1,(100+F)/2,(100+F)/2)
    print >> writer, 'brick x %11.8f y %11.8f z 200'                       % (F,100-F)       # volume 2: top-left
    print >> writer, 'move volume %u location x %11.8f y %11.8f z 100'     % (2,F/2,(100+F)/2)
    print >> writer, 'brick x %11.8f y %11.8f z 200'                       % (100-F,F)       # volume 3: bottom-right
    print >> writer, 'move volume %u location x %11.8f y %11.8f z 100'     % (3,(100+F)/2,F/2)
    print >> writer, 'group "spec.reflect" add surf 5 6'                                     # make top and right boundaries reflecting
    
    for j in range(1,N+1):
        print >> writer, 'brick x %11.8f y %11.8f z 200'                   % (j*F/N,F/N)     # volumes 4 to N+3: bricks left of the steps
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 100' % (3+j,F*j/(2*N),(N-j+0.5)*F/N)
    for j in range(1,N):
        print >> writer, 'brick x %11.8f y %11.8f z 200'                   % ((N-j)*F/N,F/N) # volumes N+4 to 2N+2: bricks right of the steps
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 100' % (3+N+j,(N+j)*F/(2*N),(N-j+0.5)*F/N)
    
    if geom == '2j':                                                                # unite the volumes so that there is one volume to the left of the steps
        if N != 1:                                                                  # and one volume to the right
            print >> writer, 'unite body       %s' % (' '.join(map(str,range(4,N+4))))
            print >> writer, 'unite body 1 2 3 %s' % (' '.join(map(str,range(N+4,2*N+3))))
        else:
            print >> writer, 'unite body 1 2 3'
        
    print >> writer, 'brick x 110 y 110 z 210'                                      # volume 2N+3: inner graveyard brick
    print >> writer, 'move volume %u location 50 50 100'            % (2*N+3)
    print >> writer, 'brick x 120 y 120 z 220'                                      # volume 2N+4: outer graveyard brick
    print >> writer, 'move volume %u location 50 50 100'            % (2*N+4)
    print >> writer, 'subtract %u from %u'                          % (2*N+3,2*N+4) # volume 2N+5: subtract inner from outer
    print >> writer, 'group "graveyard" add volume %u'              % (2*N+5)       # add the graveyard volume to the graveyard group
    if rho != 0:
        print >> writer, 'group "mat_1_rho_-%.8f" add volume all'   % (rho)         # add all volumes to the group representing material 1 and the specified
        print >> writer, 'group "mat_1_rho_-%.8f" remove volume %u' % (rho,2*N+5)   # rho; remove the graveyard volume from this group
    print >> writer, 'merge all'                                                    # merge all
    print >> writer, 'imprint body all'                                             # imprint all
    print >> writer, 'set geometry version 1902'
    print >> writer, 'set attribute on'
    print >> writer, 'export acis "%s%s.sat"'                       % (direc,fname) # export to .sat file
    
    writer.close()

# Modify LCAD files to change density
def modify_lcad(fname_in,fname_out,rho):
    
    reader = open(direc+fname_in+'.lcad','r')                                   # read from the input file
    writer = open(direc+fname_out+'.lcad','w')                                  # write to the output file

    lines = reader.readlines()
    in_str = lines[0].split()[2]                                                # input rho string
    out_str   = '-%10.8f' % (rho)                                               # output rho string
    
    line_ind = 0
    for line in lines:                                                          # replace original rho with new rho
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

# Write script to create runtpe files for native runs and do CUBIT and DAGMC-preprocessing and create runtpe files for DAG runs
def write_setup(fname,fname_in,version):
    
    writer = open(direc+fname+'.setup.sh','w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    if version == 'nat':
        print >> writer, 'mcnp5.mpi ix n=%s.i o=%s.setup.io'                                    % (fname,fname)
    elif version == 'dag':
        if fname == fname_in:
            print >> writer, 'cubit -batch -nographics -nojournal -information=off %s.jou'          % (fname)
            print >> writer, 'dagmc_preproc -f 1.0e-4 %s.sat -o %s.h5m'                             % (fname,fname)
            print >> writer, 'mcnp5.mpi ix n=%s.i g=%s.h5m o=%s.setup.io l=%s.lcad f=%s.setup.fcad' % (fname,fname,fname,fname,fname)
        else:
            print >> writer, 'mcnp5.mpi ix n=%s.i g=%s.h5m o=%s.setup.io l=%s.lcad f=%s.setup.fcad' % (fname,fname_in,fname,fname,fname)
    print >> writer, 'cp %s.ir %s.setup.ir'                                                         % (fname,fname)
        
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
    
    print >> writer, 'echo %s'                             % (fname)
    if (fid % max_concurrent) == 0:
        print >> writer, 'sh %s.setup.sh > %s.setup.out'   % (fname,fname)      # Create runtpe file, wait for completion
    else:
        print >> writer, 'sh %s.setup.sh > %s.setup.out &' % (fname,fname)      # Create runtpe file, do not wait for completion
    
    writer.close()


# Write local run script
def write_local_run(fname,version):
    
    writer = open(direc+'local_runs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s'                                                                                             % (fname)
    if version == 'nat':
        print >> writer, 'mpirun -np 22 mcnp5.mpi c i=cont.i o=%s.io r=%s.ir > %s.out'                                     % (fname,fname,fname)
    elif version == 'dag':
        print >> writer, 'mpirun -np 22 mcnp5.mpi c i=cont.i g=%s.setup.fcad o=%s.io r=%s.ir l=%s.lcad f=%s.fcad > %s.out' % (fname,fname,fname,fname,fname,fname)

    writer.close()

# Write ACI job script
def write_job(fname,version):
    
    time_extra = 1
    
    writer = open(direc+fname+'.aci.sh','w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    print >> writer, '#SBATCH --partition=univ'                                 # default "univ" if not specified
    
    print >> writer, '#SBATCH --time=0-00:%u:00' % (time_extra+ctme)            # run time in days-hh:mm:ss
    print >> writer, '#SBATCH --ntasks=16'                                      # number of CPUs
    print >> writer, '#SBATCH --mem-per-cpu=1024'                               # RAM in MB (default 4GB, max 8GB)
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err'  % (fname)       # location of error file
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out' % (fname)       # location of command output file
    print >> writer, ''
    if version == 'nat':
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi c i=cont.i o=%s.io r=%s.ir' % (fname,fname)
    elif version == 'dag':
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi c i=cont.i g=%s.setup.fcad o=%s.io r=%s.ir l=%s.lcad f=%s.fcad' % (fname,fname,fname,fname,fname)
    print >> writer, ''
    writer.close()
    
# Write master job script
def write_job_master(fname):
    
    writer = open(direc+'submit_jobs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s.sh'       % (fname)
    print >> writer, 'sbatch %s.aci.sh' % (fname)                               # execute the job script
    
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

reader.close()

min_N        =      1
max_N_nat_2s =  40000
max_N_nat_2j =   2000
max_N_dag    =  40000
min_F        =      0.1
max_F        =     99.9

fid = 0

# Write input files and job scripts
for params in list(itertools.product(versions,geoms,N_vals,F_vals,mfps)):
    
    fid = fid+1
    
    version = params[0]                                                         # loop for all MCNP versions
    geom    = params[1]                                                         # loop for all geometry configurations
    N       = params[2]                                                         # loop for all values of N
    F       = params[3]                                                         # loop for all values of F
    mfp     = params[4]                                                         # loop for all values of mfp
    
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
        max_concurrent = len(mfps)
    elif version == 'dag':
        max_concurrent = len(F_vals)
        
    fname = determine_filename(version,geom,N,F,mfp)                            # base filename (no extension)
    print fname
    fname_in = determine_filename(version,geom,N,F,mfp)
    
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
            fname_in = determine_filename(version,geom,N,F,mfp_in)
            modify_lcad(fname_in,fname,rho)
        else:
            write_cubit_2(fname,geom,N,F,rho)                                   # Write CUBIT instructions
        
    write_data_2(fname,version,geom,N,F,rho,ctme,mfp_in)                        # Write MCNP data cards

    write_setup(fname,fname_in,version)                                         # Write script to run CUBIT and DAGMC pre-processing and/or create runtpe file
    write_setup_master(fname,fname_in,fid,max_concurrent)                       # Append instructions to setup file
    write_local_run(fname,version)                                              # Append instructions to local run file
    write_job(fname,version)                                                    # Write ACI job file
    write_job_master(fname)                                                     # Append instructions to master job file
    
# Write continue run script
write_continue(ctme)

# Make all scripts executable
call('chmod 740 '+direc+'*.sh',shell=True)
