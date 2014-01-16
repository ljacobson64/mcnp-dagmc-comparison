"""

Automated MCNP5 Input File Generator

Author: Lucas Jacobson

Advisors: P.P.H. Wilson, Andrew Davis

The file "params.txt" should contain the following entries:
Line 1: Directory to place input files. Must include final slash.
Line 2: Directory to place results files. Must include final slash.
Line 3: List of versions. "nat" represents native MCNP5 and "dag" represents DAG-MCNP5.
Line 4: List of values of H (cylinder height.)
Line 5: List of values of D (cylinder diameter.)
Line 6: List of values of mfp (mfp is the number of mean-free paths of pure deuterium in
        one meter. The required density is obtained simply by multiplying mfp by 0.0078958.
Line 7: List of values of tol (faceting tolerance.)
Line 8: Computer time in minutes.

"""

from subprocess import call
import itertools
import math
import time

# Determine filename
def determine_filename(version,H,D,mfp,tol):
    
    fname = 'zCyl_%s_%.0f_%.0f_%.0e_%.0e' % (version,H,D,mfp,tol)               # base filename (no extension)
    return fname

# Write title cards
def write_title(fname,H,D,mfp,tol):
    
    # Write title cards to the specified file
    writer = open(direc+fname+'.i','w')
    
    print >> writer, 'Small geometric features test'
    print >> writer, 'c    H = %.3f' % (H)                                      # cylinder height
    print >> writer, 'c    D = %.3f' % (D)                                      # cylinder diameter
    print >> writer, 'c  rho = %.8f' % (rho)                                    # mass density of deuterium in the cube
    print >> writer, 'c  tol = %.3e' % (tol)                                    # faceting tolerance
    print >> writer, 'c'
    
    writer.close()

# Write cell cards
def write_cells(fname,H,D,rho,tol):
    
    # Determine material specification syntax
    if rho == 0:
        mat_str = '0            '                                               # if rho=0, set material to void
    else:
        mat_str = '1 -%10.8f' % (rho)                                           # otherwise, set material number to 1 with specified rho
    
    # Write cell cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, ''
    
    writer.close()

# Write surface cards
def write_surfaces(fname,H,D,rho,tol):
    
    # Write surface cards to the specified file
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c SURFACE CARDS'
    print >> writer, ''
    
    writer.close()

# Write data cards
def write_data(fname,version,H,D,rho,tol,ctme):
    
    writer = open(direc+fname+'.i','a')
    
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'                                                   # neutron transport problem
    print >> writer, 'IMP:N 0 1'                                                # neutron importance of 1 in the cylinder, 0 outside
    if version == 'nat':
        print >> writer, 'SDEF ERG=0.0000000253 X=0 Y=0 Z=%.8f' % (H/2)         # isotropic source of neutrons with an energy of 25.3 meV
    elif version == 'dag':
        print >> writer, 'SDEF ERG=0.0000000253 X=0 Y=0 Z=%.8f CELL=1' % (H/2)  # isotropic source of neutrons with an energy of 25.3 meV
    print >> writer, 'M1 1002.70C 1'                                            # material 1 is pure deuterium
    print >> writer, 'CTME %.8f' % (ctme)                                       # do not run simulation; only create runtpe file
    print >> writer, ''
    
    writer.close()

# Write continue run script
def write_continue(ctme):
    
    writer = open(direc+'cont.i','a')
    
    print >> writer, 'CONTINUE'
    print >> writer, 'CTME %.8f' % (ctme)
    
    writer.close()

# Write CUBIT commands
def write_cubit(fname,H,D,rho,tol):
    
    # Write CUBIT commands to the specified file
    writer = open(direc+fname+'.jou','w')
    
    print >> writer, 'merge all'                                                # merge all
    print >> writer, 'imprint body all'                                         # imprint all
    print >> writer, 'set geometry version 1902'
    print >> writer, 'set attribute on'
    print >> writer, 'export acis "%s%s.sat"' % (direc,fname)                   # export to .sat file
    
    writer.close()

# Write script to create runtpe files for native runs and CUBIT and DAGMC-preprocessing scripts for DAG runs
def write_setup(fname,version,H,D,rho,tol,ctme):
    
    writer = open(direc+fname+'.setup.sh','w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    if version == 'nat':
        print >> writer, 'mcnp5 ix n=%s.i o=%s.setup.io'                               % (fname,fname)
    elif version == 'dag':
        print >> writer, 'cubit -batch -nographics -nojournal -information=off %s.jou' % (fname)
        print >> writer, 'dagmc_preproc -f 1.0e-4 %s.sat -o %s.h5m'                    % (fname,fname)
        
    writer.close()

# Write script to run all setup scripts
def write_setup_master(fname,fid,max_concurrent):
    
    writer = open(direc+'setup_files.sh','a')
    
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
def write_local_run(fname,version,H,D,rho,tol,ctme):
    
    writer = open(direc+'local_runs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s'                                                % (fname)
    if   version == 'nat':
        print >> writer, 'mcnp5 c i=cont.i o=%s.io r=%s.ir > %s.out'          % (fname,fname,fname)
    elif version == 'dag':
        print >> writer, 'mcnp5 n=%s.i g=%s.h5m l=%s.lcad f=%s.fcad > %s.out' % (fname,fname,fname,fname,fname)
    
    writer.close()

# Write ACI job script
def write_job(fname,version,H,D,rho,tol,ctme):
    
    time_runtpe = 1
    
    writer = open(direc+fname+'.aci.sh','w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    print >> writer, '#SBATCH --partition=univ'                                 # default "univ" if not specified
    
    print >> writer, '#SBATCH --time=0-00:%u:00' % (time_runtpe+ctme)           # run time in days-hh:mm:ss
    print >> writer, '#SBATCH --ntasks=1'                                       # number of CPUs
    print >> writer, '#SBATCH --mem-per-cpu=2000'                               # RAM in MB (default 4GB, max 8GB)
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err'  % (fname)       # location of error file
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out' % (fname)       # location of command output file
    print >> writer, ''
    if   version == 'nat':                                                      # run MCNP
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi c i=cont.i o=%s.io r=%s.ir'          % (fname,fname)
    elif version == 'dag':
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i g=%s.h5m l=%s.lcad f=%s.fcad' % (fname,fname,fname,fname)
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
H_vals   = [  int(i) for i in params_str[3].split()[1:]]                        # list of values of N
D_vals   = [float(i) for i in params_str[4].split()[1:]]                        # list of values of F
mfps     = [float(i) for i in params_str[5].split()[1:]]                        # list of number of mean free paths in 1 m
tols     = [float(i) for i in params_str[6].split()[1:]]                        # list of faceting tolerances
ctme     =  float(            params_str[7].split()[1 ])                        # computer time

reader.close()

max_concurrent = len(tols)
fid = 0

# Write input files and job scripts
for params in list(itertools.product(versions,H_vals,D_vals,mfps,tols)):
    
    fid = fid+1
    
    version = params[0]                                                         # loop for all MCNP versions
    H       = params[1]                                                         # loop for all values of H
    D       = params[2]                                                         # loop for all values of D
    mfp     = params[3]                                                         # loop for all values of mfp
    tol     = params[4]                                                         # loop for all values of tol
    
    # Calculate mass density of deuterium from input number of mean free paths in 1 meter
    # Setting the mean free path to 1 meter results in a density of 0.0078958 g/cm^3
    rho = mfp*0.0078958
    
    fname = determine_filename(version,H,D,mfp,tol)                             # base filename (no extension)
    print fname
    
    # Write input files
    if version == 'nat':                                                        # Native MCNP
        
        write_title(fname,H,D,mfp,tol)                                          # Write title cards
        write_cells(fname,H,D,mfp,tol)                                          # Append cell cards; joined cells
        write_surfaces(fname,H,D,mfp,tol)                                       # Append surface cards
        
    elif version == 'dag':                                                      # DAG-MCNP
        
        write_cubit_2(fname,H,D,mfp,tol)                                        # Write CUBIT instructions
    
    write_data(fname,version,H,D,mfp,tol,ctme)                                  # Write MCNP data cards

    write_setup(fname,version,H,D,mfp,tol,ctme)                                 # Write script to create runtpe file or run CUBIT and DAGMC pre-processing
    write_setup_master(fname,fid,max_concurrent)                                # Append instructions to setup file

    write_local_run(fname,version,H,D,mfp,tol,ctme)                             # Append instructions to local run file
    write_job(fname,version,H,D,mfp,tol,ctme)                                   # Write ACI job file
    write_job_master(fname)                                                     # Append instructions to master job file
    
# Write continue run script
write_continue(ctme)

# Make all scripts executable
call('chmod 740 '+direc+'*.sh',shell=True)
