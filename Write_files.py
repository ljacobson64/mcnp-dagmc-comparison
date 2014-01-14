"""

Automated MCNP5 Input File Generator

Author: Lucas Jacobson

Advisors: P.P.H. Wilson, Andrew Davis

The file "params.txt" should contain the following entries:
Line 1: Directory. Must include final slash.
Line 2: List of versions. "nat" represents native MCNP5 and "dag" represents DAG-MCNP5.
Line 3: List of geometry configurations. "2s" represents two-dimensional steps with separate cell
        definitions; "2j" represents two-dimensional steps with joined cell definitions. "3s" and
        "3j" may become options in the future.
Line 4: List of values of N (integer from 1 to 40000). N is the number of steps in the geometry.
Line 5: List of values of F (0.001 to 100). F is the portion of the geometry containing the steps.
        For example, an F of 0.001 would mean that all the steps would be condensed into a very
        small space in a corner, whereas an F of 100 would mean that the steps would comprise an
        entire square diagonal.
Line 6: List of values of mfp (0 to 1000). mfp is the number of mean-free paths of pure deuterium
        in one meter. The required density is obtained simply by multiplying mfp by 0.0078958.
Line 7: Computer time in minutes. Most runs will take a total of ctme+1 minutes to run, but some
        will take an additional few minutes to create the runtpe file. The most extreme cases
        require over a day of extra time to account for runtpe file creation.
Line 8: If running in DAG-MCNP mode and mfp_in > 0, CUBIT and dagmc_preproc are not run. Rather, an
        existing .h5m file with identical geometry but different density is used, and the only
        things changed are the density entries on the LCAD card.

"""

from subprocess import call
import itertools
import math
import time

# Determine filename
def determine_filename(version,geom,N,F,mfp):
    
    fname = 'zCube_%s_%s_%u_%.0f_%.3f' % (version,geom,N,F,mfp)                 # base filename (no extension)
    return fname

# Write title cards
def write_title(fname,geom,N,F,rho):
    
    # Write title cards to the specified file
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'w')
    
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
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0                    -09000'                        # rest of universe
    print >> writer, '   12 0              09999'
    print >> writer, '   13 0              09000 -09999        -01000'
    print >> writer, '   14 0              09000 -09999  04999'
    print >> writer, '   15 0              09000 -09999  01000 -04999        -05000'
    print >> writer, '   16 0              09000 -09999  01000 -04999  08999'
    print >> writer, '    1 %s  09000 -09999  %u -04999  %u -08999' % (mat_str,10000+N,50000+N) # top-right portion of cube
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
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0              09999:-09000: 04999:-01000: 08999:-05000' # rest of universe
    print >> writer, '    1 %s  09000 -09999 -04999 -08999 ('      % (mat_str)  # bottom-left cell
    print >> writer, '                     01000  %u:'             % (50000+N)
    for j in range(1,N):
        print >> writer, '                     %u  %u:'            % (10000+j,50000+N-j)
    print >> writer, '                     %u  05000)'             % (10000+N)
    print >> writer, '    2 %s  09000 -09999  01000  05000 ('      % (mat_str)  # top-right cell
    for j in range(1,N):
        print >> writer, '                    -%u -%u:'            % (10000+j,50000+N-j+1)
    print >> writer, '                    -%u -%u)'                % (10000+N,50000+1)
    print >> writer, ''
    
    writer.close()

# Write surface cards; 2 dimensions
def write_surfaces_2(fname,geom,N,F,rho):
    
    # Write surface cards to the specified file
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c SURFACE CARDS'
    print >> writer, '01000 PX 0'                                               # left boundary of cube
    print >> writer, '05000 PY 0'                                               # front boundary of cube
    print >> writer, '09000 PZ 0'                                               # bottom boundary of cube
    print >> writer, '04999 PX 100'                                             # right boundary of cube
    print >> writer, '08999 PY 100'                                             # back boundary of cube
    print >> writer, '09999 PZ 100'                                             # top boundary of cube
    for j in range(1,N+1):
        print >> writer, '%u PX %11.8f' % (10000+j,F*j/N)                       # planes normal to the x-axis for the N steps
        print >> writer, '%u PY %11.8f' % (50000+j,F*j/N)                       # planes normal to the y-axis for the N steps
    print >> writer, ''
    
    writer.close()

# Write data cards; native MCNP; 2 dimensions
def write_data_nat_2(fname,geom,N,F,rho,ctme,mfp_in):
    
    # Determine the source point    
    x_src = 100-(200-F)/(2+math.sqrt(2))                                        # neutron source point in the x- and y-directions is the center of the circle
    y_src = 100-(200-F)/(2+math.sqrt(2))                                        # defined such that it is tangent to the top and right side of the cube and
                                                                                # tangent to the diagonal line where the N steps lie
    z_src = 50                                                                  # neutron source point in the z-direction is in the center of the cube
    
    # Write data cards to the specified file
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'                                                   # neutron transport problem
    if geom == '2s':                                                            # neutron importance of 1 in the cube, 0 in the rest of universe
        print >> writer, 'IMP:N 0 5R 1 %uR'                       % (2*N+1)
    elif geom == '2j':
        print >> writer, 'IMP:N 0 1 1'
    print >> writer, 'SDEF ERG=0.0000000253 X=%.8f Y=%.8f Z=%.8f' % (x_src,y_src,z_src) # isotropic source of neutrons with an energy of 25.3 meV
    print >> writer, 'CTME %.8f'                                  % (ctme)      # amount of computer time to run the simulation
    print >> writer, 'M1 1002.70C 1'                                            # material 1 is pure deuterium
    print >> writer, ''
    
    writer.close()

# Write data cards; DAG-MCNP; 2 dimensions
def write_data_dag_2(fname,geom,N,F,rho,ctme,mfp_in):
    
    # Determine the source point
    x_src = 100-(200-F)/(2+math.sqrt(2))                                        # neutron source point in the x- and y-directions is the center of the circle
    y_src = 100-(200-F)/(2+math.sqrt(2))                                        # defined such that it is tangent to the top and right side of the cube and
                                                                                # tangent to the diagonal line where the N steps lie
    z_src = 50                                                                  # neutron source point in the z-direction is in the center of the cube
    
    # Determine in which cell the source exists
    if geom == '2s':
        if N == 1 and F == 100:
            cell_src = 1
        elif x_src > F:
            cell_src = 2*N
        else:
            cell_src = 2*N-int(y_src*N/F)
    elif geom == '2j':
        if N == 1 and F == 100:
            cell_src = 1
        else:
            cell_src = N+1
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    # Write data cards to the specified file
    print >> writer, 'DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'SDEF ERG=0.0000000253 X=%.8f Y=%.8f Z=%.8f CELL=%u' % (x_src,y_src,z_src,cell_src)
    print >> writer, 'DAGMC CHECK_SRC_CELL=OFF'
    print >> writer, 'CTME %.8f' % (ctme)
    print >> writer, 'M1 1002.70C 1'
    print >> writer, ''
    
    writer.close()

# Write CUBIT commands; 2 dimensions
def write_cubit_2(fname,geom,N,F,rho):
    
    # Write CUBIT commands to the specified file
    fname_jou = '%s.jou' % (fname)
    writer = open(direc+fname_jou,'w')
    
    for j in range(1,N+1):
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (j*F/N,F/N)     # create bricks left of the N steps
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (j,F*j/(2*N),(N-j+0.5)*F/N)
    for j in range(1,N):
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % ((N-j)*F/N,F/N) # create bricks right of the N steps
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (N+j,(N+j)*F/(2*N),(N-j+0.5)*F/N)
    
    if F != 100:
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (100-F,100-F)   # create top-right brick
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (2*N,(100+F)/2,(100+F)/2)
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (F,100-F)       # create top-left brick
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (2*N+1,F/2,(100+F)/2)
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (100-F,F)       # create bottom-right brick
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (2*N+2,(100+F)/2,F/2)
        g_ind = 2*N+3
    else:
        g_ind = 2*N
    
    if geom == '2j':                                                            # unite the boxes so that there is one volume to the left of the steps and
        if N == 1:                                                              # one volume to the right
            if F != 100:
                print >> writer, 'unite body %s' % (' '.join(map(str,range(N+1,2*N+3))))
        elif N == 2:
            print >> writer, 'unite body %s'     % (' '.join(map(str,range(1,N+1))))
            if F != 100:
                print >> writer, 'unite body %s' % (' '.join(map(str,range(N+1,2*N+3))))
        else:
            print >> writer, 'unite body %s'     % (' '.join(map(str,range(1,N+1))))
            if F == 100:
                print >> writer, 'unite body %s' % (' '.join(map(str,range(N+1,2*N)))) 
            else:
                print >> writer, 'unite body %s' % (' '.join(map(str,range(N+1,2*N+3))))
        
        # print >> writer, 'unite body %u %u %u'      % (2*N,2*N+1,2*N+2)
        # tree_size = int(math.log(N,2))
        # for tree_level in range(1,tree_size+1):
        #     for j in range(1,N/(2**tree_level)+1):
        #         print >> writer, 'unite body %u %u' % (2**tree_level*j-(2**tree_level-1),2**tree_level*j-(2**(tree_level-1)-1))
        #         print >> writer, 'unite body %u %u' % (2**tree_level*j-(2**tree_level-1)+N,2**tree_level*j-(2**(tree_level-1)-1)+N)
        # tree_temp = N
        # end_tree = [1]
        # for tree_level in reversed(range(1,tree_size+1)):
        #     if tree_temp-2**tree_level > 0:
        #         tree_temp -= 2**tree_level
        #         end_tree.append(end_tree[-1]+2**tree_level)
        #     if tree_temp-2**tree_level == 0:
        #         break
        # for j in reversed(range(1,len(end_tree))):
        #     print >> writer, 'unite body %u %u'     % (end_tree[j-1],end_tree[j])
        #     print >> writer, 'unite body %u %u'     % (end_tree[j-1]+N,end_tree[j]+N)
    
    print >> writer, 'brick x 110 y 110 z 110'                                        # create inner graveyard brick
    print >> writer, 'move volume %u location 50 50 50'             % (g_ind)
    print >> writer, 'brick x 120 y 120 z 120'                                        # create outer graveyard brick
    print >> writer, 'move volume %u location 50 50 50'             % (g_ind+1)
    print >> writer, 'subtract %u from %u'                          % (g_ind,g_ind+1) # subtract inner from outer
    print >> writer, 'group "graveyard" add volume %u'              % (g_ind+2)       # add the graveyard volume to the graveyard group
    if rho != 0:
        print >> writer, 'group "mat_1_rho_-%.8f" add volume all'   % (rho)           # add all volumes to the group representing material 1 and the specified
        print >> writer, 'group "mat_1_rho_-%.8f" remove volume %u' % (rho,g_ind+2)   # rho; remove the graveyard volume from this group
    print >> writer, 'merge all'                                                      # merge all
    print >> writer, 'imprint body all'                                               # imprint all
    print >> writer, 'set geometry version 1902'
    print >> writer, 'set attribute on'
    # print >> writer, 'save as "%s%s.cub"'                           % (direc,fname)   # save as .cub file
    print >> writer, 'export acis "%s%s.sat"'                       % (direc,fname)   # export to .sat file
    
    writer.close()

# Modify LCAD files to change density
def modify_lcad(fname_in,fname_out,rho):
    
    fname_in  = '%s.lcad' % (fname_in)                                          # input filename
    fname_out = '%s.lcad' % (fname_out)                                         # output filename
    
    reader = open(direc+fname_in,'r')                                           # read from the input file
    writer = open(direc+fname_out,'w')                                          # write to the output file
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

# Write script to run CUBIT and DAGMC pre-processing
def write_cdp_run(fname):
    
    fname_cdp_sh = '%s.cdp.sh' % (fname)
    writer = open(direc+fname_cdp_sh,'w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    print >> writer, 'cubit -batch -nographics -nojournal -information=off %s.jou' % (fname)
    print >> writer, 'dagmc_preproc -f 1.0e-4 %s.sat -o %s.h5m'                    % (fname,fname)
    
    writer.close()

# Write script to run all CUBIT and DAGMC-preprocessing scripts in background
def write_cdp_master(fname,fid,max_concurrent):
    
    writer = open(direc+'cdp.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s'                        % (fname)
    if (fid % max_concurrent) == 0:
        print >> writer, 'sh %s.cdp.sh > %s.cdpout'   % (fname,fname)           # Run CUBIT+DAGMC_preproc script, wait for completion
    else:
        print >> writer, 'sh %s.cdp.sh > %s.cdpout &' % (fname,fname)           # Run CUBIT+DAGMC_preproc script, do not wait for completion
    
    writer.close()

# Write local run script
def write_local_run(fname,version,geom,N,F,rho,ctme,mfp_in):
    
    writer = open(direc+'local_runs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    # Write script to run CUBIT and DAGMC pre-processing
    if   version == 'nat':
        print >> writer, 'mcnp5 n=%s.i'                              % (fname)
    elif version == 'dag' and mfp_in > 0:
        fname_in = determine_filename(version,geom,N,F,mfp_in)
        print >> writer, 'mcnp5 n=%s.i g=%s.h5m l=%s.lcad f=%s.fcad' % (fname,fname_in,fname,fname)
    elif version == 'dag' and mfp_in <= 0:
        print >> writer, 'mcnp5 n=%s.i g=%s.h5m l=%s.lcad f=%s.fcad' % (fname,fname,fname,fname)
    
    writer.close()

# Write ACI job script
def write_job(fname,version,geom,N,F,rho,ctme,mfp_in):
    
    fname_aci_sh = '%s.sh' % (fname)
    writer = open(direc+fname_aci_sh,'w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, ''
    print >> writer, '#SBATCH --partition=univ'                                 # default "univ" if not specified
    
    # Determine how much time to add to ctme to account for the runtpe creation time
    if version == 'nat':
        if geom == '2s':
            if N <= 4000:
                time_runtpe = 1
            elif N <= 20000:
                time_runtpe = 10
            elif N <= 40000:
                time_runtpe = 60
        elif geom == '2j' and mfp <= 100:
            if N <= 100:
                time_runtpe = 1
            elif N <= 400:
                time_runtpe = 10
            elif N <= 1000:
                time_runtpe = 240
            elif N <= 2000:
                time_runtpe = 1200
        elif geom == '2j' and mfp <= 1000:
            if N <= 40:
                time_runtpe = 1
            elif N <= 200:
                time_runtpe = 5
            elif N <= 400:
                time_runtpe = 20
            elif N <= 1000:
                time_runtpe = 840
            elif N <= 2000:
                time_runtpe = 4200
    elif version == 'dag':
        time_runtpe = 5                                                         # placeholder; guess 5 minutes for now
    if 'time_runtpe' not in locals():
        time_runtpe = 5
    
    print >> writer, '#SBATCH --time=0-00:%u:00' % (time_runtpe+ctme)           # run time in days-hh:mm:ss
    print >> writer, '#SBATCH --ntasks=1'                                       # number of CPUs
    print >> writer, '#SBATCH --mem-per-cpu=2000'                               # RAM in MB (default 4GB, max 8GB)
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err'   % (fname)      # location of error file
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out'  % (fname)      # location of command output file
    print >> writer, ''
    if   version == 'nat':                                                      # run MCNP
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i'                              % (fname)
    elif version == 'dag' and mfp_in > 0:
        fname_in = determine_filename(version,geom,N,F,mfp_in)
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i g=%s.h5m l=%s.lcad f=%s.fcad' % (fname,fname_in,fname,fname)
    elif version == 'dag' and mfp_in <= 0:
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i g=%s.h5m l=%s.lcad f=%s.fcad' % (fname,fname,fname,fname)
    print >> writer, ''
    writer.close()
    
# Write master job script
def write_job_master(fname):
    
    writer = open(direc+'submit_jobs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s.sh'   % (fname)
    print >> writer, 'sbatch %s.sh' % (fname)                                   # execute the job script
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params_str = reader.readlines()

global direc
direc    =                    params_str[0].split()[0 ]                         # directory to place files
versions =                    params_str[1].split()[1:]                         # list of MCNP versions
geoms    =                    params_str[2].split()[1:]                         # list of geometry configurations
N_vals   = [  int(i) for i in params_str[3].split()[1:]]                        # list of values of N
F_vals   = [float(i) for i in params_str[4].split()[1:]]                        # list of values of F
mfps     = [float(i) for i in params_str[5].split()[1:]]                        # list of number of mean free paths in 1 m
ctme     =  float(            params_str[6].split()[1 ])                        # computer time
mfp_in   =  float(            params_str[7].split()[1 ])                        # mfp for pre-existing LCAD file

reader.close()

min_N        =      1
max_N_nat_2s =  40000
max_N_nat_2j =   2000
max_N_dag    =  40000
min_F_nat    =      0
min_F_dag    =      0.001
max_F        =    100

max_concurrent = len(F_vals)
fid = 0

# Write input files and job scripts
for params in list(itertools.product(versions,geoms,N_vals,F_vals,mfps)):
    
    fid = fid+1
    
    version = params[0]                                                         # loop for all MCNP versions
    geom    = params[1]                                                         # loop for all geometry configurations
    N       = params[2]                                                         # loop for all values of N
    F       = params[3]                                                         # loop for all values of F
    mfp     = params[4]                                                         # loop for all values of mfp
    
    start_time = time.time()
    
    # Calculate mass density of deuterium from input number of mean free paths in 1 meter
    # Setting the mean free path to 1 meter results in a density of 0.0078958 g/cm^3
    rho = mfp*0.0078958
    
    # print '|==============================================================================|'
    # print '|   version = %s    geom = %s    N = %5u    F = %7.3f    mfp = %8.3f   |' % (version,geom,N,F,mfp)
    # print '|==============================================================================|'
    
    # Determine whether parameters are valid
    
    valid_version = version == 'nat' or version == 'dag'                        # version must be nat or dag
    valid_geom = geom == '2s' or geom == '2j'                                   # geom must be 2s or 2j
    
    if   version == 'nat' and geom == '2s':                                     # N must be between the minimum value and maximum value
        valid_N = N >= min_N and N <= max_N_nat_2s
    elif version == 'nat' and geom == '2j':
        valid_N = N >= min_N and N <= max_N_nat_2j
    elif version == 'dag' and (geom == '2s' or geom == '2j'):
        valid_N = N >= min_N and N <= max_N_dag
    else:
        valid_N = False
    
    if   version == 'nat':                                                      # F must be between the minimum value and maximum value
        valid_F = F >= min_F_nat and F <= max_F
    elif version == 'dag':
        valid_F = F >= min_F_dag and F <= max_F
    else:
        valid_F = False

    valid_mfp = mfp >= 0                                                        # mfp must be greater than or equal to zero
    
    if not valid_N or not valid_F or not valid_mfp:
        print 'Invalid parameters'
        continue
    
    fname = determine_filename(version,geom,N,F,mfp)                            # base filename (no extension)
    
    # Write input files
    if   version == 'nat' and (geom == '2s' or geom == '2j'):                   # Native MCNP, 2 dimensions
        
        write_title(fname,geom,N,F,rho)                                         # Write title cards
        if   geom == '2s':
            write_cells_2s(fname,geom,N,F,rho)                                  # Append cell cards; separate cells
        elif geom == '2j':
            write_cells_2j(fname,geom,N,F,rho)                                  # Append cell cards; joined cells
        write_surfaces_2(fname,geom,N,F,rho)                                    # Append surface cards
        write_data_nat_2(fname,geom,N,F,rho,ctme,mfp_in)                        # Append data cards
        
    elif version == 'dag' and (geom == '2s' or geom == '2j'):                   # DAG-MCNP, 2 dimensions
        
        if mfp_in > 0:                                                          # Use an existing LCAD file with a different rho
            fname_in = 'zCube_%s_%s_%u_%.0f_%.2f' % (version,geom,N,F,mfp_in)
            modify_lcad(fname_in,fname,rho)
        else:
            write_cubit_2(fname,geom,N,F,rho)                                   # Write CUBIT instructions
            write_cdp_run(fname)                                                # Write script to run CUBIT and DAGMC pre-processing
            write_cdp_master(fname,fid,max_concurrent)                          # Append instructions to C+DP run file

        write_data_dag_2(fname,geom,N,F,rho,ctme,mfp_in)                        # Write MCNP data cards for DAG

    write_local_run(fname,version,geom,N,F,rho,ctme,mfp_in)                     # Append instructions to local run file
    write_job(fname,version,geom,N,F,rho,ctme,mfp_in)                           # Write ACI job file
    write_job_master(fname)                                                     # Append instructions to master job file
    
    # Record time taken
    write_time = time.time()-start_time
    print '%s %s %5u %7.3f %8.3f %16.8f' % (version,geom,N,F,mfp,write_time)
    writer = open(direc+'timing.txt','a')
    print >> writer, '%s %s %5u %7.3f %8.3f %16.8f' % (version,geom,N,F,mfp,write_time)
    writer.close()

# Make all scripts executable
call('chmod 770 '+direc+'*.sh',shell=True)
