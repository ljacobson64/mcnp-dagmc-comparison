from subprocess import call
import itertools
import math
import time

# Write title cards
def write_title(fname,geom,N,F,rho,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'w')
    
    print >> writer, 'Small geometric features test'
    print >> writer, 'c    N = %u'   % (N)
    print >> writer, 'c    F = %.3f' % (F)
    print >> writer, 'c  rho = %.8f' % (rho)
    print >> writer, 'c geom = %s'   % (geom)
    print >> writer, 'c'
    
    writer.close()

# Write cell cards; 2 dimensions; separate cells
def write_cells_2s(fname,geom,N,F,rho,ctme):
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0                    -09000'
    print >> writer, '   12 0              09999'
    print >> writer, '   13 0              09000 -09999        -01000'
    print >> writer, '   14 0              09000 -09999  04999'
    print >> writer, '   15 0              09000 -09999  01000 -04999        -05000'
    print >> writer, '   16 0              09000 -09999  01000 -04999  08999'
    print >> writer, '    1 1 -%.8f  09000 -09999  %u -04999  %u -08999' % (rho,10000+N,50000+N)
    print >> writer, '    2 1 -%.8f  09000 -09999  01000 -%u  %u -08999' % (rho,10000+N,50000+N)
    print >> writer, '    3 1 -%.8f  09000 -09999  %u -04999  05000 -%u' % (rho,10000+N,50000+N)
    for j in range(1,N):
        print >> writer, '%u 1 -%.8f  09000 -09999  %u -%u  %u -%u'      % (10000+j,rho,10000+j,10000+N,50000+N-j,50000+N-j+1)
        print >> writer, '%u 1 -%.8f  09000 -09999  01000 -%u  %u -%u'   % (50000+j,rho,10000+j,50000+N-j,50000+N-j+1)
    print >> writer, '%u 1 -%.8f  09000 -09999  01000 -%u  05000 -50001' % (50000+N,rho,10000+N)
    print >> writer, ''
    
    writer.close()

# Write cell cards; 2 dimensions; joined cells
def write_cells_2j(fname,geom,N,F,rho,ctme):
	
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c CELL CARDS'
    print >> writer, '   11 0              09999:-09000: 04999:-01000: 08999:-05000'
    print >> writer, '    1 1 -%.8f  09000 -09999 -04999 -08999 (' % (rho)
    print >> writer, '                     01000  %u:'             % (50000+N)
    for j in range(1,N):
        print >> writer, '                     %u  %u:'            % (10000+j,50000+N-j)
    print >> writer, '                     %u  05000)'             % (10000+N)
    print >> writer, '    2 1 -%.8f  09000 -09999  01000  05000 (' % (rho)
    for j in range(1,N):
        print >> writer, '                    -%u -%u:'            % (10000+j,50000+N-j+1)
    print >> writer, '                    -%u -%u)'                % (10000+N,50000+1)
    print >> writer, ''
    
    writer.close()

# Write surface cards; 2 dimensions
def write_surfaces_2(fname,geom,N,F,rho,ctme):
    
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

# Write data cards; native MCNP; 2 dimensions
def write_data_nat_2(fname,geom,N,F,rho,ctme):
    
    x_src = (200-F)/(2+math.sqrt(2))
    y_src = (200-F)/(2+math.sqrt(2))
    z_src = 50
    
    fname_i = '%s.i' % (fname)
    writer = open(direc+fname_i,'a')
    
    print >> writer, 'c DATA CARDS'
    print >> writer, 'MODE N'
    if geom == '2s':
        print >> writer, 'IMP:N 0 5R 1 %uR'                       % (2*N+1)
    elif geom == '2j':
        print >> writer, 'IMP:N 0 1 1'
    print >> writer, 'SDEF ERG=0.0000000253 X=%.8f Y=%.8f Z=%.8f' % (x_src,y_src,z_src)
    print >> writer, 'CTME %.8f'                                  % (ctme)
    print >> writer, 'M1 1002.70C 1'
    
    writer.close()

# Write data cards; DAG-MCNP; 2 dimensions
def write_data_dag_2(fname,geom,N,F,rho,ctme):
    
    x_src = 100-(200-F)/(2+math.sqrt(2))
    y_src = 100-(200-F)/(2+math.sqrt(2))
    z_src = 50
    
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
    
    print >> writer, 'DATA CARDS'
    print >> writer, 'MODE N'
    print >> writer, 'SDEF ERG=0.0000000253 X=%.8f Y=%.8f Z=%.8f CELL=%u' % (x_src,y_src,z_src,cell_src)
    print >> writer, 'DAGMC CHECK_SRC_CELL=OFF'
    print >> writer, 'CTME %.8f' % (ctme)
    print >> writer, 'M1 1002.70C 1'
    
    writer.close()

# Write CUBIT commands; 2 dimensions
def write_cubit_2(fname,geom,N,F,rho,ctme):
    
    fname_jou = '%s.jou' % (fname)
    writer = open(direc+fname_jou,'w')
    
    for j in range(1,N+1):
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (j*F/N,F/N)
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (j,F*j/(2*N),(N-j+0.5)*F/N)   
    for j in range(1,N):
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % ((N-j)*F/N,F/N)
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (N+j,(N+j)*F/(2*N),(N-j+0.5)*F/N)   
    
    if F != 100:
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (100-F,100-F)
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (2*N,(100+F)/2,(100+F)/2)  
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (F,100-F)
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (2*N+1,F/2,(100+F)/2)
        print >> writer, 'brick x %11.8f y %11.8f z 100'                  % (100-F,F)
        print >> writer, 'move volume %u location x %11.8f y %11.8f z 50' % (2*N+2,(100+F)/2,F/2)
        g_ind = 2*N+3
    else:
        g_ind = 2*N
    
    if geom == '2j':
        if N == 1:
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
    
    print >> writer, 'brick x 110 y 110 z 110'
    print >> writer, 'move volume %u location 50 50 50'         % (g_ind)
    print >> writer, 'brick x 120 y 120 z 120'
    print >> writer, 'move volume %u location 50 50 50'         % (g_ind+1)
    print >> writer, 'subtract %u from %u'                      % (g_ind,g_ind+1)
    print >> writer, 'group "graveyard" add volume %u'          % (g_ind+2)
    print >> writer, 'group "mat_1_rho_-%.8f" add volume all'   % (rho)
    print >> writer, 'group "mat_1_rho_-%.8f" remove volume %u' % (rho,g_ind+2)
    print >> writer, 'merge all'
    print >> writer, 'imprint body all'
    print >> writer, 'set geometry version 1902'
    print >> writer, 'set attribute on'
    # print >> writer, 'save as "%s%s.cub"'                       % (direc,fname)
    print >> writer, 'export acis "%s%s.sat"'                   % (direc,fname)
    
    writer.close()

# Write local run script
def write_local_run(fname,version):
    
    writer = open(direc+'local_runs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    if version == 'nat':
        print >> writer, 'mcnp5 n=%s.i'          % (fname)
    elif version == 'dag':
        print >> writer, 'mcnp5 n=%s.i g=%s.h5m' % (fname,fname)
    
    writer.close()

# Write HPC job script
def write_job(fname,version,geom,N,F,rho,ctme):
    
    fname_sh = '%s.sh' % (fname)
    writer = open(direc+fname_sh,'w')
    
    print >> writer, '#!/bin/bash'
    print >> writer, '#SBATCH --partition=univ'                                 # default "univ" if not specified
    print >> writer, '#SBATCH --time=0-00:06:00'                                # run time in days-hh:mm:ss
    print >> writer, '#SBATCH --ntasks=1'                                       # number of CPUs
    print >> writer, '#SBATCH --mem-per-cpu=2000'                               # RAM in MB (default 4GB, max 8GB)
    print >> writer, '#SBATCH --error=/home/ljjacobson/%s.err'                                  % (fname)
    print >> writer, '#SBATCH --output=/home/ljjacobson/%s.out'                                 % (fname)
    print >> writer, ''
    if version == 'nat':
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i'          % (fname)
    elif version == 'dag':
        print >> writer, 'srun /home/adavis23/dagmc/mcnp5v16_dag/bin/mcnp5.mpi n=%s.i g=%s.h5m' % (fname,fname)
    writer.close()
    
# Write master job script
def write_master_job(fname):
    
    writer = open(direc+'submit_jobs.sh','a')
    
    if writer.tell() == 0:
        print >> writer, '#!/bin/bash'
        print >> writer, ''
    
    print >> writer, 'echo %s.sh'   % (fname)
    print >> writer, 'sbatch %s.sh' % (fname)
    
    writer.close()

# MAIN SCRIPT

# Open file containing parameters to vary
reader = open('params.txt','r')
params_str = reader.readlines()

global direc
direc      =                    params_str[0].split()[0 ]                       # directory to place files
versions   =                    params_str[1].split()[1:]                       # list of MCNP versions
geoms      =                    params_str[2].split()[1:]                       # list of geometry configurations
N_vals     = [  int(i) for i in params_str[3].split()[1:]]                      # list of values of N
F_vals     = [float(i) for i in params_str[4].split()[1:]]                      # list of values of F
mfps       = [float(i) for i in params_str[5].split()[1:]]                      # list of number of mean free paths in 1 m
ctme       =  float(            params_str[6].split()[1 ])                      # computer time

reader.close()

min_N        =      1
max_N_nat_2s =  40000
max_N_nat_2j =   1000
max_N_dag    =  40000
min_F_nat    =      0
min_F_dag    =      0.001
max_F        =    100

# Write input files and job scripts
for params in list(itertools.product(versions,geoms,N_vals,F_vals,mfps)):
    version = params[0]                                                         # loop for all MCNP versions
    geom    = params[1]                                                         # loop for all geometry configurations
    N       = params[2]                                                         # loop for all values of N
    F       = params[3]                                                         # loop for all values of F
    mfp     = params[4]                                                         # loop for all values of mfp
    
    start_time = time.time()
    
    # Calculate mass density of deuterium from input number of mean free paths in 1 meter
    # Setting the mean free path to 1 meter results in a density of 0.0078958 g/cm^3
    rho = mfp*0.0078958
    
    print '|==============================================================================|'
    print '|   version = %s    geom = %s    N = %5u    F = %7.3f    mfp = %8.3f   |' % (version,geom,N,F,mfp)
    print '|==============================================================================|'
    
    # Determine whether parameters are valid
    
    # N must be between the minimum value and maximum value
    if   version == 'nat' and geom == '2s':
        valid_N = N >= min_N and N <= max_N_nat_2s
    elif version == 'nat' and geom == '2j':
        valid_N = N >= min_N and N <= max_N_nat_2j
    elif version == 'dag' and (geom == '2s' or geom == '2j'):
        valid_N = N >= min_N and N <= max_N_dag
    else:
        valid_N = False
    
    # F must be between the minimum value and maximum value
    if   version == 'nat' and (geom == '2s' or geom == '2j'):
        valid_F = F >= min_F_nat and F <= max_F
    elif version == 'dag' and (geom == '2s' or geom == '2j'):
        valid_F = F >= min_F_dag and F <= max_F
    else:
        valid_F = False

    # mfp must be greater than or equal to zero
    valid_mfp = mfp >= 0
    
    if not valid_N or not valid_F or not valid_mfp:
        print 'Invalid parameters'
        continue
    
    fname = 'zCube_%s_%s_%u_%.0f_%.2f' % (version,geom,N,F,mfp)                 # base filename (no extension)
    
    # Write input files
    if version == 'nat' and (geom == '2s' or geom == '2j'):                     # Native MCNP, 2 dimensions
        
        write_title(fname,geom,N,F,rho,ctme)                                    # Write title cards
        if   geom == '2s':
            write_cells_2s(fname,geom,N,F,rho,ctme)                             # Append cell cards; separate cells
        elif geom == '2j':
            write_cells_2j(fname,geom,N,F,rho,ctme)                             # Append cell cards; joined cells
        write_surfaces_2(fname,geom,N,F,rho,ctme)                               # Append surface cards
        write_data_nat_2(fname,geom,N,F,rho,ctme)                               # Append data cards
        
        write_job(fname,version,geom,N,F,rho,ctme)                              # Write HPC job file
        write_master_job(fname)                                                 # Append instructions to master job file
        write_local_run(fname,version)                                          # Append instructions to local run file
        
    if version == 'dag' and (geom == '2s' or geom == '2j'):                     # DAG-MCNP, 2 dimensions
        
        write_data_dag_2(fname,geom,N,F,rho,ctme)                               # Write MCNP data cards for DAG
        write_cubit_2(fname,geom,N,F,rho,ctme)                                  # Write CUBIT instructions; separate cells
        call('cubit -batch -nographics -nojournal -information=off '+direc+fname+'.jou',shell=True)
                                                                                # Run CUBIT from command line
        call('dagmc_preproc -f 1.0e-4 '+direc+fname+'.sat -o '+direc+fname+'.h5m',shell=True)
                                                                                # DAGMC pre-processing
        
        write_job(fname,version,geom,N,F,rho,ctme)                              # Write HPC job file
        write_master_job(fname)                                                 # Append instructions to master job file
        write_local_run(fname,version)                                          # Append instructions to local run file
    
    # Record time taken
    write_time = time.time()-start_time
    print '%.3f seconds' % (write_time)
    writer = open(direc+'timing.txt','a')
    print >> writer, '%s %s %5u %7.3f %7.3f %16.8f' % (version,geom,N,F,mfp,write_time)
    writer.close()

# Make master job script executable
call('chmod 770 '+direc+'submit_jobs.sh',shell=True)