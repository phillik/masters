#!/bin/bash
##########################
# example for an MPI job #
##########################

#SBATCH --job-name=test
#SBATCH --account=nn9264k

# 80 MPI tasks in total
# Stallo has 16 or 20 cores/node and therefore we take
# a number that is divisible by both

#SBATCH --nodes=2
#SBATCH --ntasks=32
##SBATCH --exclusive

# run for five minutes
#              d-hh:mm:ss
#SBATCH --time=0-24:00:00

# short partition should do it
#SBATCH --partition normal

# 500MB memory per core
# this is a hard limit 
#SBATCH --mem-per-cpu=1500MB


# define and create a unique scratch directory
SCRATCH_DIRECTORY=/global/work/${USER}/${SLURM_JOBID}
mkdir -p ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}

#module load VASP/5.3.5
#executable=vasp
module purge 
module load StdEnv
module load StdMod
module load intel/13.0
module load OpenMPI/1.8.5-intel-13.0
module load VASP/5.4.1.plain-intel-2016a
# we copy everything we need to the scratch directory
# ${SLURM_SUBMIT_DIR} points to the path where this script was submitted from
cp ${SLURM_SUBMIT_DIR}/{KPOINTS,INCAR,POTCAR,POSCAR} ${SCRATCH_DIRECTORY}

# we execute the job and time it
time mpirun -np 32 /global/hds/software/cpu/non-eb/VASP/5.4.1.plain-intel-2016a/bin/vasp_std > vasp_out

# after the job is done we copy our output back to $SLURM_SUBMIT_DIR
cp ${SCRATCH_DIRECTORY}/* ${SLURM_SUBMIT_DIR}

# we step out of the scratch directory and remove it
cd ${SLURM_SUBMIT_DIR}
rm -rf ${SCRATCH_DIRECTORY}

# happy end
exit 0
 
