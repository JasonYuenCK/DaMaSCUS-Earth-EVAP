#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH -o job-%j
#SBATCH --mail-type=end
#SBATCH --mail-user=lingyuxia@link.cuhk.edu.hk

if command -v module >/dev/null 2>&1; then
	module load intel/2022.2
	module load cmake/3.21.2 boost/1.77.0
fi

export LD_LIBRARY_PATH=/project/kennyng/backup_DM/lib:$LD_LIBRARY_PATH
mpirun -n 16 ./DaMaSCUS-SUN "${1:-config_Lingyu.cfg}"
