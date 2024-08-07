#!/bin/bash
#SBATCH -c 1 # Number of cores requested
#SBATCH -t 100 # Runtime in minutes
#SBATCH --mem=100000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o active_learning-%j.out # Standard out goes to this file
#SBATCH -e active_learning-%j.err # Standard err goes to this filehostname
#SBATCH --mail-user=joeri.lenaerts36@gmail.com

# load modules
module load python/3.10.13-fasrc01
module load matlab/R2022b-fasrc01

# load python packages
mamba activate optimization

# run code
export PATH=$PATH:/n/sw/helmod/apps/centos7/Core/matlab/R2022b-fasrc01/bin/matlab
python3 'active_learning.py' 

mamba deactivate
