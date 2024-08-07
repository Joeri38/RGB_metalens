#!/bin/bash
#SBATCH -c 1 # Number of cores requested
#SBATCH -t 100 # Runtime in minutes
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH -p gpu_requeue
#SBATCH --mem=100000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o jax_optimization_9_pillars-%j.out # Standard out goes to this file
#SBATCH -e jax_optimization_9_pillars-%j.err # Standard err goes to this filehostname
#SBATCH --mail-user=joeri.lenaerts36@gmail.com

# load modules
module load python/3.10.13-fasrc01

# load python packages
mamba activate optimization
pip install torch-jax-interop

# memory overload
export TF_GPU_ALLOCATOR=cuda_malloc_async

# run code
python3 'jax_optimization_9_pillars.py' phase_matching worst_case --epochs 1500 --verbose --lr 0.0001
  
mamba deactivate
