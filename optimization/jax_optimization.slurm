#!/bin/bash
#SBATCH -c 1 # Number of cores requested
#SBATCH -t 100 # Runtime in minutes
#SBATCH --gres=gpu:nvidia_a100-sxm4-80gb:1
#SBATCH -p gpu_requeue
#SBATCH --mem=100000 # Memory per node in MB (see also --mem-per-cpu)
#SBATCH --open-mode=append # Append when writing files
#SBATCH -o jax_optimization-%j.out # Standard out goes to this file
#SBATCH -e jax_optimization-%j.err # Standard err goes to this filehostname
#SBATCH --mail-user=joeri.lenaerts36@gmail.com

# load modules
module load python/3.10.13-fasrc01

# load python packages
mamba activate optimization

# memory overload
export TF_GPU_ALLOCATOR=cuda_malloc_async

#DCGM job statistics 
#group=$(dcgmi group -c joerigpus --default)
#if [ $? -eq 0 ]; then
#  groupid=$(echo $group | awk '{print $10}')
#  dcgmi stats -g $groupid -e
#  dcgmi stats -g $groupid -s $SLURM_JOBID
#fi

# run code
python3 'jax_optimization.py' phase_matching worst_case --epochs 1500 --verbose --lr 0.0001

#DCGM job statistics 
#OUTPUTDIR=$(scontrol show job $SLURM_JOBID | grep WorkDir | cut -d = -f 2)
#dcgmi stats -x $SLURM_JOBID
#dcgmi stats -v -j $SLURM_JOBID 
  
mamba deactivate
