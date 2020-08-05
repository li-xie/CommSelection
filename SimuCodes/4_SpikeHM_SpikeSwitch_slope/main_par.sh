#!/bin/bash
# matlab batch submission script

#SBATCH --partition=largenode     # more likely to get 12 cores in the full queue
#SBATCH --time=3-00:00:00        # 6 hours
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12       # request 16 cores
# #SBATCH --mail-user=youremail@fhcrc.org
#SBATCH --mail-type=END
#SBATCH --job-name=Matlab
#SBATCH --output=matlab-%j.out
#SBATCH --error=matlab-%j.err
#SBATCH --mem=217800
#SBATCH --qos=matlab

# grab environment variables for debugging
#env | grep SLURM

# set the path for the current version of matlab
ml MATLAB/R2019a

matlab -r "main_HM_SPNoise_HeriSwitch;exit" -nodisplay -nosplash -nodesktop