#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mail-type=end # send email when job ends
#SBATCH --mail-user=j.jones@tu-berlin.de
#SBATCH --partition=gpu-a100
#SBATCH --nodes=1
#SBATCH -A beb00042
#SBATCH --mem=32G
#SBATCH --ntasks=72
#SBATCH --gres=gpu:4

module load cuda/11.8
module load gaussian/16.C02
module load anaconda3/2023.09

g16 1cxy16re.com
