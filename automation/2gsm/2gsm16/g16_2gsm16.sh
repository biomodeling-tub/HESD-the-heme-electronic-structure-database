#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mail-type=begin # send email when job begins
#SBATCH --mail-type=end # send email when job ends
#SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
#SBATCH --partition=gpu-a100
#SBATCH --nodes=1
#SBATCH -A beb00042
#SBATCH --mem=32G
#SBATCH --ntasks=72
#SBATCH --gres=gpu:4

module load cuda/11.8
module load gaussian/16.C02
module load anaconda3/2023.09



g16 2gsm16.com
python ../../beb00042_agent.py "2gsm" "16"
cd ../
cp -r 2gsm16/ /../../scratch/projects/beb00042/2gsm/
