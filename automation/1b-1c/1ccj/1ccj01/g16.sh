#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --partition=gpu-a100          
#SBATCH --nodes=1
#SBATCH -A beb00042                     
#SBATCH --mem=128G                     
#SBATCH --ntasks=4             
#SBATCH --gres=gpu:1             

module load cuda/11.8
module load gaussian/16.C02

g16 1ccj01.com
