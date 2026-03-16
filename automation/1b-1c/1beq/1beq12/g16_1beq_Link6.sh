#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mail-type=begin # send email when job begins
#SBATCH --mail-type=end # send email when job ends
#SBATCH --mail-user=elizaveta.zhartovska@campus.tu-berlin.de
#SBATCH --partition=gpu-a100          
#SBATCH --nodes=1                    
#SBATCH --mem=16G                     
#SBATCH --ntasks=72             
#SBATCH --gres=gpu:4             

module load cuda/11.8
module load gaussian/16.C02

g16 1beq_Link6.com
