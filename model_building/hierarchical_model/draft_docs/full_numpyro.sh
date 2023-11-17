#!/bin/bash
#
#SBATCH --mem=40g
#SBATCH --partition=gpu
#SBATCH --gres=gpu:p100:2

#SBATCH --cpus-per-task=14
#SBATCH --time=168:00:00
#SBATCH --mail-type=END
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err

python Tox21/PyMC_model_script_numpyro.py
