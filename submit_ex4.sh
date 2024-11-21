#!/bin/bash
#SBATCH --job-name=exaust
#SBATCH --output=output_ex4.txt
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1Gb
#SBATCH --time=00:15:00
#SBATCH --partition=espec

mpirun -np 4 ./Ex4