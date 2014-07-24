#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH -n 1
#SBATCH --qos=ccmb-condo
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J fastqc_all
#SBATCH -o /users/bsiranos/projects/prostate_hic/script_outputs/2014-07-05/fastqc_all.out
#SBATCH -e /users/bsiranos/projects/prostate_hic/script_outputs/2014-07-05/fastqc_all.err

module load fastqc
cd /users/bsiranos/data/bsiranos/prostate_hic/sra
for i in $(ls *.fastq); do echo processing $i;fastqc -o /users/bsiranos/data/bsiranos/prostate_hic/sra/fastqc  $i;done

