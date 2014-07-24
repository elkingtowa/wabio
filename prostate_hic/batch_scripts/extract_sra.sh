#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH -n 1
#SBATCH --qos=ccmb-condo
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -J dump_fastq
#SBATCH -o /users/bsiranos/projects/prostate_hic/script_outputs/2014-07-05/dump_fastq.out
#SBATCH -e /users/bsiranos/projects/prostate_hic/script_outputs/2014-07-05/dump_fastq.err

module load sratoolkit/2.2.0
cd /users/bsiranos/data/bsiranos/prostate_hic/sra
for i in $(ls *.sra); do echo dumping $i;fastq-dump --split-3 $i;done

