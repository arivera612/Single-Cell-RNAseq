#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=1-00:15:00     # 1 day and 15 minutes
#SBATCH --output=demultiplex.stdout
#SBATCH --mail-user=arive019@ucr.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name="demultiplex_ddseq"
#SBATCH -p highmem # This is the default partition, you can use any of the following; intel, batch, highmem, gpu


# Print current date
date

# Load samtools
module load python/2.7.13

# Change directory to where you submitted the job from, so that relative paths resolve properly
SLURM_SUBMIT_DIR=/bigdata/messaoudilab/arivera/SVV/ddSEq/data
cd $SLURM_SUBMIT_DIR

# Concatenate BAMs
python demultiplexDDseq.py d7_33/33d7_S4_R1_001.fastq d7_33/33d7_S4_R2_001.fastq d733- Cell Barcodes Samples

# Print name of node
hostname
