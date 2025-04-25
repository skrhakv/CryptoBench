#!/bin/bash
#SBATCH --partition=mpi-homo-long           # partition you want to run job in
#SBATCH --time=167:00:00           # walltime for the job in format (days-)hours:minutes:seconds
#SBATCH --mem=100000               # memory resource
#SBATCH --job-name="extract-ahojdb"     # change to your job name
#SBATCH --output=/home/skrhakv/CryptoBench/data/A-filter-ahojdb/holo_only_data/output.txt       # stdout and stderr output file
#SBATCH --mail-user=v.skrhak@gmail.com # send email when job changes state to email address user@example.com


DATADIR=/home/skrhakv/CryptoBench/src/A-filter-ahojdb

# clean the output file
rm /home/skrhakv/CryptoBench/data/A-filter-ahojdb/holo_only_data/pairs.csv

# activate venv
source /home/skrhakv/CryptoBench/holo_venv/bin/activate

cd $DATADIR
python3 main.py
