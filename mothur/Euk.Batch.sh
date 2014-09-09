#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,walltime=08:00:00
#PBS -M mmuscare@indiana.edu
#PBS -m abe
#PBS -j oe
cd /N/dc2/projects/Lennon_Sequences/P-meso_Euk/
mothur Euk.Batch
