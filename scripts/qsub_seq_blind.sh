#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# This section specifies uger requests.  
# This is good for jobs you need to run multiple times so you don't forget what it needs.

# Memory request for 30G
#$ -l h_vmem=40G

# Cores
#$ -j y
#$ -o /broad/thechenlab/Chenlei/spotmapping/fiducial/qsub_out

# Runtime request.  
#$ -l h_rt=10:00:00


######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts

# Use your dotkit
source /broad/software/scripts/useuse

reuse -q Anaconda3

source activate /home/unix/chu/anaconda3/envs/slidelock

##################
### Run script ###
##################

cd /broad/thechenlab/Chenlei/spotmapping/fiducial

python fiducial_seq_blind_whitelist.py  -d 240311 -s H11_1_rec
python fiducial_seq_blind_whitelist.py  -d 240311 -s H11_4_rec -r2 V15
