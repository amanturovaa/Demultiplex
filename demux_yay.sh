#!/bin/bash
#script to demultiplex

#SBATCH --account=bgmp                    #REQUIRED: which account to use
#SBATCH --partition=bgmp                  #REQUIRED: which partition to use
#SBATCH --cpus-per-task=8                 #optional: number of cpus, default is 1
#SBATCH --mem=16GB                        #optional: amount of memory, default is 4GB per cpu
#SBATCH --job-name=demux_yay           #optional: job name
#SBATCH --output=demux_yay%j.out       #optional: file to store stdout from job, %j adds the assigned jobID
#SBATCH --error=demux_yay%j.err        #optional: file to store stderr from job, %j adds the assigned jobID

/usr/bin/time -v \
python demux_yay.py \
-r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
-r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
-r3 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
-r4 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
-x /projects/bgmp/shared/2017_sequencing/indexes.txt \
-q 32 \
-o demux_yay