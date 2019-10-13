#!/bin/bash

#SBATCH -N 1
#SBATCH -p RM
#SBATCH -t 08:00:00


# Removes chimeric reads
cd /pylon5/ms5fpbp/mmartin4/h3n2_human_challenge_dvg/wexo/clean_fastq
python /pylon5/ms5fpbp/mmartin4/h3n2_human_challenge_dvg/remove_chimeras.py ../raw_fastq/*.fastq 

cd /pylon5/ms5fpbp/mmartin4/h3n2_human_challenge_dvg/noexo/clean_fastq
python /pylon5/ms5fpbp/mmartin4/h3n2_human_challenge_dvg/remove_chimeras.py ../raw_fastq/*.fastq 

# Processes the wexo samples
bash pipeline.sh -f reference/h3n2/H3N2_CDS.fasta -o wexo wexo/clean_fastq/*.fastq 

# Processes the noexo samples
bash pipeline.sh -f reference/h3n2/H3N2_CDS.fasta -o noexo noexo/clean_fastq/*.fastq

# Identifies only DVGs present in both sequencing runs
bash replicates.sh -1 wexo -2 noexo 



