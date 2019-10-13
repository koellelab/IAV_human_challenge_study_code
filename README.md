# The Dynamics of Influenza A H3N2 Defective Viral Genomes from a Human Challenge Study

These are a set of simple scripts designed to identify defective viral genomes (DVGs) in short read illumina data. This code was used to process samples from an influenza human challenge study (Sobel Leonard et al. 2016, Sobel Leonard et al. 2017) and track the evolutionary trajectory of DVG populations within single hosts. Raw data for this analysis is available at <sra accession>. We were lucky to have two technical sequencing replicates (designated "wexo" and "noexo" as one run included exonuclease, all other processing steps were identical) for each sample so so some of the scripts deal with reconciling results between the two replicates.

There scripts include: 
1. pipeline_batch.sh
	This is the SLURM submission script to run the entire analysis (short of data visualization). Some of the directories are specific to where our data were located on the Bridges cluster at the Pittsburgh Supercomputer Center, you would need to change these directories to run this analysis on your data. The script also assumes you have two technical replicates for your data. It expects a FASTA reference sequence titled H3N2_CDS.fasta (easy to change in the code).

	The script performs the following steps:
	1. Removes chimeric reads (as identified by read IDs shared between) from each technical replicate using the remove_chimeras.py script
	2. Processes all samples in each technical replicate using the pipeline.sh script to identify DVGs
	3. Finds DVGs on a per-sample basis which are in both technical replicates using the replicates.sh script. 

2. remove_chimeras.py
	This script is takes a set of FASTQ files as input and identifies mispaired reads (reads whose mate is in another sample). These mispaired reads are then removed and a "clean" fastq file is written. This script was added to our analysis in response to the Xue & Bloom Nature Genetics 2019 correspondence as the data for our analysis were also sequenced at JCVI. However, the levels of mispaired reads in these data were considerably less (~6%) as compared to the Poon et al. dataset. This script requires the pysam module to run correctly. To run the script simply provide your FASTQ file names as arguments.

3. pipeline.sh
	This script does the majority of the analysis. It takes as inputs a reference sequence and a list of FASTQ files and will identify DVGs (identified by split reads) in the reads. Some of the directories in the script are specific to where our data was stored and will need to be adjusted for it to run properly in another environment. It performs a number of steps:
	1. Runs FASTQC on all input samples
	2. Uses Kraken2 with a viral database to identify only influenza A reads
	3. Trims low quality bases with trimmomatic
	4. Aligns reads to the reference sequence with STAR
	5. Calls consensus variants with samtools and creates sample-wise reference sequences
	6. Realigns reads to the sample-wise reference using STAR in basic two-pass mode
	7. Removes duplicate reads using PIcard
	8. Identifies junction spanning reads using the filter_reads.py script
	9. Tabulates DVG support (# of reads) and normalizes this by the number of mapped reads in that gene
	10. Calcualtes read depths

	This script requires a number of software pieces to run:
	* FASTQC
	* STAR
	* Kraken2 with a databse name "viral-db"
	* Trimmomatic
	* STAR
	* SAMtools/BCFtools
	* Picard
	* Pysam

	To run this script use the command: 
	```
	bash pipeline.sh -f REFERENCE_FASTA -o OUT_DIR INPUT_FILES
	```

4. filter_reads.py
	This is a Python 2 script which takes a BAM file as input and identifies split reads. It will output a seperate SAM file including ONLY split reads. It identfies split reads using Pysam and with the following criteria: 
	* "N" in the cigar string
	* A minimum of 5 consecutive alignment reference matches
	* A minimum of 15 total alignment reference matches
	* No more than 3 small indels
	* A minimum deletion size of 100

5. replicates.sh
	This script simply identifies DVGs which are present in both technical replicates and outputs a consolidated output file. It also produces reach depth and reach length files which are averaged between the two technical replicates. 



## Reference
>*In preperation*

## Contact
>michael.martin2@emory.edu