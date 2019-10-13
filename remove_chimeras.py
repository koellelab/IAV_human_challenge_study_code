#!/usr/bin/env python3
'''
 This script takes a set of fastq files as input and removes
 any reads in which the /1 and /2 reads are assigned
 to different subjects (under the asusmption that these
 are chimeric reads). Ignores index in the identification 
 of chimeric reads. 

 Author: Michael A. Martin (michael.martin2@emory.edu)

'''

# Import required modules
import csv
import pysam
import sys
import os
from collections import Counter

# Creates a list of sample names

# Converts input fastq files to tab delimited files
def fastq_tab(sample):
	os.system("grep '^@SOLEXA' -A 3 "+sample+".fastq |grep -v -- '^--$'| paste - - - - | sort -k1 > " + sample+'.tab')

# Identifies unpaired reads in a tab delimited read files
def unpaired_reads(sample,unpaired_list):
	# Because the tab files are sorted we can check to see if the current read ID matches 
	# the previous read ID to identify reads which are/are not paired
	prev_read=''
	paired_index=[]
	index=0
	with open(sample+'.tab','rb') as infile:
		reader=csv.reader(infile,delimiter='\t')
		reads=list(reader)
	for line in reads:
		read_id=line[0].split('#')[0][26:]
		if read_id==prev_read:
			paired_index.append(index-1)
			paired_index.append(index)
		index+=1
		prev_read=read_id
	unpaired_index=list(set(range(index))-set(paired_index))

	#check that the parsing here on the read IDs is correct 
	unpaired_list.extend([sample,row] for row in set([reads[item][0].split('#')[0][26:] for item in unpaired_index]))

# Identifies chimeric reads and writes a new file with just the non-chimeric reads

def chimeric_reads(sample,unpaired_list,mispaired_reads):
	with open(sample+'.tab','rb') as infile:
		reader=csv.reader(infile,delimiter='\t')
		reads=list(reader)
	these_unpaired_reads=[item[1] for item in unpaired_list if item[0]==sample]
	these_clean_reads=[item for item in reads if item[0].split('#')[0][26:] not in mispaired_reads]
	with open(sample+'_clean.fastq','w') as outfile:
		for row in these_clean_reads:
			for item in row:
				outfile.write(item+'\n')
# Begining of script
samples=[]
for i in xrange(1,len(sys.argv)):
	if 'clean' not in sys.argv[i]:
		samples.append(sys.argv[i].split('/')[-1].replace('.fastq',''))

for sample in samples:
	print sample+' tabbed'
	fastq_tab(sample)
print('Fastq files converted to .tab')

unpaired_list=[]
for sample in samples:
	print sample+' unpaired'
	unpaired_reads(sample,unpaired_list)

with open('unpaired_read_names.txt','w') as outfile:
	writer=csv.writer(outfile,delimiter='\t')
	for row in unpaired_list:
		writer.writerow(row)

print('unpaired list written')
with open('unpaired_read_names.txt','rb') as infile:
	reader=csv.reader(infile,delimiter='\t')
	unpaired_list=list(reader)
print('unpaired list read')

unpaired_reads_flat=[item[1] for item in unpaired_list]
print('unpaired list flattened')
unpaired_reads_flat=Counter(unpaired_reads_flat)
print('unpaired list counted')

mispaired_reads=set([name for name,count in unpaired_reads_flat.items() if count>1])
print('Unpaired and mispaired reads identified')
for sample in samples:
	print sample+' chimeras removed'
	chimeric_reads(sample,unpaired_list,mispaired_reads)
print('Chimeric reads removed')
# Cleans up
for sample in samples:
	os.system('rm '+sample+'.tab')
print('Done')

