#!/usr/bin/env python2.7
import pysam

infile = pysam.AlignmentFile("-","rb")
outfile = pysam.AlignmentFile("-","w",template=infile) 
for read in infile:
		min_matches=5
		min_deletions=100
		all_matches=0
		N=0
		# Key for cigar parsing:
		# M=0
		# I=1
		# D=2
		# N=3
		# S=4
		# H=5
		# P=6
		# ==7
		# X=8
		# B=9
		# if the number of matched bases is less than the
		# current value of matches, set the value of matches to the number of matched bases
		for cigar in read.cigartuples:
			if cigar[0]==0 and cigar[1]<min_matches:
				min_matches=cigar[1]
			# Tabulates the total number of matched reference bases
			if cigar[0]==0:
				all_matches+=cigar[1]
			# If there is an N in the cigar string and if
			# at least 100 bases are deleted
			if cigar[0]==3 and cigar[1]<min_deletions:
				min_deletions=cigar[1]
		# if the minimum number of matched bases is >=15 (>=5 for each individual matched segment) 
		# and the minimum number of consecutive deleted reference bases is 100 
		# and there are no more than three indels in this read
		# and this is the primary alignment 
		# print read to stdout (outfile)
		if min_matches>=5 and min_deletions>=100 and \
		(read.cigarstring.count('I')+read.cigarstring.count('D'))<=3 and \
		'N' in read.cigarstring and all_matches>=15 and read.is_secondary==False:
			#puts all the reads with only one junction in one bam file and reads with >one junction in another bam file
			outfile.write(read)
