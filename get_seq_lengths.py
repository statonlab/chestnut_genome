#!/lustre/projects/staton/software/Anaconda2-2.5.0/bin/python
import re, sys, getopt
from Bio import SeqIO

fasta_file = str(sys.argv[1])
len_file = fasta_file+'.lens.tsv'

## get lengths and print if over 10k
inhandle = open(fasta_file, 'rU')
outhandle = open(len_file, 'w')

for record in SeqIO.parse(inhandle, "fasta") :
	outhandle.write(record.id + "\t" + str(len(record.seq)) + "\n")

inhandle.close()
outhandle.close()
