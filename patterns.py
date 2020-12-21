# -*- coding: utf-8 -*-
"""
Created in Dec 2020

@author: Tom Nelson
"""
# required modules
import re
import gzip
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# function to extract the sequence in between 2 specified 
# primer or adapter sequences
def get_insert(sequence, primer1, primer2):
    pattern = primer1 + "(.+?)" + primer2
    match = re.search(pattern, sequence)
    if match:
        return match.group(1)
    else:
        return 0

# test it
sequence = 'TGAGCTCGAGACTCGGGGTGACAGCTCTTCATACATAGAGCGGGGGCGTCGAACGGTCGTGAAAGTCATAGTACCCCGGGTACCAACTTACTGAGGATAT'
get_insert(sequence, primer1 = "AGCTCTTCATACATAG", primer2 = "TCGTGAAAGTCATAGT")



# open and parse through a gzipped fastq file
file_name = 'peptideTest.fastq.gz'
count_reads = 0
count_inserts = 0
leftPrimer = 'GCTGCCCAACCAGCC'
rightPrimer = 'TAG'
with gzip.open(file_name,'rt') as file_handle:
    for name, sequence, qual_score in FastqGeneralIterator(file_handle):
        count_reads += 1
        insert = get_insert(sequence, primer1 = leftPrimer, primer2 = rightPrimer)
        if insert:
            count_inserts += 1
            print(Seq(insert).translate(to_stop=True))

print("%i reads and %i inserts in %s" % (count_reads, count_inserts, file_name))



