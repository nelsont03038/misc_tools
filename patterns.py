import re
import gzip
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
file_name = 'test.fastq.gz'
count = 0
with gzip.open(file_name,'rt') as file_handle:
    for name, sequence, qual_score in FastqGeneralIterator(file_handle):
        count += 1
        insert = get_insert(sequence, primer1 = "AGCT", primer2 = "TCGT")
        if insert:
            print(insert)

print("%i records in %s" % (count, file_name))



