# -*- coding: utf-8 -*-
"""
Created Dec 2020
@author: Tom Nelson
"""

from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO

# call external mafft for multiple sequence alignment 
mafft_exe = 'C:/Users/l006251/OneDrive - Eli Lilly and Company/Desktop/mafft-win/mafft.bat'
in_file = "cleaned_sequences.fa"
mafft_cmnd = MafftCommandline(mafft_exe, input=in_file)
print(mafft_cmnd)
stdout, stderr = mafft_cmnd()

# create a "handle" from the alignment output
my_handle = StringIO(stdout)

# create alignment object
align = AlignIO.read(my_handle, "fasta")
print(align)

# iterate through records in alignment object
for record in align:
    print("%s - %s" % (record.seq, record.id))



