# -*- coding: utf-8 -*-
"""
Created Dec 2020
@author: Tom Nelson
"""

### Cleaning up the amino acid sequences ###

# import needed modules
import re
import pandas as pd

# load some data
data = pd.read_csv('real_data_all.csv')

# function to cleanup aa sequences for alignment
# subject to change and improvement
def cleanup(seq):
    mySeq = seq
    mySeq = mySeq.replace('(1Nal)', 'A')
    mySeq = mySeq.replace('-1Nal-', 'A')
    mySeq = mySeq.replace('Nal', 'A')
    mySeq = re.sub(r'\(+[^(]*\)', '', mySeq)
    mySeq = mySeq.replace('Aib', 'A')
    mySeq = mySeq.replace('-NH2', '')
    mySeq = mySeq.replace('-OH', '')
    mySeq = mySeq.replace('-', '')
    mySeq = mySeq.replace(' ', '')
    return mySeq

data['clean'] = data['sequence'].apply(cleanup)
data.to_csv('cleaned_all.csv', index=False)

# save cleaned sequences to fasta file for alignment
myFile = open('cleaned_sequences.fa', 'w')
for index, row in data.iterrows():
    myString = ">" + str(row['LSN']) + "\n" + str(row['clean'] + "\n")
    myFile.write(myString)
myFile.close()


