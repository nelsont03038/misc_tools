# Pairwise global alignment in BioPython
# takes any character for alignment
# good for computing pairwise distances (alignment score)
# can customize scoring with penalties and algorithm choice
from Bio.pairwise2 import align, format_alignment
alignments = align.globalxx("CYPCYYLAPCM", "ATPCYYYLAPCM") 

# iterate through alignments or get max score
[a for a in alignments]
max([a.score for a in alignments])

# look at pretty alignments
for a in alignments:
    print(format_alignment(*a))

print(alignments[2].score)
print(format_alignment(*alignments[2]))



# calling an external alignment tool (MUSCLE)
from Bio.Align.Applications import MuscleCommandline
cline = MuscleCommandline(input='example_protein.fasta', 
                          out='example_protein_aligned.fasta', 
                          verbose=True)
stdout, stderr = cline()
print(stderr)




# mafft can do arbitrary text alignment
from Bio.Align.Applications import MafftCommandline
mafft_exe = "mafft.bat"
in_file = "input.fa"
mafft_cmnd = MafftCommandline(mafft_exe, input=in_file)
print(mafft_cmnd)

stdout, stderr = mafft_cmnd()
with open("aligned.fa", "w") as handle:
    handle.write(stdout)
from Bio import AlignIO
align = AlignIO.read("aligned.fa", "fasta")


# calling without Biopython do add options to command
in_file = "input.fa"
out_file = "output.fa"
cmnd = "mafft --auto --quiet " + in_file + " > " + out_file
import subprocess
output = subprocess.getoutput(cmnd)
print(output)


# this works best #
## get subprocess module 
import subprocess
p = subprocess.Popen("mafft --auto --quiet --text --clustalout input.fa", stdout=subprocess.PIPE, shell=True)
output, err = p.communicate()
output = output.decode("utf-8") 
print(output)



# Levenshtein distance in python
import numpy as np

def levenshtein(seq1, seq2):
    size_x = len(seq1) + 1
    size_y = len(seq2) + 1
    matrix = np.zeros ((size_x, size_y))
    for x in range(size_x):
        matrix [x, 0] = x
    for y in range(size_y):
        matrix [0, y] = y

    for x in range(1, size_x):
        for y in range(1, size_y):
            if seq1[x-1] == seq2[y-1]:
                matrix [x,y] = min(
                    matrix[x-1, y] + 1,
                    matrix[x-1, y-1],
                    matrix[x, y-1] + 1
                )
            else:
                matrix [x,y] = min(
                    matrix[x-1,y] + 1,
                    matrix[x-1,y-1] + 1,
                    matrix[x,y-1] + 1
                )
    print (matrix)
    return (matrix[size_x - 1, size_y - 1])

levenshtein("CYPCYYLAPCM", "ATPCYYYLAPCM")




#### Parsing out the sequences
import re
import pandas as pd

data = pd.read_csv('real_data.csv')

aa2 = 'YSKLMS(AEEE)KLM-Abi-KK(AEEE)AEEE)LSKLSKM(AEEE)MS'
aa1 = 'YAibEGTFTSDYSIAibLDKIAQK((miniPEG)2 E1 C18 COOH)A(1Nal)VQWLIAGGPSSGAPPPS-NH2'
aa1 = 'YAibEGTFTSDYSIAibLDKIAQK((miniPEG)2 E1 C20 COOH)A-1Nal-VQWLIAGGPSSGAPPPS-NH2'

mySeq = aa1
mySeq = mySeq.replace('(1Nal)', 'A')
mySeq = mySeq.replace('-1Nal-', 'A')
mySeq = re.sub(r'\(+[^(]*\)', '', mySeq)
mySeq = mySeq.replace('Aib', 'A')
mySeq = mySeq.replace('-NH2', '')
mySeq = mySeq.replace('-OH', '')
mySeq = mySeq.replace('-', '')
mySeq = mySeq.replace(' ', '')


# remove dashes
aa1 = re.sub('[-]', ' ', aa1)

# remove between parenthesis
# fails with unbalanced and nested parens
re.sub(r'\([^)]*\)', ' ', aa1)

#!!! Works with unbalanced parens
# fails with balanced nested parens
re.sub(r'\([^(]*\)', ' ', aa1)


###############################################
#!!! Works with unbalanced parens
#!!! Works with balanced that start with 2
re.sub(r'\(+[^(]*\)', ' ', aa1)
##############################################
       

#!!! Works with balanced nested parents
#!!! Works with unbalanced nested parens
# NEED TO FIND ANOTHER CHARACTER THAT IS NEVER 
# WITHIN THE PARENS LIKE THE -
re.sub(r'\([^-]*\)', ' ', aa1)

#### FINAL VERSION THAT WORKS!!!
re.sub(r'\(+[^(]*\)', '', aa1)
###############################





################################
# Sequence Graph Tree research #
################################

import numpy as np
import pandas as pd
from sgt import SGT
sgt = SGT(flatten=True)
sequence = np.array(["B","B","A","C","A","C","A","A","B","A"])
sgt.fit(sequence)

# A sample corpus of two sequences.
corpus = pd.DataFrame([[1, ["B","B","A","C","A","C","A","A","B","A"]], 
                       [2, ["C", "Z", "Z", "Z", "D"]]], 
                      columns=['id', 'sequence'])
corpus

# Learning the sgt embeddings as vector for
# all sequences in a corpus.
# mode: 'default'
sgt = SGT(kappa=1, 
          flatten=True, 
          lengthsensitive=False, 
          mode='default')
sgt.fit_transform(corpus)



data = pd.read_csv('sgt_sequences.csv')
data['seq_list'] = data['sequence'].apply(lambda x: list(x))
data.drop(['sequence'], axis=1, inplace=True)
data.columns = ['id', 'sequence']

sgt = SGT(kappa=1, 
          flatten=True, 
          lengthsensitive=False, 
          mode='default')
embedding = sgt.fit_transform(data)

embedding = embedding.set_index('id')

query_protein = 'YAQGTFTSDYSILLDEKAQAAFIQYLLAAGPSSGAPPPS'

# Step 1. Compute sgt embedding for the query protein.
query_protein_sgt_embedding = sgt.fit(list(query_protein))

# Step 2. Compute the dot product of query embedding 
# with the protein embedding database.
similarity = embedding.dot(query_protein_sgt_embedding)

# Step 3. Return the top k protein names based on similarity.
similarity = similarity.sort_values(ascending=False)



