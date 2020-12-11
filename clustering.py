# -*- coding: utf-8 -*-
"""
Created Dec 2020
@author: Tom Nelson
"""
# Graph based clutering, sequence based clustering and 
# network based clustering.

# A pair wise similarity matrix must be created and used as an
# input for graph-based clustering.
# There are various ways to do that.

# example 1 - using Biopython tools
# Pairwise global alignment, and use the alignment 
# score as distance measure 
# can customize scoring with penalties and 
# algorithm choice

# example peptide sequences
seq1 = "YAEGTFTSDYSIALDKIAQKANALVQWLIAGGPSSGAPP"
seq2 = "HHAQGTFTSDYSKYLDEKKAKEFVEWLLEGGPSSG"

# perform pairwise global alignment
from Bio.pairwise2 import align, format_alignment
alignments = align.globalxx(seq1, seq2) 

# look at alignments
for a in alignments:
    print(format_alignment(*a))

# look at details of a specific alignment
print(alignments[2])
print(format_alignment(*alignments[2]))

# get max score (score of the best alignment)
# in this example, the score is the number of identities
score = max([a.score for a in alignments])

# get length of smaller peptide
min_length = min([len(seq1), len(seq2)])

# compute distance measure (in this case fractional identity)
distance = score / min_length
print(distance)



# creation of the distance matrix
from Bio import SeqIO
from Bio.pairwise2 import align, format_alignment
distance_matrix = []
for seq_record in SeqIO.parse("cleaned_sequences.fa", "fasta"):
    outside_sequence_id = seq_record.id
    outside_sequence = seq_record.seq._data
    outside_length = len(outside_sequence)
    for seq_record in SeqIO.parse("cleaned_sequences.fa", "fasta"):
        inside_sequence_id = seq_record.id
        inside_sequence = seq_record.seq._data
        inside_length = len(outside_sequence)
        alignments = align.globalxx(outside_sequence, inside_sequence) 
        score = max([a.score for a in alignments])
        min_length = min([outside_length, inside_length])
        distance = score / min_length
        distance_matrix.append([outside_sequence_id, inside_sequence_id, distance])
        
import pandas as pd
distance_matrix_df = pd.DataFrame(distance_matrix, columns =['LSN1', 'LSN2', 'distance'])
distance_matrix_df.to_csv('distance_matrix.csv', index=False)



# could use something like Levenshtein distance as well

# Levenshtein distance in python
# sort of like Hamming distance but handles
# sequences of different lengths
# could potentially use raw, uncleaned sequences
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
    #print (matrix)
    return (matrix[size_x - 1, size_y - 1])

levenshtein(seq1, seq2)



########################################### 
# Finding most similar peptides using SGT #
# Sequence Graph Tree embedding features  #
###########################################

import numpy as np
import pandas as pd
from sgt import SGT

# load some cleaned amino acid sequence data
data = pd.read_csv('sgt_sequences.csv')

# corpus needs to be lists of amino acids
data['seq_list'] = data['sequence'].apply(lambda x: list(x))
data.drop(['sequence'], axis=1, inplace=True)
data.columns = ['id', 'sequence']

# instantiate an sgt object
sgt = SGT(kappa=1, 
          flatten=True, 
          lengthsensitive=False, 
          mode='default')

# learn feature embeddings from sequences
embedding = sgt.fit_transform(data)
embedding = embedding.set_index('id')

# test with an example peptide query
query_protein = 'YAQGTFTSDYSILLDEKAQAAFIQYLLAAGPSSGAPPPS'

# Compute sgt embedding for the query protein.
query_protein_sgt_embedding = sgt.fit(list(query_protein))

# Compute the dot product of query embedding 
# with the protein embedding database.
similarity = embedding.dot(query_protein_sgt_embedding)

# sort names based on similarity.
similarity = similarity.sort_values(ascending=False)




########################################### 
# Classification of peptides using SGT    #
# Sequence Graph Tree embedding features  #
###########################################

# Loading data
data = pd.read_csv('protein_classification.csv')

# Data preprocessing
from sklearn.preprocessing import LabelEncoder
y = data['class']
encoder = LabelEncoder()
encoder.fit(y)
encoded_y = encoder.transform(y)

# create 'corpus' of sequences (to use NLP lingo)
corpus = data.loc[:,['id','sequence']]
corpus['sequence'] = corpus['sequence'].map(list)

# Sequence embedding
sgt_ = SGT(kappa=1, 
           lengthsensitive=False, 
           mode='default')
sgtembedding_df = sgt_.fit_transform(corpus)
X = sgtembedding_df.set_index('id')

# instantiate a support vector machine classifier
from sklearn import svm
svc = svm.SVC()

# train and evaluate SVM using k-fold cross validation
from sklearn.model_selection import cross_val_score
accuracies = cross_val_score(estimator = svc, 
                             X = X, 
                             y = encoded_y, 
                             cv = 5)
print("The mean accuracy is: ", accuracies.mean())
print("The standard deviation is: ", accuracies.std())

