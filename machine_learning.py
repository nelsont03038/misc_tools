# -*- coding: utf-8 -*-
"""
Created in Dec 2020

@author: Tom Nelson
"""
# required modules
import pandas as pd
from modlamp.sequences import Random, Helices, Oblique

# Generate synthetic peptide data of 3 classes
# one class is random, other are predicted helix structure
ran = Random(3000, lenmin=10, lenmax=20)
ran.generate_sequences(proba='randnoCM')
ran_df = pd.DataFrame(ran.sequences, columns=['sequence'])
ran_df['label'] = "random"

hel = Helices(3000, lenmin=10, lenmax=20)
hel.generate_sequences()
hel_df = pd.DataFrame(hel.sequences, columns=['sequence'])
hel_df['label'] = "helix"

obl = Oblique(3000, lenmin=10, lenmax=20)
obl.generate_sequences()
obl_df = pd.DataFrame(obl.sequences, columns=['sequence'])
obl_df['label'] = "oblique"

data = pd.concat([ran_df, hel_df, obl_df], axis = 0)
data = data.sample(frac=1).reset_index(drop=True)
data['sequence'] = data['sequence'].apply(lambda x: list(x))

del hel, hel_df, ran, ran_df, obl, obl_df

y = data.label
data.drop(['label'], axis=1, inplace=True)
data.insert(0, 'id', list(data.index))

from sgt import SGT
sgt = SGT(kappa=1, 
          flatten=True, 
          lengthsensitive=True, 
          mode='default')
embedding = sgt.fit_transform(data)
embedding.drop(['id'], axis=1, inplace=True)


# Splitting the dataset into the Training set and Test set
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(embedding, 
                                                    y, 
                                                    test_size = 0.30, 
                                                    random_state = 42)

# Fitting Random Forest Classifier to the Training set
from sklearn.ensemble import RandomForestClassifier
classifier = RandomForestClassifier(n_estimators = 100, 
                                    criterion = 'gini', 
                                    random_state = 42)
classifier.fit(X_train, y_train)

# Predicting the Test set results
y_pred = classifier.predict(X_test)
print(classifier.score(X_test, y_test))

# Making the Confusion Matrix
import matplotlib.pyplot as plt
from sklearn.metrics import plot_confusion_matrix
matrix = plot_confusion_matrix(classifier, X_test, y_test, cmap=plt.cm.Blues)
plt.title('Confusion Matrix for the Random Forest Classifier')
plt.show(matrix)
plt.show()


