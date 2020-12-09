# -*- coding: utf-8 -*-
"""
Created Dec 2020
@author: Tom Nelson
"""

### Cleaning up the amino acid sequences ###

# import regular expression module
import re

# example sequence
aa1 = 'YSKLMS(AEEE)KLM-Abi-KK(AEEE)AEEE)LSKLSKM(AEEE)MS'

# remove side chain stuff (anything between parenthesis)
# Not sure if this generalizes to all cases, needs testing
# aa2 = re.sub(r'\(+[^(]*\)', '', aa1)

# replace 'Abi' with 'A'
# aa3 = aa2.replace('Abi', 'A')

# remove dashes
# aa4 = aa3.replace('-', '')

def cleanup(seq):
    mySeq = re.sub(r'\(+[^(]*\)', '', seq)
    mySeq = mySeq.replace('Abi', 'A')
    mySeq = mySeq.replace('-', '')
    return mySeq
cleanup(aa1)


