# -*- coding: utf-8 -*-
"""
Created Dec 2020
@author: Tom Nelson
"""

# mafft 
from Bio.Align.Applications import MafftCommandline
mafft_exe = "mafft.bat"
in_file = "input.fa"
mafft_cmnd = MafftCommandline(mafft_exe, input=in_file)
print(mafft_cmnd)
stdout, stderr = mafft_cmnd()

print(stdout)


