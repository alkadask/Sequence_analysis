# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 14:22:24 2020

@author: AD
"""

#script to find tandem repeats and inverted repeats in DNA sequences
import string
import numpy as np
import matplotlib.pyplot as plt

#read a sequence from a file
file1 = open("C:/Users/alaka/Google Drive/test.ape","r")
content = file1.read()

#extract the sequence only from the file, removing spaces, numbers and next lines
rawseq = content[content.index("ORIGIN")+6:]
bases = list([val for val in rawseq if val.isalpha()]) 
seq = "".join(bases) 

#get reverse complement of seq
revcomp = seq.translate(string.maketrans("atgc", "tacg"))[::-1]

#finding tandem repeats
tr=np.zeros((len(seq),len(seq)))
for x in np.arange(len(seq)):
    for y in np.arange(x, len(seq)):
        if seq[x:x+1]==seq[y:y+1]: tr[x,y]=tr[x,y]+1
        if seq[x:x+2]==seq[y:y+2]: tr[x,y]=tr[x,y]+1
        if seq[x:x+3]==seq[y:y+3]: tr[x,y]=tr[x,y]+1
        if seq[x:x+4]==seq[y:y+4]: tr[x,y]=tr[x,y]+1
        if seq[x:x+5]==seq[y:y+5]: tr[x,y]=tr[x,y]+1
        if seq[x:x+6]==seq[y:y+6]: tr[x,y]=tr[x,y]+1
        if seq[x:x+7]==seq[y:y+7]: tr[x,y]=tr[x,y]+1
        if seq[x:x+8]==seq[y:y+8]: tr[x,y]=tr[x,y]+1
        if seq[x:x+9]==seq[y:y+9]: tr[x,y]=tr[x,y]+1
        if seq[x:x+10]==seq[y:y+10]: tr[x,y]=tr[x,y]+1
plt.figure(1, figsize=(20,8))
plt.subplot(121)
plt.title('Tandem repeats')
plt.imshow(tr)

#finding inverted repeats
ir=np.zeros((len(seq),len(seq)))
for x in np.arange(len(seq)):
    for y in np.arange(x, len(revcomp)):
        if seq[x:x+1]==revcomp[y:y+1]: ir[x,y]=ir[x,y]+1
        if seq[x:x+2]==revcomp[y:y+2]: ir[x,y]=ir[x,y]+1
        if seq[x:x+3]==revcomp[y:y+3]: ir[x,y]=ir[x,y]+1
        if seq[x:x+4]==revcomp[y:y+4]: ir[x,y]=ir[x,y]+1
        if seq[x:x+5]==revcomp[y:y+5]: ir[x,y]=ir[x,y]+1
        if seq[x:x+6]==revcomp[y:y+6]: ir[x,y]=ir[x,y]+1
        if seq[x:x+7]==revcomp[y:y+7]: ir[x,y]=ir[x,y]+1
        if seq[x:x+8]==revcomp[y:y+8]: ir[x,y]=ir[x,y]+1
        if seq[x:x+9]==revcomp[y:y+9]: ir[x,y]=ir[x,y]+1
        if seq[x:x+10]==revcomp[y:y+10]: ir[x,y]=ir[x,y]+1
plt.subplot(122)
plt.title('Inverted repeats')
plt.imshow(ir)
plt.show()