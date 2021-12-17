#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Bio import SeqIO
from pyopenms import *
from urllib.request import urlretrieve
dig = ProteaseDigestion()
dig.getEnzymeName() # Trypsin
dighum = "".join([l.strip() for l in open("uniprothomo.fasta").readlines()[1:]])

for line in dighum:
    
    if line.startswith(">"): continue
    #dighum.write(line.strip()) 
    dighum = AASequence.fromString(dighum)
# create all digestion products
result = []
dig.digest(dighum, result)
for i in result :
    print(i.toString())
len(result)
 

