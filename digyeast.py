#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from Bio import SeqIO
from pyopenms import *
from urllib.request import urlretrieve
dig = ProteaseDigestion()
dig.getEnzymeName() # Trypsin
digyt = "".join([l.strip() for l in open("uniprot-yeast.fasta").readlines()[1:]])

for line in digyt:
    
    if line.startswith(">"): continue
    #dighum.write(line.strip()) 
    digyt = AASequence.fromString(digyt)
# create all digestion products
result = []
dig.digest(digyt, result)
for i in result :
    print(i.toString())
len(result)
 

