#!/usr/bin/env python
# coding: utf-8

# In[21]:


from pyopenms import *
tsg = TheoreticalSpectrumGenerator()
spec = MSSpectrum()
pep = AASequence.fromString("MLAVKDIFA")
p = Param()
p.setValue("add_b_ions", "false")
p.setValue("add_metainfo", "true")
tsg.setParameters(p)
tsg.getSpectrum(spec , pep,1,1)
spec.size()
print("Spectrum 1 of", pep, "has", spec.size(), "peaks.")
for ion, peak in zip(spec.getStringDataArrays()[0], spec):
    print(ion.decode(), "is generated at m/z", peak.getMZ())
    


# In[22]:


import matplotlib.pyplot as plt
plt.bar(spec.get_peaks()[0], spec.get_peaks()[1], snap=False) 
plt.xlabel("m/z")
plt.ylabel("intensity")


# In[24]:


from pyopenms import *
tsg = TheoreticalSpectrumGenerator()
spec = MSSpectrum()
pep = AASequence.fromString("PEPTIDERDLQMTQSPSSLSVSVGDRPEPTIDE")
p = Param()
p.setValue("add_b_ions", "false")
p.setValue("add_metainfo", "true")
tsg.setParameters(p)
tsg.getSpectrum(spec , pep,1,1)
spec.size()
print("Spectrum 1 of", pep, "has", spec.size(), "peaks.")
for ion, peak in zip(spec.getStringDataArrays()[0], spec):
    print(ion.decode(), "is generated at m/z", peak.getMZ())


# In[26]:


import matplotlib.pyplot as plt
plt.bar(spec.get_peaks()[0], spec.get_peaks()[1], snap=False) 
plt.xlabel("m/z")
plt.ylabel("intensity")


# In[31]:


from pyopenms import *
from urllib.request import urlretrieve
gh = "https://raw.githubusercontent.com/OpenMS/pyopenms-extra/master"
urlretrieve (gh + "/src/data/P02769.fasta", "bsa.fasta")

dig = ProteaseDigestion()
dig.getEnzymeName() # Trypsin
bsa = "".join([l.strip() for l in open("bsa.fasta").readlines()[1:]])
bsa = AASequence.fromString(bsa)
# create all digestion products
result = []
dig.digest(bsa, result)
print(result[4].toString())
len(result) # 82 peptides
for ii in result :
    print(ii.toString())


# In[42]:


coun =0
for i in result :
    tsg = TheoreticalSpectrumGenerator()
    spec = MSSpectrum()
    p = Param()
    p.setValue("add_b_ions", "false")
    p.setValue("add_metainfo", "true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec , i,1,1)
    coun = coun+1
    print("Spectrum" ,coun, "of", i.toString(), "has", spec.size(), "peaks.")
    for ion, peak in zip(spec.getStringDataArrays()[0], spec):
        print(ion.decode(), "is generated at m/z", peak.getMZ())
    


# In[ ]:




