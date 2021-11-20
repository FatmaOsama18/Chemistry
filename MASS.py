#!/usr/bin/env python
# coding: utf-8

# In[88]:


from pyopenms import *
seq = AASequence.fromString("VAKA")

print("The peptide", str(seq), "consists of the following amino acids:")
for a in seq:
    print(a.getName(), ":", a.getMonoWeight())
listm = [] 
for r in seq:
    listm.append(r.getMonoWeight())  
  
print(seq.getMonoWeight())  
print(listm)
print(sum(listm))

sum(listm)== seq.getMonoWeight()


# In[ ]:





# In[ ]:




