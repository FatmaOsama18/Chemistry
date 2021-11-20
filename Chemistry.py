#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pyopenms
help(pyopenms.Constants)


# In[6]:


print("Planck's num is",pyopenms.Constants.PLANCK)
print("speed of light's num is",pyopenms.Constants.SPEED_OF_LIGHT)


# In[28]:


from pyopenms import *
edb = ElementDB()
edb.hasElement("K")
edb.hasElement("V")

nitrogen = edb.getElement("N")
print(nitrogen.getName())
print(nitrogen.getSymbol())
print(nitrogen.getMonoWeight())
print(nitrogen.getAverageWeight())

hydrogen = edb.getElement("H")
print(hydrogen.getName())
print(hydrogen.getSymbol())
print(hydrogen.getMonoWeight())
print(hydrogen.getAverageWeight())
isotopes = nitrogen.getIsotopeDistribution()
print(isotopes)
print ("One mole of potassium weighs", 2*nitrogen.getAverageWeight(), "grams")
print ("One mole of 16O2 weighs", 2*nitrogen.getMonoWeight(), "grams")


# In[30]:


nitrogen_isoDist = {"mass": [], "abundance": []}
hydrogen_isoDist = {"mass": [], "abundance": []}

for iso in isotopes.getContainer():
    print ("Nitrogen isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    nitrogen_isoDist["mass"].append(iso.getMZ())
    nitrogen_isoDist["abundance"].append((iso.getIntensity() * 100))

hydrogen = edb.getElement("S")
isotopes = hydrogen.getIsotopeDistribution()
for iso in isotopes.getContainer():
    print ("Hydrogen isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    hydrogen_isoDist["mass"].append(iso.getMZ())
    hydrogen_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[34]:


import math
from matplotlib import pyplot as plt

# very simple overlappping correction of annotations
def adjustText(x1, y1, x2, y2):
    if y1 > y2:
        plt.annotate('%0.3f' % (y2), xy=(x2, y2), xytext=(x2+0.5,y2+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')
    else:
        plt.annotate('%0.3f' % (y1), xy=(x1, y1), xytext=(x1+0.5,y1+9),
                     textcoords='data',
                     arrowprops=dict(arrowstyle="->", color='r', lw=0.5),
                     horizontalalignment='right', verticalalignment='top')


def plotDistribution(distribution):
    n = len(distribution["mass"])
    for i in range(0, n):
        plt.vlines(x=distribution["mass"][i], ymin=0, ymax=distribution["abundance"][i])
        if int(distribution["mass"][i - 1]) == int(distribution["mass"][i])                 and i != 0:
            adjustText(distribution["mass"][i - 1], distribution["abundance"][i - 1],
                       distribution["mass"][i], distribution["abundance"][i])
        else:
            plt.text(x=distribution["mass"][i],
                     y=(distribution["abundance"][i] + 2),
                     s='%0.3f' % (distribution["abundance"][i]), va='center',
                     ha='center')
    plt.ylim([0, 110])
    plt.xticks(range(math.ceil(distribution["mass"][0]) - 2,
                     math.ceil(distribution["mass"][-1]) + 2))


plt.figure(figsize=(10,7))

plt.subplot(1,2,1)
plt.title("Isotopic distribution of nitrogen")
plotDistribution(nitrogen_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.subplot(1,2,2)
plt.title("Isotopic distribution of hydrogen")
plotDistribution(hydrogen_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.show()


# In[36]:


methanol = EmpiricalFormula("CH3OH")
water = EmpiricalFormula("H2O")
ethanol = EmpiricalFormula("CH2") + methanol
print("Ethanol chemical formula:", ethanol.toString())
print("Ethanol composition:", ethanol.getElementalComposition())
print("Ethanol has", ethanol.getElementalComposition()[b"H"], "hydrogen atoms")


# In[38]:


methanol = EmpiricalFormula("CH3OH")
ethanol = EmpiricalFormula("CH2") + methanol

methanol_isoDist = {"mass": [], "abundance": []}
ethanol_isoDist = {"mass": [], "abundance": []}

print("Coarse Isotope Distribution:")
isotopes = ethanol.getIsotopeDistribution( CoarseIsotopePatternGenerator(4) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("This covers", prob_sum, "probability")
for iso in isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    methanol_isoDist["mass"].append(iso.getMZ())
    methanol_isoDist["abundance"].append((iso.getIntensity() * 100))

print("Fine Isotope Distribution:")
isotopes = ethanol.getIsotopeDistribution( FineIsotopePatternGenerator(1e-3) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("This covers", prob_sum, "probability")
for iso in isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")
    ethanol_isoDist["mass"].append(iso.getMZ())
    ethanol_isoDist["abundance"].append((iso.getIntensity() * 100))


# In[40]:


plt.figure(figsize=(10,7))

plt.subplot(1,2,1)
plt.title("Isotopic distribution of methanol")
plotDistribution(methanol_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.subplot(1,2,2)
plt.title("Isotopic distribution of ethanol")
plotDistribution(ethanol_isoDist)
plt.xlabel("Atomic mass (u)")
plt.ylabel("Relative abundance (%)")

plt.savefig("methanol_ethanol_isoDistribution.png")


# In[42]:


print("Fine Isotope Distribution:")
isotopes = ethanol.getIsotopeDistribution( FineIsotopePatternGenerator(1e-6) )
prob_sum = sum([iso.getIntensity() for iso in isotopes.getContainer()])
print("This covers", prob_sum, "probability")
for iso in isotopes.getContainer():
    print ("Isotope", iso.getMZ(), "has abundance", iso.getIntensity()*100, "%")


# In[52]:


ala = ResidueDB().getResidue("Alanine")
print(ala.getName())
print(ala.getThreeLetterCode())
print(ala.getOneLetterCode())
print(ala.getAverageWeight())
print(ala.getMonoWeight())
print(ala.getPka())
print(ala.getFormula().toString())


# In[55]:


ox = ModificationsDB().getModification("Oxidation")
print(ox.getUniModAccession())
print(ox.getUniModRecordId())
print(ox.getDiffMonoMass())
print(ox.getId())
print(ox.getFullId())
print(ox.getFullName())
print(ox.getDiffFormula())
isotopes = ox.getDiffFormula().getIsotopeDistribution(CoarseIsotopePatternGenerator(5))
for iso in isotopes.getContainer():
    print (iso.getMZ(), ":", iso.getIntensity())


# In[64]:


adenosine = RibonucleotideDB().getRibonucleotide(b"A")
print(adenosine.getName())
print(adenosine.getCode())
print(adenosine.getAvgMass())
print(adenosine.getMonoMass())
print(adenosine.getFormula().toString())
print(adenosine.isModified())
methyladenosine = RibonucleotideDB().getRibonucleotide(b"m1A")
print(methyladenosine.getName())
print(methyladenosine.isModified())

