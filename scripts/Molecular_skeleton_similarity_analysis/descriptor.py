import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors,Crippen
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.cm as cm
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans

#Loading the molecules
data=pd.DataFrame()
#with open('scaffold_remove.csv','r') as file:
with open('AP-GA-SE.csv','r') as file:
    for index,line in enumerate(file):
        if 0<index:
            data.loc[index,'SMILES']=line.split(',')[0]
            data.loc[index,'LABEL']=int(line.split(',')[1].replace('\n',''))
        if index % 10000 == 0:
            print(index)
print('read done')

# calculate the descriptors and add them to dataframe
for i in data.index:
    if i % 10000 == 0:
        print(i)
    mol=Chem.MolFromSmiles(data.loc[i,'SMILES'])
    data.loc[i,'MolWt']=Descriptors.ExactMolWt (mol)
    data.loc[i,'TPSA']=Chem.rdMolDescriptors.CalcTPSA(mol) #Topological Polar Surface Area
    data.loc[i,'nRotB']=Descriptors.NumRotatableBonds (mol) #Number of rotable bonds
    data.loc[i,'HBD']=Descriptors.NumHDonors(mol) #Number of H bond donors
    data.loc[i,'HBA']=Descriptors.NumHAcceptors(mol) #Number of H bond acceptors
    data.loc[i,'LogP']=Descriptors.MolLogP(mol) #LogP
print('calculate done')
#print(data.head(10))

descriptors = data.loc[:, ['MolWt', 'TPSA', 'nRotB', 'HBD','HBA', 'LogP']].values
descriptors_std = StandardScaler().fit_transform(descriptors)

pca = PCA()
descriptors_2d = pca.fit_transform(descriptors_std)
descriptors_pca = pd.DataFrame(descriptors_2d)
descriptors_pca.index = data.index
descriptors_pca.columns = ['PC{}'.format(i+1) for i in descriptors_pca.columns]
descriptors_pca['LABEL'] = data['LABEL']
descriptors_pca['SMILES'] = data['SMILES']
descriptors_pca['MolWt'] = data['MolWt']
descriptors_pca['TPSA'] = data['TPSA']
descriptors_pca['nRotB'] = data['nRotB']
descriptors_pca['HBD'] = data['HBD']
descriptors_pca['HBA'] = data['HBA']
descriptors_pca['LogP'] = data['LogP']
descriptors_pca_0 = descriptors_pca[descriptors_pca['LABEL'] == 0]
descriptors_pca_1 = descriptors_pca[descriptors_pca['LABEL'] == 1]
descriptors_pca_2 = descriptors_pca[descriptors_pca['LABEL'] == 2]#
descriptors_pca.to_csv('./descriptors_pca_molecule.csv',index=0)#

'''plt.rcParams['axes.linewidth'] = 1.5
plt.figure(figsize=(8,6))'''
fig, ax = plt.subplots(figsize=(8,6))
var=np.cumsum(np.round(pca.explained_variance_ratio_, decimals=3)*100)
plt.plot([i+1 for i in range(len(var))],var,'k-',linewidth=2)
plt.xticks([i+1 for i in range(len(var))])
plt.ylabel('% Variance Explained',fontsize=16,fontweight='bold')
plt.xlabel('Pincipal Component (PC)',fontsize=16,fontweight='bold')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.tick_params ('both',width=2,labelsize=12)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)
ax.plot(descriptors_pca_1['PC1'],descriptors_pca_1['PC2'],'o',color='b')
ax.plot(descriptors_pca_0['PC1'],descriptors_pca_0['PC2'],'o',color='r')
ax.set_title ('Principal Component Analysis',fontsize=16,fontweight='bold',family='sans-serif')
ax.set_xlabel ('PC1',fontsize=14,fontweight='bold')
ax.set_ylabel ('PC2',fontsize=14,fontweight='bold')
plt.tick_params ('both',width=2,labelsize=12)
plt.tight_layout()
plt.show()



