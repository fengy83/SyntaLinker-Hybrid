import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from rdkit import Chem
from rdkit.Chem import MACCSkeys
from sklearn.decomposition import PCA
import csv
#%matplotlib inline

def fp2arr(fp):
    from rdkit import DataStructs
    arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def read_csv(path, line):
    '''
    path is str
    line is list
    return pd
    '''
    csv_data = pd.read_csv(path)
    need_data = pd.DataFrame(csv_data, columns=line)
    return need_data

file_name ='./chembl_transf_pca_100+.csv'
mols = []
labels = []
allsmiles = []
save_csv = []
oldsmiles = read_csv(file_name, ['SMILES','LABEL'])
for index, (smiles, label) in enumerate(oldsmiles.values):
    if index % 10000 == 0:
        print(index)
    try:
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
        mols.append(mol)
    except:
        continue
    allsmiles.append(smiles)
    labels.append(label)
print('read done')
fps = [MACCSkeys.GenMACCSKeys(mol) for mol in mols]
fpMtx = np.array([fp2arr(fp) for fp in fps])
pca = PCA(n_components=2)
pc = pca.fit_transform(fpMtx)
print('pca done')
red_x,red_y=[],[]
blue_x,blue_y=[],[]
green_x,green_y=[],[]
orange_x,orange_y=[],[]
for i in range(len(pc)):
    save_acsv = []
    save_acsv.append(allsmiles[i])
    save_acsv.append(labels[i])
    save_acsv.append(pc[i][0])
    save_acsv.append(pc[i][1])
    save_csv.append(save_acsv)
    if labels[i]==0:
        green_x.append(pc[i][0])
        green_y.append(pc[i][1])
    elif labels[i]==1:
        blue_x.append(pc[i][0])
        blue_y.append(pc[i][1])
    elif labels[i]==2:
        red_x.append(pc[i][0])
        red_y.append(pc[i][1])
    elif labels[i]==3:
        orange_x.append(pc[i][0])
        orange_y.append(pc[i][1])
print('list done')
with open('./test.csv', 'w+', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['smiles','label','pc1','pc2']) 
    writer.writerows(save_csv)
print('flie done')

plt.figure(figsize=(10, 10), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend()
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.savefig("./test_all.png")