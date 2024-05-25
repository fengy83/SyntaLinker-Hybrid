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

file_name ='./test.csv'
red_x,red_y=[],[]
blue_x,blue_y=[],[]
green_x,green_y=[],[]
orange_x,orange_y=[],[]
oldsmiles = read_csv(file_name, ['label','pc1','pc2'])
for index, (label,pc1,pc2) in enumerate(oldsmiles.values):
    if labels==0:
        green_x.append(pc1)
        green_y.append(pc2)
    elif labels==1:
        blue_x.append(pc1)
        blue_y.append(pc2)
    elif labels==2:
        red_x.append(pc1)
        red_y.append(pc2)
    elif labels==3:
        orange_x.append(pc1)
        orange_y.append(pc2)

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_all.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_dataset.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_cd.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_xd.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_cl.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_cx.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_2_cdcl.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_2_clcx.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_2_xdcx.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_3_cdxdxl.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_3_cdxdcx.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(green_x,green_y,c='g', marker=".", s=1, label='chembl_dataset')
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_3_cdclcx.tiff")

plt.figure(figsize=(8, 8), dpi=600)
plt.scatter(blue_x,blue_y,c='b', marker=".", s=2, label='xtal_dataset')
plt.scatter(red_x,red_y,c='r', marker=".", s=2, label='chembl_test')
plt.scatter(orange_x,orange_y,c='orange', marker=".", s=2, label='chembl_xtal_transfer_test')
plt.legend(fontsize=10)
plt.xlabel("PC1",fontsize=10)
plt.ylabel("PC2",fontsize=10)
plt.tick_params('both',width=2,labelsize=8)
plt.savefig("./test_3_xdclcx.tiff")








