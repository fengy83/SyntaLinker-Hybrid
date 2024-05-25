
import csv
import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import MACCSkeys

def read_csv(path, line):
    '''
    path is str
    line is list
    return pd
    '''
    csv_data = pd.read_csv(path)
    need_data = pd.DataFrame(csv_data, columns=line)
    return need_data

file_name ='./scaffold_remove.csv'
apfrag = []
linker = []
allsmiles = read_csv(file_name, ['SMILES','LABEL'])
for smiles, label in allsmiles.values:
    if label == 0:
        apfrag.append(smiles)
    elif label == 1:
        linker.append(smiles)
print(len(apfrag),len(linker))


save_date = []
for index, apf in enumerate(apfrag):
    print(index)
    a_data = []
    label_smi = ['AP']
    try:
        #fp1 = Chem.RDKFingerprint(Chem.MolFromSmiles(apf))
        fp1 = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(apf))
        a_data.append(apf)
    except:
        continue
    for lin in linker:
        try:
            #fp2 = Chem.RDKFingerprint(Chem.MolFromSmiles(lin))
            fp2 = MACCSkeys.GenMACCSKeys(Chem.MolFromSmiles(lin))
            simi = DataStructs.FingerprintSimilarity(fp1, fp2)
            a_data.append(simi)
            label_smi.append(lin)
        except:
            continue
    load_label = label_smi
    save_date.append(a_data)
    
    
    

with open('./similarity_frag.csv', 'w+', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(load_label)
    writer.writerows(save_date)
print('flie done')