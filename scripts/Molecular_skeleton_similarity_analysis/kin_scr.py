import csv
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
from rdkit.Chem.Scaffolds.MurckoScaffold import MurckoScaffoldSmiles

def remove_dummys(smi_string):
    return Chem.MolToSmiles(Chem.RemoveHs(AllChem.ReplaceSubstructs(Chem.MolFromSmiles(smi_string), \
                                                                    Chem.MolFromSmiles('*'), \
                                                                    Chem.MolFromSmiles('[H]'), True)[0]))

data = pd.read_csv('./KIN-SE.csv',names = ['fragment_ids','bond_ids','subpocket_key','smiles','ap_frag'])
fragment_ids = data['fragment_ids']
bond_ids = data['bond_ids']
subpocket_key = data['subpocket_key']
smiles = data['smiles']
ap_frag = data['ap_frag']

save_data = []
for index, a_ap in enumerate(ap_frag):
    a_data = []
    try:
        ap_sca = MurckoScaffoldSmiles(remove_dummys(a_ap))
        kin_sca = MurckoScaffoldSmiles(smiles[index])
        a_data.append(fragment_ids[index])
        a_data.append(bond_ids[index])
        a_data.append(subpocket_key[index])
        a_data.append(Chem.MolToSmiles(Chem.MolFromSmiles(smiles[index])))
        a_data.append(Chem.MolToSmiles(Chem.MolFromSmiles(a_ap)))
        a_data.append(ap_sca)
        a_data.append(kin_sca)
        save_data.append(a_data)
    except:
        continue
    if index % 10000 == 0:
        print(index)

with open('./KIN-SE_sca.csv', 'w+', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['fragment_ids','bond_ids','subpocket_key','smiles','ap_frag','ap_sca','kin_sca'])
    writer.writerows(save_data)
print('done')

