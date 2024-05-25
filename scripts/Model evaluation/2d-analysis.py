from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem import rdMolAlign

from joblib import Parallel, delayed
import pandas as pd
import numpy as np
import os, sys
import re

# script for calculate the 2D-metrics
# input file: a csv file (starting points fragments---original molecules---generated molecules)

# laod data
filename = "random_L_site_beam250_1000.csv"
data = pd.read_csv(filename)
#print(data.columns.values) # ['tgt' 'linker' 'linker_sites' 'frags_site' 'frags' 'gens']
tgt = list(data["molecules"])
frags_site = list(data["frags_site"])
linker_sites = list(data["linker_sites"])
linker = list(data["linker"])
frags = list(data["frags"])
gens = list(data["gens"])

# get_linker from generated molecuels
def get_linker(gen, frags):
    m = Chem.MolFromSmiles(gen)
    match = m.GetSubstructMatch(Chem.MolFromSmiles(frags))
    atoms = m.GetNumAtoms()
    atoms_list = list(range(atoms))
    for i in match:
        atoms_list.remove(i)

    linker_list = atoms_list.copy()

    for i in atoms_list:
        atom = m.GetAtomWithIdx(i)
        for j in atom.GetNeighbors():
            linker_list.append(j.GetIdx())

    linker_list = list(set(linker_list))
    sites = list(set(linker_list).difference(set(atoms_list)))

    bonds = []
    for i in sites:
        atom = m.GetAtomWithIdx(i)
        for j in atom.GetNeighbors():
            if j.GetIdx() in atoms_list:
                b = m.GetBondBetweenAtoms(i, j.GetIdx())
                bonds.append(b.GetIdx())
    bonds = list(set(bonds))

    bricks = Chem.FragmentOnBonds(m, bonds)  # dummyLabels=labels
    smi = Chem.MolToSmiles(bricks)

    pattern = re.compile(r"\[.*?\]")
    for s in smi.split("."):
        count = pattern.findall(s)
        if len(count) == 2:
            s = s.replace(count[0], "[*]")
            linker_smi = s.replace(count[1], "[*]")

    return linker_smi

# part1
# 2D metrics -- Check match (if two frags all in generated molecules) and valid smiles (if the gen mol is available)
# ['tgt' 'linker' 'linker_sites' 'frags_site' 'frags' 'gens']
matches = []
linker_gen = []
linker_tgt = linker_sites
valid = len(gens) # if smiles is valid smiles

for f, f_s, t, g, l, l_s in zip(frags, frags_site, tgt, gens,linker, linker_sites):
    try:
        m = Chem.MolFromSmiles(g)
        if len(Chem.MolFromSmiles(g).GetSubstructMatch(Chem.MolFromSmiles(f)))>0:
            matches.append([f, f_s, t, g, l, l_s])
            linker_gen.append(get_linker(g, f))
    except:
        valid = valid - 1
n_matches = len(matches)
N = len(gens)
print("Number of generated SMILES: \t%d" % len(gens))
print("Number of valid SMILES: \t%d" % valid)
print("Number of match SMILES: \t%d" % n_matches)
print("The ratio of matched gens: %.3f%%" % (n_matches/N))
print("\n")
# novel linker count

n_linker_tgt = len(set(linker_sites))
n_linker_gen = len(set(linker_gen))
unique_linker = [i for i in linker_gen if i not in linker_sites]
n_novel_linker = len(unique_linker)
print("\n")
print("Number of tgt linker: \t%d" % n_linker_tgt)
print("Number of gen unique linker: \t%d" % n_linker_gen)
print("Number of novel linker: \t%d" % n_novel_linker)
print("The ratio of novel linkers: %.3f%%" % (n_novel_linker/N))
print("\n")
# unique generated molecules
# Unique identifier - starting fragments and original molecule

identifiers = []
for m in matches: # unique 在matches中计算
    #f, f_s, t, g, l, l_s
    f_s = m[1]
    t = m[2]
    g = m[3]
    #identifier = f + "." + t //n_best * count of mmps
    identifiers.append(f_s + "." + t + "." + g)
    n_set_identifier = len(list(set(identifiers)))

    unique = n_set_identifier/len(identifiers) * 100
print("\n")
print("Number of tgt data: \t%d" % len(identifiers))
print("Number of gen unique data: \t%d" % n_set_identifier)
print("The ratio of unique molecules: %.3f%%" % unique)
print("\n")

# recovered generated molecules vs tgt molecules

# Unique identifier - starting fragments and tgt molecule vs starting fragments and gen molecule

count = 0
for m in matches: # unique 在matches中计算
    #f, f_s, t, g, l, l_s
    #f_s = m[1]
    t = m[2]
    g = m[3]
    #id_t.append(f_s + "." + t)
    #id_g.append(f_s + "." + g)
    #for i in range(len(id_t)):
    if t == g:
        count = count + 1
recover = count/len(matches)
print("\n")
print("The number of repeated count: \t%d" % count)
print("The ratio of recovered: %.3f%%" % recover)
print("\n")
