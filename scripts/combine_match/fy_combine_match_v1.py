import argparse
import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MolStandardize
from rdkit.Chem import rdMolAlign
#from mmps_utils import *
from onmt.utils.logging import logger

def remove_dummys(smi_string):
    try:
        smi = Chem.MolToSmiles(Chem.RemoveHs(AllChem.ReplaceSubstructs(Chem.MolFromSmiles(smi_string), \
                                                                    Chem.MolFromSmiles('*'), \
                                                                    Chem.MolFromSmiles('[H]'), True)[0]))
    except:
        smi = ""

    return smi
#
def load_data(src_file, predfile, beam_size):
    f = open(predfile, "r")
    gens = f.readlines()
    gens = [i.split("\n")[0] for i in gens]
    gens = ["".join(i.split(" ")) for i in gens]

    #
    f2 = open(src_file, "r")
    srcs = f2.readlines()
    srcs = [i.split("\n")[0] for i in srcs]
    #srcs = ["".join(i.split(" ")) for i in srcs]
    #
    srcs_beam = [src for src in srcs for i in range(beam_size)]

    d = dict(zip(srcs_beam, gens))
    keys = list(d.keys())
    #
    for key in keys:
        if d[key] == "":
            del d[key]

    logger.info("Totally %d molecules to filter "%(len(d)))

    return d

#
def get_all(gen, frags):
    "input generated molecules and the starting fragments of original molecules \
     return to the generated linker and  the two linker sites in fragments"

    m = Chem.MolFromSmiles(gen)
    matches = m.GetSubstructMatches(Chem.MolFromSmiles(frags))
    atoms = m.GetNumAtoms()
    smis = []
    for match in matches:
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
        smis.append(Chem.MolToSmiles(bricks))
        
        smis = [s for s in smis if len(s.split(".")) == 3]
        
        tuples = []
        for smi in smis:
            
            pattern = [re.compile(r"[\[][0-9][0-9][*][\]]"), re.compile(r"[\[][0-9][*][\]]")]
            frags_smi = []
            for s in smi.split("."):
                count = [j for p in pattern for j in p.findall(s)]
                
                if len(count) == 2:
                    s = s.replace(count[0], "[*]")
                    linker_smi = s.replace(count[1], "[*]")
                else:
                    frags_smi.append(s.replace(count[0], "[*]"))
            frags_site = frags_smi[0] + "." + frags_smi[1]
            
            tuples.append([linker_smi, frags_site])

            return tuples


def dummy(frag_sites):
    frag_sites_list = []
    frag_sites_list.append(frag_sites)
    frag_sites_list.append(frag_sites.split(".")[1] + "." + frag_sites.split(".")[0])
    frag_sites_list = [Chem.MolToSmiles(Chem.MolFromSmiles(f)) for f in frag_sites_list]
    frags = remove_dummys(frag_sites)

    return frags, frag_sites_list

    return

def match(g, f, frag_sites_list):
    print("list ^^^^^^^^^^^",frag_sites_list)
    #linker, frags = get_all(g, f)
    tuples = get_all(g, f)
    
    matches = []
    for t in tuples:
        frags_mol = Chem.MolFromSmiles(t[1])   # get the two terminal ,t[0] is the linker
        smi = Chem.MolToSmiles(frags_mol)
        #for ff in frag_sites_list:
         #   if smi
        

        if smi in frag_sites_list:
            m = "match"
        else:
            m = ""
        matches.append(m)
    matches = [m for m in matches if m !=""]
    
    return matches


def save_match(smis, filename):
    f = open(filename, "w")
    smis = list(set(smis))
    for k in smis:
        f.write(k)
        f.write("\n")
    f.close()

    return 0





def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-gens_file", "-i", required=True,
                        help="File of the generative molecules")

    parser.add_argument("-src_file", "-src", required=True,  
                        help="The original source file")   # re_frags_2.py outputfile, output the two terminal pairs

    parser.add_argument("-beam_size", "-bs", required=True,
                        help="The original source file", type=int)

    parser.add_argument("-output", "-o", required=True,
                        help="Output file (matched molecules)")

    parser.add_argument("-contain", "-c", required=True,
                        help=" yes, molecules contail both two terminals, not only terminal matching molecues")
    opt = parser.parse_args()

    # main function
    d = load_data(opt.src_file, opt.gens_file, opt.beam_size)
    match_smi = []
    for k in d.keys(): # k example:L_11 * C C C C C C N ( C C c 1 c * c 2 c c c c c 1 2 ) C c 1 c c c ( C = C C ( = O ) N O ) c c 1 . * c 1 c c c ( C = C C ( = O ) N c 2 c c c c c 2 N ) c c 1
        src = k
        print('KKKKKKKKKKKKKK,    ',k)
        SLBD = src.split(" ")[0]
        frag_sites = "".join(src.split(" ")[1:])  # two terminal is the frag_sites,example:*CCCCCCN(CCc1c*c2ccccc12)Cc1ccc(C=CC(=O)NO)cc1.*c1ccc(C=CC(=O)Nc2ccccc2N)cc1
        
        f_sites = dummy(frag_sites)[1]
        f = dummy(frag_sites)[0]

        gen = d[k]  # gen is the predict molecues 


        if opt.contain == "yes":
            frag_sites_1 = frag_sites.split(".")[0] # first terminal
            frag_sites_2 = frag_sites.split(".")[1] # second terninal  
            
            if frag_sites_1 in gen and frag_sites_2 in gen 
        else:
            try:
                m = match(gen, f, f_sites)
                if "match" in m:
                    match_smi.append(SLBD + "," + f + "," + gen)
            except:
                print("gen molfile error!")
                continue
            

#        try:
#            m = match(gen, f, f_sites)
#            if "match" in m:
#                match_smi.append(SLBD + "," + f + "," + gen)
#        except:
#            print("gen molfile error!")
#            continue

    logger.info("Totally %d molecules remaining " % (len(match_smi)))
    save_match(match_smi, opt.output)

if __name__ == "__main__":
    main()
