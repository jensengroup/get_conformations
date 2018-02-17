# written by Jan Jensen
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
import itertools
import random

def get_subset(number_of_conformations,max_number_of_conformations,angle_combinations):
    subset = []
    subset_index = random.sample(range(number_of_conformations), max_number_of_conformations)
    for index in subset_index:
        subset.append(angle_combinations[index])
        
    return subset

def get_conformations(mol,max_number_of_conformations):
    raw_rot_bonds =  mol.GetSubstructMatches(Chem.MolFromSmarts("[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]"))
    raw_rot_bonds += mol.GetSubstructMatches(Chem.MolFromSmarts("[*]~[*]-[O,S]-[#1]"))
    raw_rot_bonds += mol.GetSubstructMatches(Chem.MolFromSmarts("[*]~[*]-[NX3;H2]-[#1]"))

    bonds = []
    rot_bonds = []
    for k,i,j,l in raw_rot_bonds:
        if (i,j) not in bonds:
            bonds.append((i,j))
            rot_bonds.append((k,i,j,l))

    conf = mol.GetConformer()

    angles = []
    for k,i,j,l in rot_bonds:
        theta = rdMolTransforms.GetDihedralDeg( conf, k,i,j,l )
        angles.append([theta, theta+120., theta-120.])

    angle_combinations = list(itertools.product(*angles))
    
    number_of_conformations = len(angle_combinations)
    
    if number_of_conformations > max_number_of_conformations:
        angle_combinations = get_subset(number_of_conformations,max_number_of_conformations,angle_combinations)

    for angle_combination in angle_combinations:
        for (k,i,j,l),angle in zip(rot_bonds,angle_combination):
            rdMolTransforms.SetDihedralDeg( conf, k,i,j,l, angle )
        mol.AddConformer(conf,assignId=True)

    #print mol.GetNumConformers()
    
    return mol

if __name__ == "__main__":
    smiles = "C(C)(C)CCN"
    #smiles = "CCCC"
    #smiles = "C1(C2=CC=CN2)=NC=CC=C1"
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol,explicitOnly=False)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    
    max_conf = 10

    mol = get_conformations(mol,max_conf)

    w = Chem.SDWriter('conformers.sdf') 
    for i,conf in enumerate(mol.GetConformers()):
        tm = Chem.Mol(mol,False,conf.GetId())
        w.write(tm)
    
