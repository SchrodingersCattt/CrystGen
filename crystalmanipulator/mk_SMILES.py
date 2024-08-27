import numpy as np
from pymatgen.core.structure import Structure, Molecule
from rdkit import Chem
from multiprocessing import Pool
from tqdm import tqdm
import itertools
import copy
import time

from crystalmanipulator.mol_identifier import MoleculeIdentifier


MAX_DISTANCE_THRESHOLD = 3.5

def search_molecules_in_crystal(struct_file):
    """
    input: cif file or POSCAR
    output: tuple: (list of pymatgen Molecule, list of original site indices)
    """
    mm = MoleculeIdentifier()
    struct = Structure.from_file(struct_file)
    mols = []
    original_indices = []
    molecules_site_indices = mm.find_molecules(struct)
    
    for mol_indices in molecules_site_indices:
        imaged_sites = []
        added_indices = set()
        all_coords = np.array([struct[i].coords for i in mol_indices])
        mol_original_indices = []
        
        for ii, jj in itertools.combinations(range(len(mol_indices)), 2):
            if ii not in added_indices:
                imaged_sites.append(struct[mol_indices[ii]])
                added_indices.add(ii)
                mol_original_indices.append(mol_indices[ii])
                
            distance, image = struct[mol_indices[ii]].distance_and_image(struct[mol_indices[jj]])
            
            jj_site = struct[mol_indices[jj]]  # shallow copy
            jj_site.frac_coords += image
            
            if jj not in added_indices:
                imaged_sites.append(jj_site)
                added_indices.add(jj)
                mol_original_indices.append(mol_indices[jj])
        
        coords = [s.coords for s in imaged_sites]
        coords = np.array(coords)
        numbers = [s.specie.number for s in imaged_sites]
        mols.append(Molecule(numbers, coords))
        original_indices.append(mol_original_indices)
    
    return mols, original_indices

def mol2smiles(mol):
    mol_file = save_mol(mol)
    stringWithMolData = open(mol_file, 'r').read()
    rdkit_mol = Chem.MolFromMolBlock(stringWithMolData, sanitize=False)
    try:
        smiles = Chem.MolToSmiles(rdkit_mol)
        return smiles

    except ValueError as e:
        print(f"Error converting molecule to SMILES: {e}")
        return None

def save_mol(mol, mol_file="temp.mol"):
    mol.to(mol_file)  
    return mol_file

def calculate_molecule_com(struct, mol_indices):
    """
    Calculate the center of mass for a molecule in the original structure.
    
    :param struct: pymatgen Structure
    :param mol_indices: list of site indices for the molecule
    :return: numpy array of shape (3,) representing the center of mass
    """
    total_mass = 0
    weighted_coords = np.zeros(3)
    
    for idx in mol_indices:
        site = struct[idx]
        mass = site.specie.atomic_mass
        total_mass += mass
        weighted_coords += float(mass) * site.coords
    
    return weighted_coords / total_mass

def get_molecules_com(struct_file):
    struct = Structure.from_file(struct_file)
    mols, original_indices = search_molecules_in_crystal(struct_file)
    
    centers_of_mass = []
    for mol_indices in original_indices:
        com = calculate_molecule_com(struct, mol_indices)
        centers_of_mass.append(com)
    
    return centers_of_mass, mols, original_indices

if __name__ == '__main__':
    struct_file = "/123/guomingyu/CrystalPhilately/20240701/PAP-H5/PAP-H5_disorder_unravelled/_replica_0/POSCAR"
    mols = search_molecules_in_crystal(struct_file)
    
    for idx, mol in enumerate(mols):
        save_mol(mol)
        smiles = mol2smiles(mol)
