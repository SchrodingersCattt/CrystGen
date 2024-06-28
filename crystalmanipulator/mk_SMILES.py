import numpy as np
from pymatgen.core.structure import Structure, Molecule
from rdkit import Chem
from multiprocessing import Pool
from tqdm import tqdm
import itertools
import copy
import time

from crystalmanipulator.mol_identifier import MoleculeIdentifier

struct_file = "/123/PX/px483/PX483.cif"
struct = Structure.from_file(struct_file)

#struct = Structure.from_file("/123/gmy_data/TOOLS_IN_DEV/CrystGen/example/Antipyrin_912076.cif")
crys_coords = struct.cart_coords
crys_lattice = struct.lattice
crys_species = struct.species
mm = MoleculeIdentifier()

MAX_DISTANCE_THRESHOLD = 3.8

def search_molecules_in_crystal(struct):
    ss1 = time.time()
    mols = []
    molecules_site_indices = mm.find_molecules(struct)
    ss2 = time.time()
    for mol_indices in molecules_site_indices:
        imaged_sites = []
        added_indices = set()
        all_coords = np.array([struct[i].coords for i in mol_indices])
        for ii, jj in itertools.combinations(range(len(mol_indices)), 2):
            if ii not in added_indices:
                imaged_sites.append(struct[mol_indices[ii]])
                added_indices.add(ii)
            
            distance, image = struct[mol_indices[ii]].distance_and_image(struct[mol_indices[jj]])
            
            jj_site = struct[mol_indices[jj]] ## shallow copy
            jj_site.frac_coords += image
            
            if jj not in added_indices:
                imaged_sites.append(jj_site)
                added_indices.add(jj)
        # print(imaged_sites)
        coords = [s.coords for s in imaged_sites]
        coords = np.array(coords)
        numbers = [s.specie.number for s in imaged_sites]
        mols.append(Molecule(numbers, coords))
        
    ee = time.time()
    print(ee - ss1, ee - ss2)   
    print([len(ll) for ll in mols])
    return mols

mols = search_molecules_in_crystal(struct)


for idx, mol in enumerate(mols):   

    mol.to(f"tmp.mol")
    stringWithMolData = open(f"tmp.mol", 'r').read()
    rdkit_mol = Chem.MolFromMolBlock(stringWithMolData, sanitize=False)
    # Check if rdkit_mol is None before calling Chem.MolToSmiles
    if rdkit_mol is None:
        print("Error: Invalid or empty molecule object")
    else:
        try:
            smiles = Chem.MolToSmiles(rdkit_mol)
            print(f"SMILES: {smiles}")
        except ValueError as e:
            print(f"Error converting molecule to SMILES: {e}")