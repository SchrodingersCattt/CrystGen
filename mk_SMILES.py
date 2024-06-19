import numpy as np
from pymatgen.core.structure import Structure, Molecule
from pymatgen.core.bonds import CovalentBond
from mol_identifier import MoleculeIdentifier
from rdkit import Chem
from multiprocessing import Pool
from tqdm import tqdm

struct_file = "/123/PX/px483/PX483.cif"
struct = Structure.from_file(struct_file)
crys_coords = struct.cart_coords
crys_lattice = struct.lattice
crys_species = struct.species
mm = MoleculeIdentifier()

def search_molecules_in_crystal(struc, tol=0.2, once=False, ignore_HH=True):
    def check_one_layer(struc, sites0, visited):
        new_members = []
        for site0 in sites0:
            sites_add, visited = check_one_site(struc, site0, visited)
            new_members.extend(sites_add)
        return new_members, visited

    def check_one_site(struc, site0, visited, rmax=2.8):
        neigh_sites = struc.get_neighbors(site0, rmax)
        ids = [m.index for m in visited]
        sites_add = []
        ids_add = []
        pbc = isinstance(struc, Structure)
        for site1 in neigh_sites:
            if site1.index not in ids + ids_add:
                try:
                    if CovalentBond.is_bonded(site0, site1, tol):
                        if pbc:
                            (d, image) = site0.distance_and_image(site1)
                        else:
                            d = site0.distance(site1)
                        val1, val2 = site1.specie.value, site0.specie.value
                        key = "{:s}-{:s}".format(val1, val2)
                        if key == 'H-H':
                            if not ignore_HH:
                                if pbc:
                                    site1.frac_coords += image
                                sites_add.append(site1)
                                ids_add.append(site1.index)
                        else:
                            # Calculate bond distance using atomic radii
                            r1 = mm.get_radii()[site1.specie.value]
                            r2 = mm.get_radii()[site0.specie.value]
                            bond_dist = r1 + r2 + tol
                            if d < bond_dist:
                                if pbc:
                                    site1.frac_coords += image
                                sites_add.append(site1)
                                ids_add.append(site1.index)
                except ValueError:
                    # Calculate bond distance using atomic radii
                    if pbc:
                        (d, image) = site0.distance_and_image(site1)
                    else:
                        d = site0.distance(site1)
                    val1, val2 = site1.specie.value, site0.specie.value
                    r1 = mm.get_radii()[val1]
                    r2 = mm.get_radii()[val2]
                    bond_dist = r1 + r2 + tol
                    if d < bond_dist:
                        if pbc:
                            site1.frac_coords += image
                        sites_add.append(site1)
                        ids_add.append(site1.index)
        if len(sites_add) > 0:
            visited.extend(sites_add)
        return sites_add, visited

    molecules = []
    visited_ids = []
    for id, site in enumerate(struc.sites):
        if id not in visited_ids:
            first_site = site
            visited = [first_site]
            first_site.index = id
            n_iter, max_iter = 0, len(struc) - len(visited_ids)
            while n_iter < max_iter:
                if n_iter == 0:
                    new_sites, visited = check_one_site(struc, first_site, visited)
                else:
                    new_sites, visited = check_one_layer(struc, new_sites, visited)
                n_iter += 1
                if len(new_sites) == 0:
                    break
            coords = [s.coords for s in visited]
            coords = np.array(coords)
            numbers = [s.specie.number for s in visited]
            molecules.append(Molecule(numbers, coords))
            visited_ids.extend([s.index for s in visited])
            # print(molecules[-1].to(fmt='xyz')); import sys; sys.exit()
        if once and len(molecules) == 1:
            break
    return molecules


    molecules_site_indices = find_molecules(struc)
    molecules = []
    for molecule_site_indices in molecules_site_indices:
        coords = [struc.sites[i].coords for i in molecule_site_indices]
        coords = np.array(coords)
        numbers = [struc.sites[i].specie.number for i in molecule_site_indices]
        molecules.append(Molecule(numbers, coords))
        if once and len(molecules) == 1:
            break

    return molecules

mols = search_molecules_in_crystal(struct)


for idx, mol in enumerate(mols):   

    mol.to("tmp.mol")
    stringWithMolData = open("tmp.mol", 'r').read()
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