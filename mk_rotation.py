from pymatgen.core.structure import Structure
import numpy as np
import os
import copy
from mol_identifier import MoleculeIdentifier
from chk_bonding import ChkBonding

class MoleculeRotator:
    def __init__(self, 
                 input_file, 
                 angle, 
                 rotate_axes, 
                 num_replicas,
                 tol_factor=1.1, 
                 intra_crit_dist=0.60,
                 intra_crit_lo_ratio=0.95, 
                 intra_crit_hi_ratio=1.05
                 ):
        self.input_file = input_file
        self.struct = Structure.from_file(input_file)
        self.angle = angle
        self.rotate_axes = rotate_axes
        self.num_replicas = num_replicas
        self.mol_identifier = MoleculeIdentifier(tol_factor)
        self.chk_bonding = ChkBonding(self.mol_identifier, intra_crit_dist, intra_crit_lo_ratio, intra_crit_hi_ratio)
        self.deg_to_rad = np.pi / 180

    def get_anchor(self, mol_site_idx):
        positions = []
        for idx in mol_site_idx:
            pos = self.struct[idx].coords
            positions.append(pos)
        anchor = np.array(positions).mean(axis=0)
        return anchor

    def rot_mol(self, mol_idx, mol_site_idx, axis, anchor):
        for aa in range(int(-1 * self.angle), int(self.angle + 1), int(self.angle / (self.num_replicas / 2))):
            if aa == 0:
                continue
            angle = aa * self.deg_to_rad
            new_struct = copy.deepcopy(self.struct)
            new_struct.rotate_sites(mol_site_idx, angle, axis, anchor, to_unit_cell=False)
            file_name = self.input_file.split('.cif')[0]
            axis_str = "".join(str(x) for x in axis)
            directory = f"{file_name}_rot/mol_{mol_idx}__axis_{axis_str}__rot_{aa}deg"
            print(directory, self.chk_bonding.pass_checking(new_struct, self.molecules_site_indices))
            if self.chk_bonding.pass_checking(new_struct, self.molecules_site_indices):
                if not os.path.exists(directory):
                    os.makedirs(directory)
                new_struct.to(f"{directory}/POSCAR")

    def rotate_molecules(self):
        self.molecules_site_indices = self.mol_identifier.find_molecules(self.struct)
        for mol_idx, mol_site_idx in enumerate(self.molecules_site_indices):
            anchor = self.get_anchor(mol_site_idx)
            for axis in self.rotate_axes:
                self.rot_mol(mol_idx, mol_site_idx, axis, anchor)