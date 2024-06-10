from pymatgen.core.structure import Structure
import numpy as np
import os
import copy
from mol_identifier import MoleculeIdentifier

class MoleculeRotator:
    def __init__(self, input_file, angle, rotate_axes, num_replicas, radii_file,
                 pm_to_angstrom=0.01, tol_factor=1.1, intra_crit_dist=0.60,
                 intra_crit_lo_ratio=0.95, intra_crit_hi_ratio=1.05):
        self.input_file = input_file
        self.struct = Structure.from_file(input_file)
        self.angle = angle
        self.rotate_axes = rotate_axes
        self.num_replicas = num_replicas
        self.mol_identifier = MoleculeIdentifier(radii_file, pm_to_angstrom, tol_factor)
        self.intra_crit_dist = intra_crit_dist
        self.intra_crit_lo_ratio = intra_crit_lo_ratio
        self.intra_crit_hi_ratio = intra_crit_hi_ratio
        self.deg_to_rad = np.pi / 180

    def get_anchor(self, mol_site_idx):
        positions = []
        for idx in mol_site_idx:
            pos = self.struct[idx].coords
            positions.append(pos)
        anchor = np.array(positions).mean(axis=0)
        return anchor

    def check_dirty_struct_pass(self, new_struct, mol_site_idx):
        overall_set = set([num for sublist in self.molecules_site_indices for num in sublist])
        sel_mol_set = set(mol_site_idx)
        exclude_set = overall_set.difference(sel_mol_set)
        exclude_idx = list(exclude_set)
        sel_mol_idx = list(sel_mol_set)

        check_results = []
        for sel_idx in sel_mol_idx:
            for exl_idx in exclude_idx:
                chk_res_inter = self.mol_identifier.calc_bond(new_struct, exl_idx, sel_idx)
                check_results.append(chk_res_inter)

            for other_sel_idx in sel_mol_idx[sel_mol_idx.index(sel_idx) + 1:]:
                intra_crit_1 = new_struct.get_distance(sel_idx, other_sel_idx) < self.intra_crit_dist
                intra_crit_2 = new_struct.get_distance(sel_idx, other_sel_idx) < \
                               self.intra_crit_lo_ratio * self.struct.get_distance(sel_idx, other_sel_idx)
                intra_crit_3 = new_struct.get_distance(sel_idx, other_sel_idx) > \
                               self.intra_crit_hi_ratio * self.struct.get_distance(sel_idx, other_sel_idx)
                chk_res_intra = bool(intra_crit_1 or intra_crit_2 or intra_crit_3)
                check_results.append(chk_res_intra)

        if True in check_results:
            return False
        else:
            return True

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
            print(directory, self.check_dirty_struct_pass(new_struct, mol_site_idx))
            if self.check_dirty_struct_pass(new_struct, mol_site_idx):
                if not os.path.exists(directory):
                    os.makedirs(directory)
                new_struct.to(f"{directory}/POSCAR")

    def rotate_molecules(self):
        self.molecules_site_indices = self.mol_identifier.find_molecules(self.struct)
        for mol_idx, mol_site_idx in enumerate(self.molecules_site_indices):
            anchor = self.get_anchor(mol_site_idx)
            for axis in self.rotate_axes:
                self.rot_mol(mol_idx, mol_site_idx, axis, anchor)