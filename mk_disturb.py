import os
import random
import numpy as np
import copy
from pymatgen.core import Structure
from mol_identifier import MoleculeIdentifier
from chk_bonding import ChkBonding

class StructurePerturber:
    def __init__(
        self,
        input_file,
        perturb,
        scale_percentage,
        max_tilt,
        num_replicas,
        tol_factor=1.1,
        intra_crit_dist=0.60,
        intra_crit_lo_ratio=0.95,
        intra_crit_hi_ratio=1.05
    ):
        self.input_file = input_file
        self.orig_structure = Structure.from_file(input_file)
        self.perturb = perturb
        self.scale_percentage = scale_percentage
        self.max_tilt = max_tilt
        self.num_replicas = num_replicas
        self.mol_identifier = MoleculeIdentifier(tol_factor)
        self.chk_bonding = ChkBonding(self.mol_identifier, intra_crit_dist, intra_crit_lo_ratio, intra_crit_hi_ratio)

    def perturb_atoms(self):
        sys_name = self.input_file.split('.cif')[0]
        self.molecules_site_indices = self.mol_identifier.find_molecules(self.orig_structure)
        for ii in range(self.num_replicas):
            # Create a copy of the original structure
            new_struct = copy.deepcopy(self.orig_structure)
            new_struct.perturb(self.perturb)
            if self.chk_bonding.pass_checking(new_struct, self.molecules_site_indices):
                directory = f'{sys_name}_disturbed/{sys_name}_perturb_{self.perturb}_{ii}'
                if not os.path.exists(directory):
                    os.makedirs(directory)
                output = f'{directory}/POSCAR'
                new_struct.to(output)
                print(f'Done making perturbations and output to {output}')
            else:
                print(f"Replica {ii} of atomic perturbing not passed.")

    def scale_cells(self):
        volume = self.orig_structure.volume
        sys_name = self.input_file.split('.cif')[0]
        self.molecules_site_indices = self.mol_identifier.find_molecules(self.orig_structure)
        step_size = (2 * self.scale_percentage) / (self.num_replicas - 1)
        for ii in range(self.num_replicas):
            scale = 1.0 - self.scale_percentage + ii * step_size
            print(scale)
            # Create a copy of the original structure
            new_struct = copy.deepcopy(self.orig_structure)
            new_struct.scale_lattice(volume * scale)
            if self.chk_bonding.pass_checking(new_struct, self.molecules_site_indices):
                directory = f'{sys_name}_disturbed/{sys_name}_scale_{scale:.2f}'
                if not os.path.exists(directory):
                    os.makedirs(directory)
                output = f'{directory}/POSCAR'
                new_struct.to(output)
                print(f'Done making scales and output to {output}')
            else:
                print(f"Replica of cell scaling at {scale:.2f} not passed.")

    def tilt_cells(self):
        sys_name = self.input_file.split('.cif')[0]
        lattice_matrix = self.orig_structure.lattice.matrix
        self.molecules_site_indices = self.mol_identifier.find_molecules(self.orig_structure)
        for ii in range(self.num_replicas):
            tilt_xy = random.uniform(-self.max_tilt, self.max_tilt)
            tilt_xz = random.uniform(-self.max_tilt, self.max_tilt)
            tilt_yx = random.uniform(-self.max_tilt, self.max_tilt)
            tilt_yz = random.uniform(-self.max_tilt, self.max_tilt)
            tilt_zx = random.uniform(-self.max_tilt, self.max_tilt)
            tilt_zy = random.uniform(-self.max_tilt, self.max_tilt)
            tilt_matrix = np.array([[1, tilt_xy, tilt_xz],
                                    [tilt_yx, 1, tilt_yz],
                                    [tilt_zx, tilt_zy, 1]])
            new_lattice = lattice_matrix @ tilt_matrix
            # Create a copy of the original structure
            new_struct = copy.deepcopy(self.orig_structure)
            new_struct.lattice = new_lattice
            if self.chk_bonding.pass_checking(new_struct, self.molecules_site_indices):
                directory = f'{sys_name}_disturbed/{sys_name}_tilt_{ii}'
                if not os.path.exists(directory):
                    os.makedirs(directory)
                output = f'{directory}/POSCAR'
                new_struct.to(output)
                print(f'Done making tilt and output to {output}')
            else:
                print(f"Replica {ii} of cell tilting not passed.")