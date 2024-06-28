import numpy as np
from collections import defaultdict
from ase.io import read, write
from ase import Atoms
from ase.spacegroup import get_spacegroup
import os
from CifFile import CifFile
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

class DisorderUnraveller:
    def __init__(self, input_file):
        self.input_file = input_file
        self.struct = Structure.from_file(input_file)
        self.cif_file = CifFile(input_file)
        self.system_name = self.cif_file.keys()[0]
        system_name = self.system_name
        self.cif = self.cif_file[system_name]

    def get_lattice_parameters(self):
        a, b, c, alpha, beta, gamma = \
            float(self.cif["_cell_length_a"].split('(')[0]), \
            float(self.cif["_cell_length_b"].split('(')[0]), \
            float(self.cif["_cell_length_c"].split('(')[0]), \
            np.radians(float(self.cif["_cell_angle_alpha"].split('(')[0])), \
            np.radians(float(self.cif["_cell_angle_beta"].split('(')[0])), \
            np.radians(float(self.cif["_cell_angle_gamma"].split('(')[0]))
            
        lattice_matrix = np.array([
            [a, b * np.cos(gamma), c * np.cos(beta)],
            [0, b * np.sin(gamma), c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)],
            [0, 0, c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))]
        ])
        return {
            "lattice_params": [a, b, c, np.degrees(alpha), np.degrees(beta), np.degrees(gamma)],
            "lattice_matrix": lattice_matrix
        }

    def get_occupancies(self):
        occupancies = []
        for idx, ss in enumerate(self.struct):
            occ = self.struct[idx].species.num_atoms
            occupancies.append(occ)
        return occupancies

    def get_symbols(self):
        symbols = []
        for idx, ss in enumerate(self.struct):
            sym = self.struct[idx].species.chemical_system
            symbols.append(sym)
        return symbols

    def get_atomic_coords(self):
        positions = []
        for idx, ss in enumerate(self.struct):
            pos = self.struct[idx].coords
            positions.append(pos)
        return positions

    def get_spacegroup_from_cif(self):
        return self.cif["_space_group_IT_number"]

    def group_elements(self, my_list):
        result = defaultdict(list)
        indices = defaultdict(list)
        for i, item in enumerate(my_list):
            result[item].append(item)
            indices[item].append(i)
        grouped_elements = list(map(list, result.values()))
        grouped_indices = list(indices.values())
        return grouped_elements, grouped_indices

    def unravel_disorder(self):
        symbols = self.get_symbols()
        occ = self.get_occupancies()
        coords = self.get_atomic_coords()
        splitted_occ, indices = self.group_elements(occ)
        order_part_idx = []
        disorder_parts_indices = []
        for ii, sublist in enumerate(splitted_occ):
            if 1.0 in sublist:
                order_part_idx = indices[ii]
            else:
                disorder_parts_indices.append(indices[ii])
        return {
            "order": order_part_idx,
            "disorder": disorder_parts_indices
        }

    def modify_struct(self):
        order = self.unravel_disorder()["order"]
        disorder = self.unravel_disorder()["disorder"]
        sg = int(self.get_spacegroup_from_cif())
        ase_ss = []
        for idx, frag in enumerate(disorder):
            fragment_idx = frag + order
            pos = [self.get_atomic_coords()[ii] for ii in fragment_idx]
            ele = [self.get_symbols()[ii] for ii in fragment_idx]
            ss = Atoms(ele, pos)
            ss.set_cell(self.get_lattice_parameters()["lattice_params"])
            ase_ss.append(ss)
        return ase_ss

    def write_structure(self, output_filename):
        ase_ss = self.modify_struct()
        sys_name = self.input_file.split('.cif')[0]
        for ii, ss in enumerate(ase_ss):
            directory = f'{sys_name}_disorder_unravelled/{output_filename}_replica_{ii}'

            if not os.path.exists(directory):
                os.makedirs(directory)

            output_file = f'{directory}/POSCAR'
            ss.write(output_file)
            print(f"Replica {ii} has been written to {output_file}")