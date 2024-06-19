import os
import numpy as np
from tqdm import tqdm
from pymatgen.core import Structure

class CellExpander:
    def __init__(self, input_file, supercell_scale, supercell_list):
        self.input_file = input_file
        self.structure = Structure.from_file(input_file)
        self.supercell_scale = supercell_scale
        self.supercell_list = supercell_list

    @staticmethod
    def scale_to_list(n):
        return [[x, y, z] 
                for x in [n, 1] 
                for y in [n, 1] 
                for z in [n, 1] 
                if [x, y, z] != [1, 1, 1]]

    def output_file(self, sys_name, supercell_name):
        directory = f'{sys_name}_supercell/{sys_name}_supercell_{supercell_name}'
        if not os.path.exists(directory):
            os.makedirs(directory)
        output = f'{directory}/POSCAR'
        return output 

    def expand_cells(self):
        sys_name = os.path.splitext(self.input_file)[0] 
        if self.supercell_scale is not None:
            supercell_lists = self.scale_to_list(self.supercell_scale)
            for ll in tqdm(supercell_lists):
                temp_struct = self.structure.copy()
                temp_struct.make_supercell(ll)
                supercell_name = "_".join(map(str, ll))
                output_file_path = self.output_file(sys_name, supercell_name)
                temp_struct.to(fmt="poscar", filename=output_file_path)

        if self.supercell_list is not None:
            self.structure.make_supercell(self.supercell_list)
            supercell_name = "_".join(map(str, self.supercell_list))  
            output_file_path = self.output_file(sys_name, supercell_name)
            self.structure.to(fmt="poscar", filename=output_file_path)