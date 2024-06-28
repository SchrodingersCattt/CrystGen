from pymatgen.core import Structure
from multiprocessing import Pool
import json
import os
import time

class MoleculeIdentifier:
    def __init__(self, tol_factor=1.1):
        _script_path = os.path.dirname(os.path.abspath(__file__))
        radii_file = os.path.join(_script_path, "BASIC_DATA", "atomic_radius.json")
        self.pm_to_angstrom = 0.01
        self.radii = self.load_radii(radii_file)
        self.tol_factor = tol_factor

    def load_radii(self, filename):
        with open(filename, 'r') as f:
            radii = json.load(f)
        return radii

    def get_radii(self):
        return self.radii

    def calc_distance(self, struct, ii, jj):
        return struct.get_distance(ii, jj)

    def calc_bond(self, struct, i, j):
        dist = self.calc_distance(struct, i, j)
        elem_i = struct.species[i].name
        elem_j = struct.species[j].name
        cov_radius_i = self.radii[elem_i] * self.pm_to_angstrom
        cov_radius_j = self.radii[elem_j] * self.pm_to_angstrom
        bond_length = cov_radius_i + cov_radius_j
        return dist < bond_length * self.tol_factor

    def find_molecules(self, struct):
        ss = time.time()
        molecules_site_indices = []
        visited = set()

        for i in range(len(struct)):
            if i not in visited:
                molecule = [i]
                visited.add(i)
                stack = [i]

                while stack:
                    j = stack.pop()
                    for k in range(len(struct)):
                        if k not in visited and self.calc_bond(struct, j, k):
                            molecule.append(k)
                            visited.add(k)
                            stack.append(k)

                molecules_site_indices.append(molecule)
        ee = time.time()
        print("find mol", ee - ss)
        return molecules_site_indices