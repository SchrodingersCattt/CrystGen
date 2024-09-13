import numpy as np
from pymatgen.core import Structure
import json
import os
import time


class MoleculeIdentifier:
    def __init__(self, tol_factor=1.1):
        _script_path = os.path.dirname(os.path.abspath(__file__))
        radii_file = os.path.join(_script_path, "../BASIC_DATA", "atomic_radius.json")
        self.pm_to_angstrom = 0.01
        self.radii = self.load_radii(radii_file)
        self.tol_factor = tol_factor

    def load_radii(self, filename):
        with open(filename, 'r') as f:
            radii = json.load(f)
        return radii

    def get_radii(self):
        return self.radii

    def calc_distance_matrix(self, struct):
        
        ss = time.time()
        # Get all fractional coordinates in the structure
        #frac_coords = np.array([site.frac_coords for site in struct])

        # Calculate and return the pairwise distance matrix considering PBC
        dist_matrix = np.array(
            [[struct.get_distance(i, j) for j in range(len(struct))] for i in range(len(struct))]
        )
        ee = time.time()
        print("find mollll", ee - ss)
        return dist_matrix

    def calc_bond_matrix(self, struct):
        # Calculate the distance matrix
        dist_matrix = self.calc_distance_matrix(struct)
        # Get the radii for each element in the structure
        cov_radii = np.array([self.radii[elem.name] * self.pm_to_angstrom for elem in struct.species])

        # Calculate bond length matrix
        bond_length_matrix = self.tol_factor * (cov_radii[:, None] + cov_radii[None, :])

        # Determine where bonds are formed
        bond_matrix = dist_matrix < bond_length_matrix
        return bond_matrix

    def find_molecules(self, struct):

        # Initialize an empty list to store molecular site indices
        molecules_site_indices = []

        # Initialize an all-false array indicating visited sites
        visited = np.zeros(len(struct), dtype=bool)

        # Get the bond matrix
        bond_matrix = self.calc_bond_matrix(struct)

        # Identify molecules using bond matrix
        for i in range(len(struct)):
            if not visited[i]:
                # List to accumulate the atoms belonging to the same molecule
                molecule = []

                # Use a queue to explore connected molecules
                queue = [i]
                visited[i] = True

                while queue:
                    current = queue.pop(0)
                    molecule.append(current)

                    # Find all unvisited neighbors
                    neighbors = np.where(bond_matrix[current] & ~visited)[0]
                    queue.extend(neighbors)

                    # Mark those neighbors as visited
                    visited[neighbors] = True

                molecules_site_indices.append(molecule)

        return molecules_site_indices