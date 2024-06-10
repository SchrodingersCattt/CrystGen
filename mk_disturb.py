from pymatgen.core import Structure
import os
from mol_identifier import MoleculeIdentifier

class StructurePerturber:
    def __init__(self, input_file, perturb, scale_list, num_replicas, radii_file):
        self.input_file = input_file
        self.structure = Structure.from_file(input_file)
        self.perturb = perturb
        self.scale_list = scale_list
        self.num_replicas = num_replicas
        self.mol_identifier = MoleculeIdentifier(radii_file)

    def perturb_and_scale(self):
        for ii in range(self.num_replicas):
            sys_name = self.input_file.split('.cif')[0]
            self.structure.perturb(self.perturb)
            volume = self.structure.volume

            for scale in self.scale_list:
                scale = float(scale)
                self.structure.scale_lattice(volume * scale)
                directory = f'{sys_name}_disturbed/{sys_name}_{self.perturb}_{scale}_{ii}'

                if not os.path.exists(directory):
                    os.makedirs(directory)

                output = f'{directory}/POSCAR'
                self.structure.to(output)
                print(f'Done making perturbations and output to {output}')