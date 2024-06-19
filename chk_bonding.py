from pymatgen.core.structure import Structure
from tqdm import tqdm
import multiprocessing as mp

class ChkBonding:
    def __init__(self, mol_identifier, intra_crit_dist, intra_crit_lo_ratio, intra_crit_hi_ratio):
        self.mol_identifier = mol_identifier
        self.intra_crit_dist = intra_crit_dist
        self.intra_crit_lo_ratio = intra_crit_lo_ratio
        self.intra_crit_hi_ratio = intra_crit_hi_ratio

    def pass_checking(self, struct, molecules_site_indices):
        overall_set = set([num for sublist in molecules_site_indices for num in sublist])
        check_results = []

        pool = mp.Pool(processes=16)  # Create a pool of worker processes
        results = [pool.apply_async(self._process_molecule, args=(struct, mol_site_idx, overall_set)) for mol_site_idx in molecules_site_indices]
        pool.close()
        pool.join()

        for result in results:
            check_results.extend(result.get())  # Collect the results from worker processes

        if True in check_results:
            return False
        else:
            return True

    def _process_molecule(self, struct, mol_site_idx, overall_set):
        check_results = []
        sel_mol_set = set(mol_site_idx)
        exclude_set = overall_set.difference(sel_mol_set)
        exclude_idx = list(exclude_set)
        sel_mol_idx = list(sel_mol_set)

        for sel_idx in tqdm(sel_mol_idx):
            for exl_idx in exclude_idx:
                chk_res_inter = self.mol_identifier.calc_bond(struct, exl_idx, sel_idx)
                check_results.append(chk_res_inter)

            for other_sel_idx in sel_mol_idx[sel_mol_idx.index(sel_idx) + 1:]:
                intra_crit_1 = struct.get_distance(sel_idx, other_sel_idx) < self.intra_crit_dist
                intra_crit_2 = struct.get_distance(sel_idx, other_sel_idx) < \
                               self.intra_crit_lo_ratio * struct.get_distance(sel_idx, other_sel_idx)
                intra_crit_3 = struct.get_distance(sel_idx, other_sel_idx) > \
                               self.intra_crit_hi_ratio * struct.get_distance(sel_idx, other_sel_idx)
                chk_res_intra = bool(intra_crit_1 or intra_crit_2 or intra_crit_3)
                check_results.append(chk_res_intra)

        return check_results