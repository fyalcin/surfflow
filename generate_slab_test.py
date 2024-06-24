from mp_api.client import MPRester

from htflow_utils.shaper import Shaper
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# mpid = "mp-2657"

with MPRester() as mpr:
    data = mpr.summary.search(material_ids="mp-2657")

structure = data[0].structure

structure_conv = SpacegroupAnalyzer(structure).get_conventional_standard_structure()
structure_conv.to("POSCAR_conv_bulk.vasp", "poscar")

sg_params = {'miller': (1, 1, 0),
             'max_index': None,
             'symmetrize': False,
             'slab_thick': 8,
             'vac_thick': 30.0,
             'primitive': True,
             'lll_reduce': True,
             'tol': 0.05,
             'max_normal_search': 'max',
             'resize': True,
             'preserve_terminations': True,
             'min_thick_A': 10.0,
             'minimize_structures': False,
             'match_ouc_lattice': True,
             'calculate_bonds': True,
             'center_slab': True,
             'max_nsites': 100,
             'filter_polar': True}

slabs_dict, _ = Shaper.generate_slabs(bulk_conv=structure_conv, sg_params=sg_params)

for miller, slabs in slabs_dict.items():
    for index, slab in enumerate(slabs):
        millerstr = "".join([str(i) for i in miller])
        slab.to(f"POSCAR_{millerstr}_{index}.vasp", "poscar")