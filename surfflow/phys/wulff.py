import numpy as np
from pymatgen.analysis.wulff import WulffShape


def get_wulff_shape(lattice, miller_list, surfen_list):
    miller_list = np.asarray(miller_list)
    surfen_list = np.asarray(surfen_list)
    # get the indices that would sort the surfen_list in reverse order
    sort_idx = np.argsort(surfen_list)[::-1]
    # sort the surfen_list in reverse order
    surfen_list = surfen_list[sort_idx]
    # sort the miller_list
    miller_list = miller_list[sort_idx]
    miller_list = [tuple(m) for m in miller_list]

    miller_surfen_dict = {m: s for m, s in zip(miller_list, surfen_list)}

    ws = WulffShape(
        lattice=lattice,
        miller_list=list(miller_surfen_dict.keys()),
        e_surf_list=list(miller_surfen_dict.values()),
    )
    return ws
