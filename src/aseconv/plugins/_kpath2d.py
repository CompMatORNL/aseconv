"""
    2D Kpath generating functions. Modified from a following source.
    https://github.com/aiidateam/aiida-core/blob/b7e37ce50ba1d6a5de62ad95357befdedda0de65/aiida/tools/data/array/kpoints/legacy.py
"""
from __future__ import annotations
from ase import Atoms
from aseconv import geoutil as gu

def kpath_2d(args: argparse.Namespace, oatom: Atoms, lnums: list, ksidx: int) -> tuple(Atoms,dict):
    """2D Kpath generator. 
    
    Args:
        args:  Parsed arguments.
        oatom: Input atom image.
        lnums: List of the number of atoms per each element.
        ksidx: Slab index.
        
    Returns:
        ratoms: Rotated atom.
        kp: Dict of kpoints. 
    
    """
    import spglib as spg
    import numpy as np

    arrs = {
        1: {  # square
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "M": (0.5, 0.5, 0),
                "X": (0.5, 0, 0),
                "Y": (0, 0.5, 0),
            },
            "path": [
                ("M", "GAMMA"),
                ("GAMMA", "X"),
                ("X", "M"),
                ("M", "Y"),
                ("Y", "GAMMA"),
            ],
        },
        2: {  # (primitive) rectangular
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "X": (0.5, 0, 0),
                "Y": (0, 0.5, 0),
                "S": (0.5, 0.5, 0),
            },
            "path": [
                ("S", "X"),
                ("X", "GAMMA"),
                ("GAMMA", "S"),
                ("S", "Y"),
                ("Y", "GAMMA"),
            ],
        },
        3: {  # centered rectangular (rhombic)
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "X": (0.5, 0.5, 0),
                "Y1": (0.25, 0.75, 0),
                "Y": (-0.25, 0.25, 0),
                "C": (0, 0.5, 0),
            },
            "path": [
                ("Y1", "X"),
                ("X", "GAMMA"),
                ("GAMMA", "Y"),
                ("Y", "C"),
            ],
        },
        4: {  # hexagonal
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "M": (0.5, 0, 0),
                "K": (1 / 3, 1 / 3, 0),
            },  # , "M'": (0,0.5,0)
            "path": [("K", "M"), ("M", "GAMMA"), ("GAMMA", "K")],
        },
        -4: {  # hexagonal P1
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "M": (0.5, 0, 0),
                "K": (1 / 3, 1 / 3, 0),
                "M'": (0, 0.5, 0),
            },
            "path": [
                ("K", "M"),
                ("M", "GAMMA"),
                ("GAMMA", "K"),
                ("K", "M'"),
                ("M'", "GAMMA"),
            ],
        },
        6: {  # hexagonal2
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "M": (0.5, 0, 0),
                "K": (2 / 3, 1 / 3, 0),
            },
            "path": [
                ("M", "K"),
                ("K", "GAMMA"),
                ("GAMMA", "M"),
            ],
        },
        -6: {  # hexagonal2 P1
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "M": (0.5, 0, 0),
                "K": (2 / 3, 1 / 3, 0),
                "M'": (0, 0.5, 0),
                "K'": (1 / 3, 2 / 3, 0),
            },
            "path": [
                ("GAMMA", "M"),
                ("M", "K"),
                ("K", "GAMMA"),
                ("GAMMA", "M'"),
                ("M'", "K'"),
                ("K'", "GAMMA"),
            ],
        },
        5: {  # oblique acute
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "X": (0.5, 0, 0),
                "Y": (0, 0.5, 0),
                "A": (0.5, 0.5, 0),
            },
            "path": [
                ("X", "GAMMA"),
                ("GAMMA", "Y"),
                ("A", "GAMMA"),
            ],
        },
        -5: {  # oblique obtuse
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "X": (0.5, 0, 0),
                "Y": (0, 0.5, 0),
                "A": (0.5, -0.5, 0),
            },
            "path": [
                ("A", "GAMMA"),
                ("GAMMA", "X"),
                ("GAMMA", "Y"),
            ],
        },
        7: {  # rhombus acute
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "X": (0.5, 0, 0),
                "Y": (0, 0.5, 0),
                "A": (0.5, 0.5, 0),
                "C": (0, 0, 0),
            },
            "path": [
                ("C", "X"),
                ("X", "GAMMA"),
                ("GAMMA", "C"),
                ("C", "A"),
                ("A", "GAMMA"),
            ],
        },
        -7: {  # rhombus obtuse
            "point_coords": {
                "GAMMA": (0, 0, 0),
                "X": (0.5, 0, 0),
                "Y": (0, 0.5, 0),
                "A": (0.5, -0.5, 0),
                "C": (0, 0, 0),
            },
            "path": [
                ("C", "X"),
                ("X", "GAMMA"),
                ("GAMMA", "C"),
                ("C", "A"),
                ("A", "GAMMA"),
            ],
        },
    }

    ratom = gu.align_slabaxistoz(oatom, ksidx)  # Always must be rotated for 2D path

    spcell = (ratom.cell, ratom.get_scaled_positions(), lnums)
    data = spg.get_symmetry_dataset(spcell, symprec=0.5)
    if data is None:
        raise Exception(" - [kpath-2d] Could not get symmetry info..")

    grp = data["international"]
    grpn = data["number"]
    brav = _kpidx_bravais(ratom, grpn)
    bidx = brav["index"]
    angs = ratom.cell.angles()

    print(" - [kpath-2d] %s(%d) (%.2f degree)" % (brav["extended_name"], bidx, angs[2]))
    kp = arrs[bidx].copy()
    if bidx == 7 or bidx == -7:
        sign = np.sign(bidx)
        sind = lambda degrees: np.sin(np.deg2rad(degrees))
        cosd = lambda degrees: np.cos(np.deg2rad(degrees))

        alpha = (180 - angs[2]) / 2  # angle between A and b1
        k = 1 / (2 * (sign * cosd(180 - 2 * alpha) + 1))
        cvec = (1 - k, sign * k, 0)
        kp["point_coords"].update({"C": cvec})

    kp["primitive_types"] = spcell[0]
    kp.update(
        {
            "primitive_lattice": spcell[0],
            "primitive_positions": spcell[1],
            "primitive_types": spcell[2],
            "spacegroup_international": grp,
            "spacegroup_number": data["number"],
        }
    )

    """
	cr=kp['point_coords']
	cr2d={}
	islab=2 # rotated 
	for k,v in cr.items():
		if (v[islab] != 0): continue
		vc=list(v)
		del vc[islab]
		cr2d.update({k:vc})
	print(cr)
	print(cr2d)
	kp['point_2dcoords']=cr2d
	"""

    return ratom, kp



def _kpidx_bravais(atom: Atoms, grpn: int) -> dict:
    """Obtain KPath index from Bravis info.
    
    Args:
        atom: Atom image.
        grpn: Symmetry group number from spg library.
        
    Returns:
        2d kpath index.
        
    """
    epsilon_length = 0.1
    epsilon_angle = 0.1

    def l_are_equals(a, b):
        return abs(a - b) <= epsilon_length

    def a_are_equals(a, b):
        return abs(a - b) <= epsilon_angle

    out_of_plane_index = 2  # rotated to z axis
    in_plane_indexes = list(set(range(3)) - set([out_of_plane_index]))
    cell = atom.cell
    vectors = [cell[i] for i in in_plane_indexes]
    lens = [cell.lengths()[i] for i in in_plane_indexes]
    angs = cell.angles()
    phi = angs[out_of_plane_index]
    equa_len = l_are_equals(lens[0], lens[1])
    equa_90 = a_are_equals(phi, 90)
    if equa_90 and equa_len:
        ret = {"index": 1, "short_name": "sq", "extended_name": "square"}
    elif equa_90:
        ret = {"index": 2, "short_name": "rec", "extended_name": "rectangular"}
    elif equa_len and a_are_equals(phi, 120):
        ret = {
            "index": 4,
            "short_name": "hex",
            "extended_name": "hexagonal",
        }
        if grpn == 1 or grpn == 2:
            ret = {
                "index": -4,
                "short_name": "hexP1",
                "extended_name": "hexagonal P1",
            }
    elif equa_len and a_are_equals(phi, 60):
        ret = {"index": 6, "short_name": "hex2", "extended_name": "hexagonal2"}
        if grpn == 1 or grpn == 2:
            ret = {"index": -6, "short_name": "hex2P1", "extended_name": "hexagonal2P1"}
    elif equa_len and phi > 90:
        ret = {"index": -7, "short_name": "rhoo", "extended_name": "rhombus optuse"}
    elif equa_len and phi < 90:
        ret = {"index": 7, "short_name": "rhoa", "extended_name": "rhombus acute"}
    elif phi > 90:
        ret = {
            "index": -5,
            "short_name": "oblo",
            "extended_name": "oblique optuse",
        }
    else:
        ret = {
            "index": 5,
            "short_name": "obla",
            "extended_name": "oblique acute",
        }
    return ret

