"""
 Geometry utilities
"""

import numpy as np
from ase import Atoms
from typing import Union


def basis_projected_pos_sorted(atom: Atoms, insbottom: bool = False) -> tuple:
    """Caclculate basis-vector-projected atom positions and sort.

    Args:
        atom: An atom image.
        insbottom:  Whether to insert perodic image of the top atom to the bottom. True for yes, False for no.

    Returns:
        A ``tuple`` containing
        
         - isorted_ppos (numpy.ndarray): The atom indices of sorted projected positions.
         - sorted_ppos (numpy.ndarray): The array of the sorted projected positions.
         - sorted_pgap (numpy.ndarray): The array of the sorted projected gaps beween layers (0th:1st-0th, 1st:2nd-1st, ...).
    """
    cell = atom.cell
    vol = cell.volume
    if vol == 0:  # When no cell info
        ppos = atom.get_positions()
        cpar = [1e10] * 3
    else:
        cpar = cell.cellpar()
        mcell = [[cpar[0], 0, 0], [0, cpar[1], 0], [0, 0, cpar[2]]]
        rpos = atom.get_scaled_positions(wrap=False)
        ppos = np.matmul(rpos, mcell)  # projected position

    isorted_ppos = ppos.argsort(axis=0)
    sorted_ppos = np.take_along_axis(ppos, isorted_ppos, axis=0)
    ret_ppos = sorted_ppos.copy()
    if insbottom:
        ret_ppos = np.insert(sorted_ppos, 0, [sorted_ppos[-1] - cpar[:3]], axis=0)
    ret_ppos = np.append(ret_ppos, [sorted_ppos[0] + cpar[:3]], axis=0)
    sorted_pgap = ret_ppos[1:] - ret_ppos[:-1]

    return isorted_ppos, ret_ppos, sorted_pgap


def identify_layers(atom: Atoms, sidx: int, mingap: float = 1.9) -> tuple:
    """Indentify layers based on the projected length along the axis ``sidx``.

    The layers are collection of atoms whose positions are differed less than the ``mingap``. This doesn't recognize layers lying on periodic boundaries. Please use ``--wslab`` to wrap a separated slab.

    Args:
        atoms:  An atom image
        sidx:   A slab axis index. 1-3(a-c).

    Returns:
        A ``tuple`` containing
        
         - layeredge_hi (numpy.ndarray): The array of sorted projected position of the highest atom in a layer.
         - layeredge_lo (numpy.ndarray): The array of sorted projected position of the lowest atom in a layer.
         - ilayers  (numpy.ndarray): The array of indices of atoms in a layer.

    """
    isorted_ppos, sorted_ppos, sorted_pgap = basis_projected_pos_sorted(atom)
    isorted_sppos = isorted_ppos.transpose()[sidx - 1]
    sorted_spgap = sorted_pgap.transpose()[sidx - 1]

    isorted_hi = np.where(sorted_spgap > mingap)[0]
    isorted_lo = np.array([0, *(isorted_hi[:-1] + 1)])

    ilayers = []
    for i, _ in enumerate(isorted_lo):
        ilayers.append(
            [isorted_sppos[x] for x in range(isorted_lo[i], isorted_hi[i] + 1)]
        )

    layeredge_hi = sorted_ppos[isorted_hi]
    layeredge_lo = sorted_ppos[isorted_lo]
    return layeredge_hi, layeredge_lo, ilayers


def _axis_to_int(ax):
    if type(ax) == str:
        return ord(ax) - ord("a") + 1
    return ax


def roll_axes(atom: Atoms, src: Union[str, int], tgt: Union[str, int]) -> Atoms:
    """Roll axes of the ``atom`` image from the ``src`` axis to the ``tgt`` axis.

    Example:
        Rolling 'a' axis to 'c' axis results the following axis changes: a->c, b->a, c->b.

        >>> roll_axes(atom, 'a', 'c')

    Args:
        atom: An atom image.
        src: The source axis. '{a-c}' or 1-3.
        tgt: The tatget axis. '{a-c}' or 1-3.

    Returns:
        The axes rolled atom image.

    """

    isrc = _axis_to_int(src)
    itgt = _axis_to_int(tgt)
    if isrc == itgt:
        return atom

    syms = atom.get_chemical_symbols()

    shift = (itgt + 3 - isrc) % 3
    rcell = np.roll(atom.cell, shift, axis=0)
    rpos = np.roll(atom.get_scaled_positions(), shift, axis=1)  # should be scaled
    ratom = Atoms(symbols=syms, scaled_positions=rpos, cell=rcell, pbc=True)

    return ratom


def align_slabaxistoz(atom: Atoms, sidx: int) -> Atoms:
    """Align the slab axis (``sidx``) of the ``atom`` to z axis.

    Args:
        atom: An atom image.
        sidx: A slab axis to be aligned. 1-3.

    Returns:
        The Z-alined atom image.

    """
    from scipy.spatial.transform import Rotation as R

    ratom = atom
    if sidx == 1 or sidx == 2:
        print(f" - Slab axis({chr(ord('a')+sidx-1)}) is not c axis, rolling ...")
        if len(atom.constraints) > 0:
            print("  | constraints are not preservered...")
        ratom = roll_axes(atom, sidx, 3)
    print(" - aligning c->[001] and a->[100]")
    if len(ratom.constraints) > 0:
        print("  | constraints is not preservered...")

    cell = ratom.cell
    ralign = R.align_vectors(
        [[0, 0, 1], [1, 0, 0]], [cell[2], cell[0]], weights=[1, 0.1]
    )
    # ralign=R.align_vectors([[0,1,0],[1,0,0]],[cell[2],cell[0]],weights=[1,0.1])
    rr = ralign[0]
    rrcell = rr.apply(cell)
    rrpos = rr.apply(ratom.get_positions())
    rratom = Atoms(
        symbols=ratom.get_chemical_symbols(), positions=rrpos, cell=rrcell, pbc=True
    )
    return rratom


def identify_slabaxis(atom: Atoms, isslabwrap: bool = False, mingap: float = 12):
    """Identify slab axis of the ``atom`` image.

    If a gap of the projected positions in an axis is greater than ``mingap``, the axis is identifed as the slab axis. The search order is c,b,a axes, and the first indentified axis is returned.

    Args:
        atom: An atom image
        isslabwrap: Whether to wrap the separated slab. True for yes, False for no.
        mingap: A minimum position gap to be a slab.

    Returns:
        0 for non-slab, otherwise 1-3 (a-c axes).

    """
    cell = atom.cell
    vol = cell.volume
    if vol == 0:
        return 3

    slabidx = 0
    cpar = cell.cellpar()

    isorted_ppos, sorted_ppos, sorted_pgap = basis_projected_pos_sorted(atom)
    imaxgap = sorted_pgap.argmax(axis=0)
    maxgap = sorted_pgap.max(axis=0)
    natom = atom.get_global_number_of_atoms()
    # print(imaxgap)
    for i in range(2, -1, -1):
        # TODO mingap to user input
        if maxgap[i] >= mingap:
            slabidx = i + 1
            if isslabwrap and imaxgap[i] < natom - 1:
                slen = cpar[i]
                shift = slen - sorted_ppos[imaxgap[i] + 1][i] + 5
                print(
                    f" - Separated slab, shifting {shift:.3f} A in {chr(ord('a')+i)} axis ..."
                )
                sar = [0] * 3
                sar[i] = shift / slen  # Projected Fractional
                sar = np.matmul(cell, sar)  # Cartesian
                atom.translate(sar)
                atom.wrap()
            break
            
    return slabidx
