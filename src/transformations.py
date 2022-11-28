#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from MDAnalysis.analysis import align

from . import _ctransformations

"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""


def make_whole(ts, bonds, sel):
    box = ts.triclinic_dimensions
    ts.positions[sel] = _ctransformations.make_whole(
        ts.positions[sel], bonds, box)

    return ts


class Unwrapper:
    """ Make molecules in ag whole over the pbc. Only considers
        continuously bonded molecules, so if molecules are in many parts,
        the nonbonded parts will be ignored and can be broken.
            Starters is a list (or atom group) of atoms to use as starting points.
        They will be put into the box if they are outside and each consecutive
        bonded atom will be moved by a box vector if it is more than half the
        length of the box.
            Virts can help with virtual sites, as MDAnalysis does not return their
        connections as bonds. A new bond is made to a random atom in the residue,
        so as long as virtual sites are within residues with real atoms and any 
        single residue is not larger than half the box vector in any direction,
        this should work just fine. If it doesn't, you can manually add bonds to
        each virtual site before calling this initialiser, e.g. with the add_virt_bonds
        function. 
        parameters:
            ag:       Atom group of molecules to make whole
            starters: Atom goup (or list) of atoms that are guaranteed to stay
                      in box. If multiple atoms are part of same molecule, only
                      first is guaranteed.
            virts:    Boolean of whether to search manually for virtual sites and
                      to connect them. Virtual sites are assumed to be any atom that
                      is in a residue with more than one atom, but has no bonds. Bonds
                      are added first to any non-virtual atom in the residue, or if not
                      found, to any virtual atom. The bonds are not added to the original
                      universe, only considered internally.
        returns:
            transformation function
    """

    def __init__(self, ag, starters=[], virts=False):
        self.sel = ag.indices
        bond_ind = ag.universe.bonds.to_indices()
        if (virts):
            newbonds = _ctransformations.bond_unbonded(
                self.sel,
                bond_ind,
                ag.resindices
            )
            bond_ind = np.concatenate((bond_ind, newbonds))

        self.bonds = _ctransformations.traverse_mol(
            self.sel,
            bond_ind,
            np.array([s.index for s in starters], dtype=int)
        )
        self.sel = ag.indices

    def __call__(self, ts):
        box = ts.triclinic_dimensions
        ts.positions[self.sel] = _ctransformations.make_whole(
            ts.positions[self.sel], self.bonds, box)
        return ts


class MolWrapper:
    """
    Put centre of mass of molecules in selection to box
    parameters:
        ag:       Atom group of molecules to put in box
        virts:    Boolean of whether to search manually for virtual sites and
                  to connect them. Virtual sites are assumed to be any atom that
                  is in a residue with more than one atom, but has no bonds. Bonds
                  are added first to any non-virtual atom in the residue, or if not
                  found, to any virtual atom. This does not work for virtual sites
                  that are between residues.
    returns:
        transformation function
    """

    def __init__(self, ag, virts=False):
        self.selection = ag.indices
        self.weights = ag.masses.astype(ag.positions.dtype)
        bond_ind = ag.universe.bonds.to_indices()
        if (virts):
            newbonds = _ctransformations.bond_unbonded(
                self.sel,
                bond_ind,
                ag.resindices
            )
            bond_ind = np.concatenate((bond_ind, newbonds))
        self.mols, self.nmols = _ctransformations.find_frags(
            self.selection,
            bond_ind
        )

        assert not np.any(np.array(self.mols) < 0), "FFuuuu, %d" % np.sum(
            np.array(self.mols) < 0)
        assert not np.any(np.array(self.mols) >= self.nmols), "FFuuuu2, %d, %d" % (
            np.sum(np.array(self.mols) >= self.nmols, self.nmols))

    def get_frags(self):
        mols = []
        for i in range(self.nmols):
            mols.append(self.selection[np.array(self.mols) == i])
        return mols

    def __call__(self, ts):
        box = ts.triclinic_dimensions
        ts.positions[self.selection] = _ctransformations.wrap_mols(
            ts.positions[self.selection],
            self.weights,
            self.mols,
            self.nmols,
            box
        )
        return ts


class Superpos:
    """
    Superposition for optimal mass weighted rmsd
    parameters:
        ag:            Atom group of atoms to fit, from the reference universe
        centre:        Boolean of whether to centre the selection
        superposition: Boolean of whether to centre the selection and fit rotationally
    returns:
        transformation function
    """

    def __init__(self, ag, centre, superposition):
        self.seli = ag.indices.copy()
        self.ref = ag.positions.copy()
        self.ref_com = ag.center_of_mass()
        self.w = ag.masses.copy()
        self.totw = self.w.sum()
        if (superposition):
            self.func = self.superpos
        elif (centre):
            self.func = self.centre
        else:
            self.func = self.nothing

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def superpos(self, ts):
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        sel -= sel_com
        R, rmsd = align.rotation_matrix(
            sel, self.ref-self.ref_com, weights=self.w)
        sel = sel @ R.T
        ts.positions[self.seli] = sel+self.ref_com
        return ts

    def centre(self, ts):
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        ts.positions[self.seli] += self.ref_com-sel_com
        return ts

    def nothing(self, ts):
        return ts
