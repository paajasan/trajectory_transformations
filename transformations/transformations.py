#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from MDAnalysis.analysis import align

from . import _ctransformations

"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""


def make_whole(ts, bonds, sel):
    """ Makes molecules in sel whole using bond info from bonds.
        Coordinates are first switched to reciprocal space, then passed to
        _iterate_bonds to actually make them whole. The fixed positions are
        switched back to normal space and set as the current positions for sel.
    """
    box = ts.triclinic_dimensions
    # Inverse box
    invbox = np.linalg.inv(box.astype(np.float64))
    # Transfer coordinates to unit cell and put in box
    unitpos = (ts.positions[sel].astype(np.float64) @ invbox) % 1
    # Fix box in reciprocal space
    unitpos = _ctransformations.iterate_bonds(unitpos, bonds)

    ts.positions[sel] = unitpos @ box

    return ts


class Unwrapper:
    """ Make molecules in ag whole over the pbc. only considers
        continuosly bonded molecules, so if molecules are in many parts,
        the nonbonded parts will be ignored and can be broken.
            Starters is a list (or atom group) of atoms to use as starting points.
        They will be put into the box if they are outside and each consecutive
        bonded atom will be moved by a box vector if it is more than half the
        length of the box.
        parameters:
            ag:       Atom group of molecules to make whole
            starters: Atom goup (or list) of atoms that are guaranteed to stay
                      in box. If multiple atoms are part of same molecule, only
                      first is guaranteed.
        returns:
            transformation function
    """
    def __init__(self, ag, starters=[]):
        self.sel = ag.indices
        self.bonds = _ctransformations.traverse_mol(
                            self.sel,
                            ag.bonds.to_indices(),
                            np.array([s.index for s in starters],dtype=int)
                        )
        self.sel = ag.indices

    def __call__(self, ts):
        return make_whole(ts, self.bonds, self.sel)




class MolWrapper:
    """
    Put centre of mass of molecules in selection to box
    parameters:
        ag:       Atom group of molecules to put in box
    returns:
        transformation function
    """
    def __init__(self,ag):
        self.selection = ag.indices
        self.mols = []
        self.weights = []
        for frag in _ctransformations.find_frags(ag.indices, ag.bonds.to_indices()):
            self.mols.append(frag)
            self.weights.append(ag.universe.atoms[frag].masses)

        self.totw = np.array([np.sum(w) for w in self.weights])

    def __call__(self,ts):
        box = ts.triclinic_dimensions
        # Inverse box
        invbox = np.linalg.inv(box)
        for i,m in enumerate(self.mols):
            pos = np.sum(self.weights[i]*ts.positions[m].T,axis=-1)/self.totw[i]
            # Transfer coordinates to unit cell and put in box
            unitpos = (pos @ invbox) % 1
            newpos = unitpos @ box
            trans = newpos-pos
            ts.positions[m] += trans
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
    def __init__(self,ag, centre, superposition):
        self.seli = ag.indices.copy()
        self.ref = ag.positions.copy()
        self.ref_com = ag.center_of_mass()
        self.w   = ag.masses.copy()
        self.totw = self.w.sum()
        if(superposition):
            self.func = self.superpos
        elif(centre):
            self.func = self.centre
        else:
            self.func = self.nothing

    def __call__(self, *args, **kwargs):
        return self.func(*args, **kwargs)

    def superpos(self,ts):
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        sel -= sel_com
        R, rmsd = align.rotation_matrix(sel, self.ref-self.ref_com, weights=self.w)
        sel = sel @ R.T
        ts.positions[self.seli] = sel+self.ref_com
        return ts

    def centre(self,ts):
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        ts.positions[self.seli] += self.ref_com-sel_com
        return ts

    def nothing(self,ts):
        return ts
