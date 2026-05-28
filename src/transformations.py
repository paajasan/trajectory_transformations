#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from MDAnalysis.analysis import align
from MDAnalysis.lib.mdamath import triclinic_vectors

import warnings

try:
    import parmed
    pmd_imported = True

except ImportError as e:
    pmd_msg = e
    pmd_imported = False

# Type hints
from typing import Optional
from MDAnalysis.core.groups import AtomGroup, Atom
from MDAnalysis.coordinates.base import Timestep


from . import _ctransformations

"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""


def add_bonds(univ, topfile):
    """
    Add bonds to the universe by making a ParmEd structure from topfile and
    adding all the bonds. No checking is performed between the universe and
    ParmEd structure, so as long as all the indices are found in univ, new
    bonds are made even if they make no sense in reality.

    This can be used to add connections for virtual sites or to add bond info
    into an mdanalysis universe that's missing them.

    Parameters:
        univ:       MDAnalysis universe to add the bonds to
        topfile:    String or path to the topology file that parmed.load_file
                    can read and get bonds from.
    """
    if (not pmd_imported):
        raise ImportError("ParmEd has not be imported. Make sure it is installed\n"
                          "Original message: {pmd_msg}")

    pmd_struct = parmed.load_file(topfile)
    bonds = [(b.atom1.idx, b.atom2.idx) for b in pmd_struct.bonds]

    univ.add_bonds(bonds)


def make_whole(ts: Timestep, bonds: np.ndarray, sel: np.ndarray):
    box = ts.triclinic_dimensions
    ts.positions[sel] = _ctransformations.make_whole(
        ts.positions[sel], bonds, box)

    return ts


class Unwrapper:
    """ 
    Make molecules in ag whole over the pbc. Only considers
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
        starters: Atom group or list of atoms that are guaranteed to stay
                    in box. If multiple atoms are part of same molecule, only
                    first is guaranteed. [default: []]
        virts:    Boolean of whether to search manually for virtual sites and
                    to connect them. Virtual sites are assumed to be any atom that
                    is in a residue with more than one atom, but has no bonds. Bonds
                    are added first to any non-virtual atom in the residue, or if not
                    found, to any virtual atom. The bonds are not added to the original
                    universe, only considered internally. [default: False]
        initsetup: If True, setup is run when initializing object, otherwise only on first
                    execution. [default: False]
    returns:
        transformation function
    """

    def __init__(self, ag: AtomGroup, starters: AtomGroup | list[Atom] = [], virts: bool = False, initsetup: bool = False):
        self.starters = starters
        self.ag = ag
        self.sel = self.ag.indices
        self._virts = virts
        self.__setup_run = False
        if (initsetup):
            self._setup()

    def _setup(self) -> None:
        if (self.__setup_run):
            warnings.warn("Unwrapper setup being run after __setup_run is already. "
                          "Continuing without new setup.", RuntimeWarning)
            return
        bond_ind = self.ag.universe.bonds.to_indices()
        if (self._virts):
            newbonds = _ctransformations.bond_unbonded(
                self.sel,
                bond_ind,
                self.ag.resindices
            )
            bond_ind = np.concatenate((bond_ind, newbonds))

        self.bonds = _ctransformations.traverse_mol(
            self.sel,
            bond_ind,
            np.array([s.index for s in self.starters], dtype=int)
        )
        self.__setup_run = True
        # No need to remember atom group or starters after setup
        del self.ag
        del self._virts
        del self.starters

    def __call__(self, ts: Timestep) -> Timestep:
        if (not self.__setup_run):
            self._setup()
        return make_whole(ts, self.bonds, self.sel)


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
        initsetup: If True, setup is run when initializing object, otherwise only on first
                    execution. [default: False]
    returns:
        transformation function
    """

    def __init__(self, ag: AtomGroup, virts: bool = False, initsetup: bool = False):
        self.ag = ag
        self.selection = ag.indices
        self.weights = ag.masses.astype(ag.positions.dtype)
        self._virts = virts
        self.__setup_run = False
        if (initsetup):
            self._setup()

    def _setup(self) -> None:
        if (self.__setup_run):
            warnings.warn("MolWrapper setup being run after __setup_run is already True."
                          "Continuing without new setup.", RuntimeWarning)
            return

        bond_ind = self.ag.universe.bonds.to_indices()
        if (self._virts):
            newbonds = _ctransformations.bond_unbonded(
                self.selection,
                bond_ind,
                self.ag.resindices
            )
            bond_ind = np.concatenate((bond_ind, newbonds))
        self.mols, self.nmols = _ctransformations.find_frags(
            self.selection,
            bond_ind
        )

        del self.ag
        del self._virts

    def get_frags(self):
        if (not self.__setup_run):
            self._setup()
        mols = []
        for i in range(self.nmols):
            mols.append(self.selection[np.array(self.mols) == i])
        return mols

    def __call__(self, ts: Timestep):
        if (not self.__setup_run):
            self._setup()
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
        centre:        Boolean of whether to centre the selection [default: True]
        superposition: Boolean of whether to centre the selection and fit rotationally [default: True]
        subselection:  The atom group to move and/or rotate, None to use ag. [default: None]
    returns:
        transformation function
    """

    def __init__(self,
                 ag: AtomGroup,
                 centre: bool = True,
                 superposition: bool = True,
                 subselection: Optional[AtomGroup] = None):
        self.seli = ag.indices.copy()
        self.ref = ag.positions.copy()
        self.ref_com = ag.center_of_mass()
        self.w = ag.masses.copy()
        self.totw = self.w.sum()
        if (subselection is None):
            self.subsel = self.seli
        else:
            self.subsel = subselection.indices
        if (superposition):
            self.func = self.superpos
        elif (centre):
            self.func = self.centre
        else:
            self.func = self.nothing

    def __call__(self, ts: Timestep) -> Timestep:
        return self.func(ts)

    def superpos(self, ts: Timestep) -> Timestep:
        sel = ts.positions[self.seli].copy()
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        sel -= sel_com
        R, rmsd = align.rotation_matrix(
            sel, self.ref-self.ref_com, weights=self.w)
        sel = sel @ R.T
        ts.positions[self.subsel] = (
            (ts.positions[self.subsel]-sel_com) @ R.T)+self.ref_com
        return ts

    def centre(self, ts: Timestep) -> Timestep:
        sel = ts.positions[self.seli]
        sel_com = np.sum(self.w*sel.T, axis=-1)/self.totw
        ts.positions[self.subsel] += self.ref_com-sel_com
        return ts

    def nothing(self, ts: Timestep) -> Timestep:
        return ts


class Precenter:
    """
    Center a single atom before putting molecules back in box 

    parameters:
        ag:             Atom group of atoms to move, from the reference universe
        centre_atom:    The Atom to centre. If None the closest atom to box centre
                        is used. [default: None]
        subselection:   The atom group to move, None to use ag. [default: None]

    returns:
        transformation function
    """

    def __init__(self,
                 ag: AtomGroup,
                 centre_atom: Optional[Atom] = None,
                 subselection: Optional[AtomGroup] = None):
        self.seli = ag.indices.copy()
        if (centre_atom is None):
            box = triclinic_vectors(ag.universe.dimensions)
            box_center = (box/2).sum(axis=0)
            dists = np.linalg.norm(ag.positions-box_center, axis=-1)
            self.centre_atom = ag[dists.argmin()].index
        else:
            self.centre_atom = centre_atom.index

        if (subselection is None):
            self.sel = ag.indices
        else:
            self.sel = subselection.indices

        self.diff = None

    def __call__(self, ts: Timestep) -> Timestep:
        box = ts.triclinic_dimensions
        box_center = (box/2).sum(axis=0)
        self.diff = box_center-ts.positions[self.centre_atom]
        ts.positions[self.sel] = ts.positions[self.sel]+self.diff
        return ts

    def uncenter(self, ts: Timestep):
        if (self.diff is None):
            raise ValueError(
                "Precentering has to be called every time before uncentering"
            )

        ts.positions[self.sel] = ts.positions[self.sel]-self.diff
        self.diff = None
        return ts
