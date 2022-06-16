import cython
import numpy as np
cimport numpy as np
import warnings
from libcpp.unordered_set cimport unordered_set as cunset
from libcpp.unordered_map cimport unordered_map as cunmap
from libcpp.map cimport map as cmap
from libcpp.stack cimport stack as cstack
from libcpp.vector cimport vector as cvector
from cython.operator cimport dereference

"""
This module includes functions to make on the fly transformations for MDAnalysis trajectories
"""

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def find_frags(np.intp_t[:] ag, int[:,:] bond_ind):
    """ Traverse all molecules that atoms in ag are part of.
        Parameters:
            ag:       Indices of atoms to check
            bonds:    int(n,2) array of atom indices for n bonds.

        Returns:
            List of indices for each molecule, only including indices in ag.
    """

    cdef cvector[cvector[int]]   mols
    cdef cunset[int]             done,agi
    cdef cunmap[int,cunset[int]] bonds
    cdef cvector[int]            mol
    cdef int                     i,a,start
    cdef cstack[int]             stack

    # A dict to hold the bond info
    for i in range(bond_ind.shape[0]):
        bonds[bond_ind[i,0]].insert(bond_ind[i,1])
        bonds[bond_ind[i,1]].insert(bond_ind[i,0])

    # A set of the selection
    for i in range(ag.shape[0]):
        agi.insert(ag[i])


    # We iterate until there are no more atoms to search
    while(not agi.empty()):
        # A set to hold info on whether we have "visited" atom
        done.clear()
        # Take any atom as the starting point for next search
        start = dereference(agi.begin())
        agi.erase(start)

        # Put start in stack
        stack.push(start)
        mol.push_back(start)
        done.insert(start)

        # Continue until stack is empty
        while(not stack.empty()):
            # Pop from top of stack
            curr = stack.top()
            stack.pop()
            # Iterate atoms bonded to popped atom
            for a in bonds[curr]:
                # Ignore if atom visited
                if(done.count(a)):
                    continue

                # Add atom to list if it is in selection
                if(agi.count(a)):
                    mol.push_back(a)
                    agi.erase(a)
                # Add atom on top of stack and mark it as visited
                stack.push(a)
                done.insert(a)


        mols.push_back(mol)
        mol.clear()

    return mols


@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
def traverse_mol(np.intp_t[:] ag, int[:,:] bond_ind, np.intp_t[:] starters):
    """ Traverse all graphs of bonded atoms in ag, using atoms in starters (or the
        smallest indexes) as starting points for a DFS.
        Parameters:
            ag:       Indices of atoms to check
            bonds:    int(n,2) array of atom indices for n bonds.
            starters: Indices of atoms to use as starting points.
        Returns:
            Indexes for each bonds of the bonded atoms within the selection
            (i.e. indexes within selection, not within whole universe). int[n,2].
    """
    cdef list                     bonds
    cdef cunmap[int,int]          agi_indexes
    cdef cunset[int]              done
    cdef cunmap[int,cvector[int]] bond_map
    cdef int                      i, a, next_starter=0,next_first=0
    cdef cstack[int]              stack

    # A dict to hold the bond info
    for i in range(bond_ind.shape[0]):
        bond_map[bond_ind[i,0]].push_back(bond_ind[i,1])
        bond_map[bond_ind[i,1]].push_back(bond_ind[i,0])
        
    # A dict to map the indices between system and selection
    for i in range(ag.shape[0]):
        agi_indexes[ag[i]] = i

    # A list to hold the bonds
    bonds = []

    # We iterate until there are no more atoms to search
    while(not agi_indexes.empty()):
        # A set to hold info on whether we have "visited" atom
        done.clear()
        # Decide next starting atom
        start = -1
        # Go thorugh starter atoms until done or we choose one
        while(start < 0 and next_starter < starters.shape[0]):
            a = starters[next_starter]
            next_starter += 1
            if(agi_indexes.count(a)):
                start = a

        # If we didn't choose one take the next from ag still in agi_indexes
        while(start < 0):
            a = ag[next_first]
            next_first += 1
            if(agi_indexes.count(a)):
                start = a

        # Make a stack to hold atoms to visit and add start atom
        stack.push(start)
        done.insert(start)
        # Continue until stack is empty
        while(not stack.empty()):
            # Pop from top of stack
            curr = stack.top()
            stack.pop()
            # Iterate atoms bonded to popped atom
            for a in bond_map[curr]:
                # Ignore if atom visited or if it is not in selection
                if(done.count(a) or not agi_indexes.count(a)):
                    continue

                # Add bond to list
                bonds.append((agi_indexes[curr], agi_indexes[a]))
                # Add atom on top of stack and mark it as visited
                stack.push(a)
                done.insert(a)

        # Remove visited atoms from selction to keep track of whether we are done
        for a in done:
            agi_indexes.erase(a)

    if(not bonds):
        warnings.warn("No bonds found in unwrap selection")
        return np.zeros((0,2),dtype=int)
    return np.array(bonds)



