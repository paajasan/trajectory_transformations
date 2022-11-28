# Trajectory transformations #


## What? ###

[On-the-fly](https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/) transformations for MDAnalysis trajectories.

### Numba vs cython ###

The project has two version of the same tools. The numba-accelerated one is a bit slower in iteration and a lot slower in setup than the cython accelerated one, but does not need to be built. If you plan on transforming whole systems (inluding all waters), I seriously suggest using the cython version. With small proteins without hydrogens the differences less noticeable.

### unwrap ####

Makes "fragments" whole over the PBC. Fragments are groups of bonded atoms in the selection. If the molecules are continuous in the selection, these are the same as molecules.

### wrap ####

Put fragment COM back into the box.

### center and/or superposition ####

Centres the selection COM to reference COM and optionally does a rotational superposition onto the reference.

## Why? ###

The built-in unwrap method was super slow. With one system the built-in method took ~8 minutes, while our method took 2.6 s to transform the same frames. This was only 0.6 s overhead compared to iterating the frames without transformation.

**Note**: there is a [pull request](https://github.com/MDAnalysis/mdanalysis/pull/3169#issue-831915405) for MDAnaysis, that should make its transformations faster by many orders of magnitude. This project will most likely still be faster for large systems, but of course by a smaller margin than currently.

## How? ###


### unwrap ####

As a setup does a Depth-First-Search (DFS) on the graph of bonded atoms to gather bond information. Then for each frame makes these bonds unbroken over the PBC, in the order they came up in the search.

### wrap ####

As before, uses DFS to find "fragments" (groups of bonded atoms). If the selection is continuous, these are the same as molecules.

### superpositioning ####

Centering by simply translating, so that COM is in (0,0,0), after it supersitioning with `MDAnalysis.analysis.align.rotation_matrix` to get optimal rotation matrix. Finally translates by refrerence COM.




## Requirements ###

1. Python 3 with NumPy
1. [MDAnalysis](https://docs.mdanalysis.org/stable/index.html)
1. [Numba](https://numba.pydata.org/) (only for the numba-accelerated case)

If you have conda, you can make sure all dependencies are met by running

```
conda install -c conda-forge numpy numba mdanalysis
```


## Installation ##

### Simple (numba accelerated) ###

Simply, copy the `transformations_numba.py` from the project root to  the current directory (or anywhere in the PYTHONPATH). 

The transformations can then be imported with
```
from transformations_numba import Unwrapper, Wrapper, Superpos
```

### Slightly less simple (cython accelerated) ###

This way the project will be installed as a package, so make sure you have the correct environment activated. **Optional**: If you use conda, you can make a new environment with the requirements installed with

```
conda create env myTransformEnv -c conda-forge numpy mdanalysis
conda activate myTransformEnv
```

The actual command to build and install the project is simply
```
pip install .
```
which of course has to be run in the project folder. This will take care of any buld-time dependencies (like cython) even if you do not have them installed. If you want to update the tool, you can simply rerun the command with the new version. To uninstall the tool run

```
pip uninstall trajectory_transformations
```


The transformations can then be imported with
```
from trajectory_transformations import Unwrapper, Wrapper, Superpos
```


## Usage ##


Let's set up an example case, like in the [MDAnalysis documentation](https://userguide.mdanalysis.org/stable/trajectories/transformations.html), assuming you have imported the transformations as shown above

```
import MDAnalysis as mda
from MDAnalysis.tests.datafiles import TPR, XTC

u = mda.Universe(TPR, XTC)
protein = u.select_atoms('protein')
```

The transformations are callable classes, taking at least the atom group as parameter in the constructor. 

```
unwrap = Unwrapper(protein)
```

Then we add it to the trajectory with 

```
u.trajectory.add_transformations(unwrap)
```

Now you can iterate over the trajectory just as you would normally


```
for ts in u.trajectory:
    # Do something with transformed positions
```

**Optionally**, it can be called on each timestep of the trajectory when iterating it as

```
for ts in u.trajectory:
    # Do something with positions before transformation
    ts = unwrap(ts)
    # Do something with transformed positions
```

In that case the positions can also be accessed before the transformation.

All three transformations work with the same principle.


### Virtual site handling ###

These methods have been tested on a simple TIP4P system. They should in theory work with different systems, but further testing is still needed.

As MDAnalysis does not add the connections of virtual sites to bonds, virtual sites are not connected to the rest of the molecule. Two methods are included to work around this.

#### Guessing connections within residues ####

If the virtual sites are all within residues and MDAnalysis correctly understands labels them, you can use the virts-option of the Unwrapper and Wrapper initialisers. It will guess virtual sites as atoms that are part of a residue with other atoms while having no bonds. Then it will connect these to any (effectively the first) non-virtual site within the same residue (if the residue is only made up of virtual sites, it will take one of them). As long as the virtual sites connect to the real atoms this way and any single residue is not larger than half the box length, this should work.

The intialisation then simply becomes:

```
unwrap = Unwrapper(protein, virts=True)
```

The connections are only considered internally, and are not added to the universe object.

#### Reading them straight from a topology ####

The more robust option is to read the data from a topology file with [ParmEd](https://parmed.github.io/ParmEd/html/index.html). It is only an optional dependency and can be installed any time before or after this package. If it is not installed, and ImportError will be raised upon calling the function.

This method should work with any system that ParmEd can handle, but it does add some overhead (in the order of a few seconds to tens of seconds depending on your system size) from reading the system in. If you work on many systems successively, this can end up being quite a lot of slower. You also need to have the topology file handy, though you can always make one from teh tpr file if you have gromacs, with `gmx dump`.

You can import the function as

```
from trajectory_transformations import add_bonds
```

The initialisation then works (assuming you have the topology in `"topol.top"`) as

```
add_bonds(u, "topol.top")
unwrap = Unwrapper(protein, virts=True)
```