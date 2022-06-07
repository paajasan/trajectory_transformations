# Trajectory transformations #


## What? ###

Trajectory [on-the-fly](https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/) transformations for MDAnalysis

### `unwrap` ####

Makes "fragments" whole over the PBC. Fragments are groups of bonded atoms in the selection. If the molecules are continuous in the selection, these are the same as molecules.

### `wrap` ####

Put fragment COM back into the box.

### center and/or superposition ####

Centres the selection COM to reference COM and optionally does a rotational superposition onto the reference.

## Why? ###

The built-in unwrap method was super slow. With one system the built-in method took ~8 minutes, while our method took 2.6 s to do transform the same frames. This was only 0.6 s overhead compared to iterating the frames without transformation.

## How? ###


### unwrap ####

Basically as a setup does a DFS on the graph of bonded atoms to save bond information. Then for each frame makes these bonds unbroken over the PBC.

### wrap ####

Uses same search as before to find "fragments": groups of bonded atoms. If The selection is continuous, these are the same as molecules.

### superpositioning ####

Centering by simply translating, so that COM is in (0,0,0), after it supersitioning with `MDAnalysis.analysis.align.rotation_matrix` to get optimal rotation matrix. Finally translates by refrerence COM.

## Usage ##

First, move the `transformations.py` somewhere in the PYTHONPATH, for example in the current directory. Then just import the transformations with 

```
from transformations import Unwrapper, Wrapper, Superpos
```

Then, let's set up an example case, like in the [MDAnalysis documentation](https://userguide.mdanalysis.org/stable/trajectories/transformations.html) 

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




## Requirements ###

1. Python 3 with NumPy
1. [MDAnalysis](https://docs.mdanalysis.org/stable/index.html)
1. [Numba](https://numba.pydata.org/)

If you have conda, you can make sure all dependencies are met by running

```
conda install -c conda-forge numpy numba mdanalysis
```
