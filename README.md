# Trajectory transformations #


## What? ###

Trajectory [on-the-fly](https://www.mdanalysis.org/2020/03/09/on-the-fly-transformations/) transformations for MDAnalysis

### unwrap ####

Makes "fragments" whole over the PBC. Fragments are groups of bonded atoms in the selection. If the molecules are continuous in the selection, these are the same as molecules.

### wrap ####

Put fragment COM back into the box.

### center and/or superposition ####

Centres the selection COM to (0,0,0) and optionally does a rotational superposition onto the reference.

## Why? ###

The built-in unwrap method was super slow.

## How? ###


### unwrap ####

Basically as a setup does a DFS on the graph of bonded atoms to save bond information. Then for each frame makes these bonds unbroken over the PBC.

### wrap ####

Uses same search as before to find "fragments": groups of bonded atoms. If The selection is continuous, these are the same as molecules.

## Usage ##
TODO