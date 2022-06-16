import sys
sys.path.append("../")

import transformations
print(transformations.__file__)
import numba_transform
import timeit
import numpy as np
import MDAnalysis as mda


rep = 3
num = 5

print("-------------- UnWrapper ------------------------")

u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")

bonds1 = transformations.Unwrapper(sel).bonds
bonds2 = numba_transform.Unwrapper(sel).bonds
print(bonds1.shape,bonds2.shape)
print(np.array_equal(bonds1, bonds2))

u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")
print(timeit.timeit(lambda: transformations.Unwrapper(sel), number =1))
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")
print(timeit.timeit(lambda: numba_transform.Unwrapper(sel), number =1))


#quit()

print("-------------- MolWrapper -----------------------")

u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")

frags1 = transformations.MolWrapper(sel).mols
frags1 = [np.sort(frags1[i]) for i in np.argsort([a[0] for a in frags1])]
frags2 = numba_transform.MolWrapper(sel).mols
frags2 = [np.sort(frags2[i]) for i in np.argsort([a[0] for a in frags2])]
print(len(frags1),len(frags2), [(len(f1),len(f2)) for f1,f2 in zip(frags1,frags2)])
print([np.array_equal(f1,f2) for f1,f2 in zip(frags1,frags2)])

u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")
print("cython    ", timeit.timeit(lambda: transformations.MolWrapper(sel), number =1))
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")
print("MDAnalysis", timeit.timeit(lambda: numba_transform.MolWrapper_mdafrags(sel), number =1))
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("protein")
print("Python    ", timeit.timeit(lambda: numba_transform.MolWrapper(sel), number =1))

#u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
#sel = u.select_atoms("protein")

#print(np.min(timeit.repeat(lambda: transformations.MolWrapper(sel), number =num, repeat=rep))/num)

#u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
#sel = u.select_atoms("protein")


#print(np.min(timeit.repeat(lambda: numba_transform.MolWrapper(sel), number =num, repeat=rep))/num)

print("-------------- Timings -----------------------")

print("Unwrapper whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("all")
print(np.min(timeit.repeat(lambda: transformations.Unwrapper(sel), number =num, repeat=rep))/num)

print("MolWrapper whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("all")
print(np.min(timeit.repeat(lambda: transformations.MolWrapper(sel), number =num, repeat=rep))/num)
