import sys
import trajectory_transformations as transformations
sys.path.insert(0,"../")
import transformations_numba as numba_transform
print(transformations.__file__)
print(numba_transform.__file__)
import timeit
import numpy as np
import MDAnalysis as mda
import subprocess as subp
import os
import pstats, cProfile

def notransform(ts):
    return ts

class NoTransform():
    def __init__(self,sel):
        pass
    def __call__(self,ts):
        return ts

def write_transform(sel,transform):
    with mda.Writer("tmp.xtc", sel.n_atoms) as writer:
        for ts in sel.universe.trajectory:
            ts = transform(ts)
            writer.write(sel.atoms)

def write_traj(struct, traj, selstr="all", trajout="tmp.xtc", transform=transformations.Unwrapper):
    u = mda.Universe(struct, traj)
    sel = u.select_atoms(selstr)
    trns = transform(sel)
    with mda.Writer(trajout, sel.n_atoms) as writer:
        for ts in sel.universe.trajectory:
            ts = trns(ts)
            writer.write(sel.atoms)


def write_trjconv(struct, traj, *args, input=b"0\n", trajout="tmp.xtc", **kwargs):
    command = ["gmx", "trjconv", "-f", traj, "-s", struct] + \
              ["-"+a for a in args]
    for k in kwargs:
        command += ["-"+k, kwargs[k]]

    with open("output_trjconv.txt", "w") as fout:
        compProc = subp.run(command, stdout=fout, stderr=subp.STDOUT, input=input)



#cProfile.run('write_traj("../rsc/struct.tpr","../rsc/traj.xtc", transform=transformations.MolWrapper)', "profile.prof")

#s = pstats.Stats("profile.prof")
#s.strip_dirs().sort_stats("time").print_stats("transformations")


#cProfile.run('write_traj("../rsc/struct.tpr","../rsc/traj.xtc")', "profile.prof")


#s = pstats.Stats("profile.prof")
#s.strip_dirs().sort_stats("time").print_stats("transformations")


print("Make whole protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
u.atoms.positions += (50,-50,-50)
transform_cyt = transformations.Unwrapper(sel)
transform_num = numba_transform.Unwrapper(sel)
p_orig = u.atoms.positions.copy()
transform_cyt(u.trajectory[0])
cyt_pos = u.atoms.positions.copy()
u.atoms.positions = p_orig
transform_num(u.trajectory[0])
num_pos = u.atoms.positions
print(np.all((num_pos-cyt_pos)<1e-4))
print(np.max(np.linalg.norm(num_pos-cyt_pos, axis=-1)))

rep = 1
num = 1
write_trjconv("../rsc/struct.tpr","../rsc/traj.xtc", "nobackup",pbc="whole",ur="tric")
print()
print(np.min(timeit.repeat(lambda: write_traj("../rsc/struct.tpr","../rsc/traj.xtc", transform=NoTransform), number =num, repeat=rep))/num)
print(np.min(timeit.repeat(lambda: write_traj("../rsc/struct.tpr","../rsc/traj.xtc"), number =num, repeat=rep))/num)
print(np.min(timeit.repeat(lambda: write_trjconv("../rsc/struct.tpr","../rsc/traj.xtc", "nobackup",pbc="whole",ur="tric"), number =num, repeat=rep))/num)


print("\nIteration\n")

print("notransform protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
print(np.min(timeit.repeat(lambda: write_transform(sel, notransform), number =num, repeat=rep))/num)

print("cython protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
transform = transformations.Unwrapper(sel)
print(np.min(timeit.repeat(lambda: write_transform(sel, transform), number =num, repeat=rep))/num)

print("Numba protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
transform = numba_transform.Unwrapper(sel)
print(np.min(timeit.repeat(lambda:  write_transform(sel, transform), number =num, repeat=rep))/num)


print("notransform whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("all")
print(np.min(timeit.repeat(lambda: write_transform(sel, notransform), number =num, repeat=rep))/num)

print("cython whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("all")
transform = transformations.Unwrapper(sel)
print(np.min(timeit.repeat(lambda: write_transform(sel, transform), number =num, repeat=rep))/num)


#quit()

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

frags1 = transformations.MolWrapper(sel).get_frags()
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

#quit()

rep = 5
num = 10
print("\nSetup\n")

print("Unwrapper whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("all")
print(np.min(timeit.repeat(lambda: transformations.Unwrapper(sel), number =num, repeat=rep))/num)

print("MolWrapper whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/struct.gro")
sel = u.select_atoms("all")
print(np.min(timeit.repeat(lambda: transformations.MolWrapper(sel), number =num, repeat=rep))/num)




rep = 5
num = 10
print("\nIteration\n")

print("notransform protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
print(np.min(timeit.repeat(lambda: write_transform(sel, transform), number =num, repeat=rep))/num)

print("cython protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
transform = transformations.Unwrapper(sel)
print(np.min(timeit.repeat(lambda: write_transform(sel, transform), number =num, repeat=rep))/num)

print("Numba protein")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("protein")
transform = numba_transform.Unwrapper(sel)
print(np.min(timeit.repeat(lambda:  write_transform(sel, transform), number =num, repeat=rep))/num)


print("notransform whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("all")
print(np.min(timeit.repeat(lambda: write_transform(sel, notransform), number =num, repeat=rep))/num)

print("cython whole")
u = mda.Universe("../rsc/struct.tpr","../rsc/traj.xtc")
sel = u.select_atoms("all")
transform = transformations.Unwrapper(sel)
print(np.min(timeit.repeat(lambda: write_transform(sel, transform), number =num, repeat=rep))/num)

print("\nStart to finish")
print("cython whole")
print(np.min(timeit.repeat(lambda: write_traj("../rsc/struct.tpr","../rsc/traj.xtc"), number =num, repeat=rep))/num)

print("trjconv whole")
print(np.min(timeit.repeat(lambda: write_trjconv("../rsc/struct.tpr","../rsc/traj.xtc", "nobackup",pbc="whole",ur="tric"), number =num, repeat=rep))/num)

print("cython prot")
print(np.min(timeit.repeat(lambda: write_traj("../rsc/struct.tpr","../rsc/traj.xtc", "protein"), number =num, repeat=rep))/num)

print("trjconv prot")
print(np.min(timeit.repeat(lambda: write_trjconv("../rsc/struct.tpr","../rsc/traj.xtc", "nobackup", input=b"Protein\n",pbc="whole",ur="tric"), number =num, repeat=rep))/num)
