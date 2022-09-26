import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import distances

topol =  # .tpr
traj =  # .xtc
u = mda.Universe(topol, traj)
print("Is this the simulation with protein: 0 = no / 1 = yes")
protein = input()
print("Which SMA do you want to mesure?")
n_sma = input()
f = open(arch, "x")
f.write("@  yaxis label \"distance\"\n")
f.write("@  xaxis label \"time (ns)\"\n")
f.write("@  title \"SMA-Fusion Pore distance\"\n")

if int(protein) == 0:
    first = 24993 + (int(n_sma) - 1) * 195
    t = 24993 + int(n_sma) * 195 - 1
else:
    first = 25104 + (int(n_sma) - 1) * 195
    t = 25104 + int(n_sma) * 195 - 1

sma = u.atoms[first:t]
for ts in u.trajectory:
    mem_center = u.select_atoms(
        "resname POPC o r resname POPS or resname POP2 ").center_of_geometry()
    c4b_selection = "name C4B and prop abs z <= " + \
        str(mem_center[2] + 10) + " and prop abs z >= " + \
        str(mem_center[2] - 10)
    if len(u.select_atoms(c4b_selection)) > 10:
        cog_pore = u.select_atoms(c4b_selection).center_of_geometry()
        dist = distances.distance_array(
            cog_pore, sma.positions, box=u.dimensions)
        linea = str(u.trajectory.time/1000) + " " + str(np.min(dist)) + "\n"
        f.write(linea)

f.close()
