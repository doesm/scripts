import numpy as np
import MDAnalysis as mda

topol =  # archivo .tpr
traj =  # archivo .xtc

u = mda.Universe(topol, traj)
f = open(arch, "x")
f.write("@     yaxis label \"Count C4B Beads \"\n")
f.write("@     xaxis label \"Time ( ns ) \"\n")
f.write("@     title \"C4B between bilayers \"\n")

for ts in u.trajectory:
    mem_center = u.select_atoms(
        "resname POPC or resname POPS or resname POP2").center_of_geometry()
    c4b_selection = "name C4B and prop abs z <= " + \
        str(mem_center[2]+10) + " and prop abs z >= " + str(mem_center[2]-10)
    linea = str(u.trajectory.time/1000) + " " + \
        str(len(u.select_atoms(c4b_selection))) + "\n"
    f.write(linea)

f.close()
