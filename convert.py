import mdtraj as md

frame0 = md.load_frame('trajectory.h5', 0)
frame0.save("traj.pdb")

# comment by Saul

traj = md.load('trajectory.h5')
traj.save("traj.xtc")

