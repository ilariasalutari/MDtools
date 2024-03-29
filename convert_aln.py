import mdtraj as md
import argparse
import os

parser = argparse.ArgumentParser(
                    prog='Convert_traj',
                    description='Convert the trajectory from .h5 to .xtc',
                    epilog='Thank you !')

parser.add_argument('-f', '--traj_file', help='The input trajectory', required=True)
parser.add_argument('-top', '--topology_file', help = 'Optional (depending on the traj format): a reference topology (.pdb or .gro format)', required=False)
parser.add_argument('-o', '--out_file', help='Add a name for output', default='traj.xtc')
parser.add_argument('-skip', '--stride', help = 'Save every n-th frame', type = int,  default = 1)
parser.add_argument('-sel', '--selection', help = 'Select what to visualize using the MDtraj synthax between quotes (default: no water or ions)', default = 'not (resname HOH or name NA or name CL)')


args = parser.parse_args()
dic_input = vars(args)

## Setting the path and filename variables:

out_dir = os.path.dirname(dic_input['out_file'])

xtc_name = os.path.basename(dic_input['out_file'])
pdb_name = xtc_name[:-4] + ".pdb"
out_dir_pdb = os.path.join(out_dir, pdb_name)

## Convert and save:

if dic_input['traj_file'].endswith(".h5"):

    frame0 = md.load_frame(dic_input['traj_file'], index=0)

    protein_ind = frame0.topology.select(dic_input['selection'])
    frame0 = frame0.atom_slice(protein_ind)
    frame0.save(out_dir_pdb)

    traj = md.load(dic_input['traj_file'], atom_indices = protein_ind, stride = dic_input['stride'])

    prot_CA = traj.topology.select("protein and name CA")
    traj.superpose(frame0, frame=0, atom_indices=prot_CA, ref_atom_indices=prot_CA)

    traj.save(dic_input['out_file'])

else:

    frame0 = md.load_frame(dic_input['traj_file'], top = dic_input['topology_file'], index=0)

    protein_ind = frame0.topology.select(dic_input['selection'])
    frame0 = frame0.atom_slice(protein_ind)
    frame0.save(out_dir_pdb)

    traj = md.load(dic_input['traj_file'], top = dic_input['topology_file'], atom_indices = protein_ind, stride = dic_input['stride'])

    prot_CA = traj.topology.select("protein and name CA")
    traj.superpose(frame0, frame=0, atom_indices=prot_CA, ref_atom_indices=prot_CA)

    traj.save(dic_input['out_file'])

