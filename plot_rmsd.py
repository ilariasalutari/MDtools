import os
import warnings
import argparse
import math
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib 
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so


parser = argparse.ArgumentParser(
                    prog = 'Calculate and plot RMSD in trajectory',
                    description = 'Calculate and plot RMSD of protein and ligand. NOTE: the ligand must be called LIG.',
                    epilog = 'Thank you !')

parser.add_argument('-f', '--traj_file', help = 'The input trajectory', required=True)
parser.add_argument('-top', '--topology_file', help = 'The first frame or reference topology (.pdb or .gro format)', default='traj.pdb', required = True)
parser.add_argument('-o', '--out_file', help = 'Add a name for image output', default = 'rmsd.png')
parser.add_argument('-aln', '--aln_atoms', help = 'Optional: select atoms of the protein for RMSD alignment', default='protein and name CA')
parser.add_argument('-skip', '--stride', help = 'Optional: use every n-th frame. Recommended at least 10 for long simulations.', type = int, default = 1)
parser.add_argument('-frame_freq', '--freq', help = 'Optional: saving frequency in picoseconds (to have a nice plot)', type = int)
parser.add_argument('-sim_time', '--length', help = 'Optional: length of the simulation in nanoseconds (to have a nice plot)', type = int)
parser.add_argument('-sel', '--system', help = 'Optional: select parts of the system (default: no water or ions)', default = 'not (resname HOH or name NA or name CL)')

args = parser.parse_args()
### keep all the args of the program in a dictionary:
dic_input = vars(args)


## Function to compute RMSD:
def get_rmsd(traj, ref, indices):

    if type(indices) == str:
        ind_aln = traj.topology.select(indices)
    
    elif type(indices) in [list, np.ndarray, np.array, tuple]:   
        ind_aln  = indices
    
    else:
        raise ValueError("Not a valid atom index or string for protein RMSD selection")

    ## output:
    return md.rmsd(traj, ref, frame=0, atom_indices=ind_aln)*10      ### RMSD in Angstrom


## Load the traj:
pdb =  md.load(dic_input['topology_file'], 0)
traj = md.load(dic_input['traj_file'], top=pdb, stride=dic_input['stride'])
prot_lig = traj.atom_slice(traj.topology.select(dic_input['system']))

#prot_ind = prot_lig.topology.select("protein")
lig_ind = prot_lig.topology.select("resname LIG")

## Calculate the RMSDs:
rmsd_prot = get_rmsd(prot_lig, pdb, dic_input['aln_atoms'])   
rmsd_lig = get_rmsd(prot_lig, pdb, lig_ind)

##  Use pandas dataframe for a nice plot:
rmsd_df = pd.concat([pd.DataFrame(rmsd_prot), pd.DataFrame(rmsd_lig)], axis=1)
rmsd_df.columns = ["Protein", "Ligand"]           ## will be used in the plot legend
rmsd_df['frame_ind'] = [ind for ind in rmsd_df.index]


## Calculate the time for the plot:
if dic_input['freq'] != None:
    rmsd_df['time_ns'] = (rmsd_df.index * dic_input['freq'] * dic_input['stride'])/1000  

elif dic_input['freq'] == None and dic_input['length'] != None :
    frame_freq = (dic_input['length']*1000) / (len(rmsd_df) * dic_input['stride'])
    rmsd_df['time_ns'] = (rmsd_df.index * frame_freq * dic_input['stride'])/1000 

else:
    print("\n", 'WARNING: to have Time (ns) on the x-axis, provide either the saving frequency (in ps) or the total simulation length (in ns).', "\n")
    rmsd_df['time_ns'] = np.nan


rmsd_to_plot = pd.melt(rmsd_df, ['time_ns'])
rmsd_to_plot.columns = ['time_ns', 'RMSD of:', 'RMSD']  


## Make the timeseries plot:

if pd.notnull(rmsd_df['time_ns']).all() == True:      ## true only if all elements of time_ns exist (are not NaN)
    
    rmsd_df = rmsd_df.drop('frame_ind', axis=1)
    rmsd_to_plot = pd.melt(rmsd_df, ['time_ns'])
    rmsd_to_plot.columns = ['time_ns', 'RMSD of:', 'RMSD'] 

    sns.set_style("whitegrid")      
    p = sns.lineplot(data=rmsd_to_plot, x="time_ns", y="RMSD", hue="RMSD of:", linewidth=1)
    p.set_xlim(0, max(rmsd_to_plot['time_ns']))
    plt.xlabel("Time (ns)", fontsize=12)

else:
    
    rmsd_df = rmsd_df.drop('time_ns', axis=1)
    rmsd_to_plot = pd.melt(rmsd_df, ['frame_ind'])
    rmsd_to_plot.columns = ['frame_ind', 'RMSD of:', 'RMSD'] 
    
    sns.set_style("whitegrid")   
    p = sns.lineplot(data=rmsd_to_plot, x="frame_ind", y="RMSD", hue="RMSD of:", linewidth=1)
    p.set_xlim(0, max(rmsd_to_plot['frame_ind']))
    plt.xlabel("Frame index", fontsize=12)


## other plot parameters:
p.set_ylim(0, math.ceil(max(rmsd_to_plot['RMSD'])+2))
plt.ylabel("RMSD (Å)", fontsize=12)
plt.tick_params(axis='x', labelsize=12)
plt.tick_params(axis='y', labelsize=12)
p.set_title(("Alignment to "+ dic_input['aln_atoms']), fontsize=12)
sns.move_legend(p, "upper right")
plt.setp(p.get_legend().get_texts(), fontsize='10') # for legend text
plt.setp(p.get_legend().get_title(), fontsize='11')

## save the plot:
fig = p.get_figure()
fig.savefig(dic_input['out_file'], dpi=300) 

