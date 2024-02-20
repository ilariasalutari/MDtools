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
import matplotlib.colors as col


parser = argparse.ArgumentParser(
                    prog = 'Calculate and plot DSSP',
                    description = 'Calculate and plot DSSP annotation in a trajectory',
                    epilog = 'Thank you !')

parser.add_argument('-f', '--traj_file', help = 'The input trajectory', required=True)
parser.add_argument('-top', '--topology_file', help = 'The first frame or reference topology (.pdb or .gro format)', default='traj.pdb', required = True)
parser.add_argument('-o', '--out_file', help = 'Add a name for image output', default = 'distances.png')
parser.add_argument('-sel', '--res_sel', help = 'Optional: select stretch of residues for plotting')
parser.add_argument('-ref', '--xtal_ref', help = 'Optional: reference structure to plot the crystal residues. The residue sequence needs to match the topology pdb.')
parser.add_argument('-skip', '--stride', help = 'Optional: use every n-th frame. Recommended at least 10 for long simulations.', type = int, default = 1)
parser.add_argument('-frame_freq', '--freq', help = 'Optional: saving frequency in picoseconds (to have a nice plot)', type = int)
parser.add_argument('-sim_time', '--length', help = 'Optional: length of the simulation in nanoseconds (to have a nice plot)', type = int)
parser.add_argument('-sys', '--system', help = 'Optional: select parts of the system (default: no water or ions)', default = 'protein and not (resname HOH or name NA CL)')

args = parser.parse_args()
### keep all the args of the program in a dictionary:
dic_input = vars(args)


## Function to compute distances:
def get_dssp(traj, ref):

    dssp_assigned = md.compute_dssp(traj, simplified=True)
    ## H helix (alpha, 3-helix or 5-helix); E extended strands or residue in beta bridge; C coil (H-bonded turn, bend or loop); other options if simplified = False (see Mdtraj) }
    legend_simple = {"H":1, "E":2, "C":3, "NA":4}

    legend_list = [legend_simple[element] for element in dssp_assigned.flatten()]
    dssp_mapped = np.array(legend_list).reshape(dssp_assigned.shape)
     
    table_traj, bonds = traj.atom_slice(traj.topology.select("protein")).topology.to_dataframe()
    un_resnumbs = (sorted(table_traj['resSeq'].unique()))
      
    ## Dictionary to map the traj residue numbers to a numbering starting from 0:
    map_indices = dict()
    for c, i in enumerate(un_resnumbs):
        map_indices[c] = i
    
    resnames = []
    for i in range(len(un_resnumbs)):
        resnames.append(traj.topology.residue(i))

    ## Map_indices from 0 is always used, unless another -ref is provided.
    ## Not a problem if the top or ref structures have gaps, because everything is mapped to 0-numbering and the corresponding names are used only in the y-labels.
    
    ## Switch to crystal residue labels:
    if ref != None:
        ## IMP: the sequence of residues needs to be identical between the two structures, so need to exclude before the ACE and NME?
        xtal = md.load(dic_input['xtal_ref'])
        table_xtal, bonds = xtal.atom_slice(xtal.topology.select("protein")).topology.to_dataframe()

        ## Map the xtal residue numbers to a numbering starting from 0:        
        map_indices = dict()
        for c, i in enumerate(sorted(table_xtal['resSeq'].unique())): 
            map_indices[c] = i

        ## Save crystal residue numbers and names:
        resnumb = []
        for  a in [xtal.topology.residue(index_protein).resSeq for index_protein in map_indices.keys()]:
            resnumb.append(a)

        resnames = list(zip([str(traj.topology.residue(index_protein).name) for index_protein in map_indices.keys()], [str(n) for n in resnumb] ))   
        names_numb = ["-".join(n) for n in resnames] 

    else:
        names_numb = [str(prot_lig.topology.residue(index_protein)) for index_protein in range(0, len(dssp_mapped.T))]


    ## output:
    return map_indices, names_numb, dssp_mapped       


## Load the traj:
pdb =  md.load(dic_input['topology_file'], 0)
traj = md.load(dic_input['traj_file'], top=pdb, stride=dic_input['stride'])
prot_lig = traj.atom_slice(traj.topology.select(dic_input['system']))

## Use the above function:
mapped, res_names, dssp_values = get_dssp(prot_lig, dic_input['xtal_ref'])

##  Use pandas dataframe for a nice plot:
dssp_df = pd.DataFrame(dssp_values)

## Choose:
colors = {"#c1737c":1, "#73a3c1":2, "#a3c173":3}
cMap = col.ListedColormap(colors)


## Plot the DSSP vs. time (ns)
time_ticks = 5

fig = plt.figure(figsize = (12, 7))
ax = sns.heatmap(dssp_df.T, cmap=cMap)

## To add residues on y-axis:

if type(dic_input['res_sel']) == str:  
    skip_resi = 1

    limit_resi = [int(s) for s in dic_input['res_sel'].split(", ")]
    print(limit_resi)

    ## Not necessary, if selected limits are consitent with the topology file provided:
    ## map the limits provided to the 0-numbering of before:
    # invert_map = {b: a for (a, b) in mapped.items()}
    # convert_limits = sorted([invert_map[i] for i in limit_resi])

    resi_stride = np.arange(min(limit_resi), max(limit_resi)+1, skip_resi, dtype=int)
    resi_plot = [res_names[index_prot] for index_prot in resi_stride]
    print(resi_plot)

    ax.set_ylim(min(resi_stride), max(resi_stride))      
    ### all major tick labels have to be shifted by one half up to be more readable
    ax.set_yticks(np.arange(min(resi_stride), max(resi_stride)+1, skip_resi))      
    ax.set_yticks(np.arange(min(resi_stride), max(resi_stride)+1, skip_resi)-.5, minor=True) 
    ax.set_yticklabels(resi_plot, rotation=0, size=10,  va='bottom')
    ax.tick_params(which="major", bottom=True, left=False)
    ax.tick_params(which="minor", bottom=True, left=True)
    ax.set_ylabel("Residue", size=16)

elif dic_input['res_sel'] == None:

    skip_resi = 10
    ### If no reference numbering is provided, use the trajectory numering:
    resi_stride = np.arange(min(mapped.keys()), max(mapped.keys())+1, skip_resi, dtype=int)  
    resi_plot = [res_names[index_prot] for index_prot in resi_stride]

    ax.set_ylim(min(resi_stride), max(resi_stride))  
    ### all major tick labels have to be shifted by one half up to be more readable
    ax.set_yticks(np.arange(min(resi_stride), max(resi_stride)+1, skip_resi))      
    ax.set_yticks(np.arange(min(resi_stride), max(resi_stride)+1, skip_resi)-.5, minor=True) 
    ax.set_yticklabels(resi_plot, rotation=0, size=10,  va='bottom')
    ax.tick_params(which="major", bottom=True, left=False)
    ax.tick_params(which="minor", bottom=True, left=True)
    ax.set_ylabel("Residue", size=16)
    
else:
    raise ValueError("Not a valid residue selection for the y-axis.")


## Calculate the time for the plot:
if dic_input['freq'] != None:
    time_ns = (dssp_df.index * dic_input['freq'] * dic_input['stride'])/1000  
    tot_sim = max(time_ns)

elif dic_input['freq'] == None and dic_input['length'] != None:
    frame_freq = (dic_input['length']*1000) / (len(dssp_df) * dic_input['stride'])   ## saving frequency in ps
    tot_sim = dic_input['length']    # time

else:
    raise ValueError("\n", 'WARNING: to have Time (ns) on the x-axis, provide either the saving frequency (in ps) or the total simulation length (in ns).', "\n")
    #time_ns = np.nan
    #tot_sim = np.nan

## Uses the time (ns):
ax.set_xticks(np.linspace(start=min(dssp_df.index), stop=max(dssp_df.index), num=10))
seq_xlab = np.linspace(start=0, stop=tot_sim, num=10)
ax.set_xticklabels([int(i) for i in seq_xlab], rotation=0, size=10)
ax.set_xlabel("Time (ns)", size=16)


# Manually specify colorbar labelling after it's been generated
colorbar = ax.collections[0].colorbar
ori_ticks = colorbar.get_ticks()
new_ticks = np.linspace(start=int(min(ori_ticks))+.35, stop=int(max(ori_ticks))-.35, num=3)
colorbar.set_ticks(new_ticks)
colorbar.set_ticklabels(['Helix', 'Beta-bridge', 'Turn or Loop'], rotation=90, va='center', size=16) 

## save the plot:
fig = ax.get_figure()
fig.savefig(dic_input['out_file'], dpi=300)

