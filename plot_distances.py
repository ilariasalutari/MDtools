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
                    prog = 'Calculate distances',
                    description = 'Measure and plot specific distances in a trajectory',
                    epilog = 'Thank you !')

parser.add_argument('-f', '--traj_file', help = 'The input trajectory', required=True)
parser.add_argument('-top', '--topology_file', help = 'The first frame or reference topology (.pdb or .gro format)', default='traj.pdb', required = True)
parser.add_argument('-o', '--out_file', help = 'Add a name for image output', default = 'distances.png')
parser.add_argument('-g1', '--atoms_1', help = 'First group of atom(s) with MDtraj synthax, individual selections separated by comma', required = True)
parser.add_argument('-g2', '--atoms_2', help = 'Second group of atom(s) with MDtraj synthax, individual selections separated by comma', required = True)
parser.add_argument('-pair', '--pair_by_order', help = 'yes or no: calculate pairs of distances according to input order. Default: no (all distance combinations)', default = 'no')
parser.add_argument('-skip', '--stride', help = 'Optional: use every n-th frame. Recommended at least 10 for long simulations.', type = int, default = 1)
parser.add_argument('-frame_freq', '--freq', help = 'Optional: saving frequency in picoseconds (to have a nice plot)', type = int)
parser.add_argument('-sim_time', '--length', help = 'Optional: length of the simulation in nanoseconds (to have a nice plot)', type = int)
parser.add_argument('-sel', '--system', help = 'Optional: select parts of the system (default: no water or ions)', default = 'not (resname HOH or name NA or name CL)')

args = parser.parse_args()
### keep all the args of the program in a dictionary:
dic_input = vars(args)


## Function to compute distances:
def get_distance(traj, sel1, sel2, paired):

    if type(sel1) == str:
        ind1 = np.concatenate([traj.topology.select(s) for s in sel1.split(", ")])      ## keeps the order of the atoms as in input command
    
    elif type(sel1) == list or type(sel1) == np.array or type(sel1) == tuple:
        ind1 = sel1
    
    else:
        raise ValueError("Not a valid atom index or string")

    if type(sel2) == str:
        ind2 = np.concatenate([traj.topology.select(s) for s in sel2.split(", ")])
    
    elif type(sel2) == list or type(sel2) == np.array or type(sel2) == tuple:
        ind2 = sel2
    
    else:
        raise ValueError("Not a valid atom index or string")


    ## Switch to only ordered pairs = same order as in input (with 4 atoms, 2 output distances)
    if paired == "yes":
        dist_all = list(zip(ind1, ind2))

        labels_1 = [prot_lig.topology.atom(i) for i in ind1]
        labels_2 = [prot_lig.topology.atom(j) for j in ind2]
        labels_all = list(zip(labels_1, labels_2))

    else:
        ### use all the possible combinations of atom pairs (with 4 atoms, 4 output distances)
        dist_all = []
        labels_all = []
        for i in ind1:
            for j in ind2:
                dist_all.append([i, j])

                labels_1 = prot_lig.topology.atom(i)
                labels_2 = prot_lig.topology.atom(j)
                labels_all.append([labels_1, labels_2])


    ## output:
    return dist_all, labels_all, md.compute_distances(traj, dist_all)*10       ### distances in Angstrom


## Load the traj:
pdb =  md.load(dic_input['topology_file'], 0)
traj = md.load(dic_input['traj_file'], top=pdb, stride=dic_input['stride'])
prot_lig = traj.atom_slice(traj.topology.select(dic_input['system']))

## Calculate the distances:
atoms, names, dist_values = get_distance(prot_lig, dic_input['atoms_1'], dic_input['atoms_2'], dic_input['pair_by_order'])


##  Use pandas dataframe for a nice plot:
distances = pd.DataFrame(dist_values)  
distances.columns = [str(c) for c in names]                ## will be used in the plot legend
distances['frame_ind'] = [ind for ind in distances.index]


## Calculate the time for the plot:
if dic_input['freq'] != None:
    distances['time_ns'] = (distances.index * dic_input['freq'] * dic_input['stride'])/1000  

elif dic_input['freq'] == None and dic_input['length'] != None :
    frame_freq = (dic_input['length']*1000) / (len(distances) * dic_input['stride'])
    distances['time_ns'] = (distances.index * frame_freq * dic_input['stride'])/1000 

else:
    print("\n", 'WARNING: to have Time (ns) on the x-axis, provide either the saving frequency (in ps) or the total simulation length (in ns).', "\n")
    distances['time_ns'] = np.nan


distances_to_plot = pd.melt(distances, ['time_ns'])
distances_to_plot.columns = ['time_ns', 'Atoms', 'Distance']  


## Make the timeseries plot:

if pd.notnull(distances['time_ns']).all() == True:      ## true only if all elements of time_ns exist (are not NaN)
    
    distances = distances.drop('frame_ind', axis=1)
    distances_to_plot = pd.melt(distances, ['time_ns'])
    distances_to_plot.columns = ['time_ns', 'Atoms', 'Distance'] 

    sns.set_style("whitegrid")      
    p = sns.lineplot(data=distances_to_plot, x="time_ns", y="Distance", hue="Atoms", linewidth=1)
    p.set_xlim(0, max(distances_to_plot['time_ns']))
    plt.xlabel("Time (ns)", fontsize=12)

else:
    
    distances = distances.drop('time_ns', axis=1)
    distances_to_plot = pd.melt(distances, ['frame_ind'])
    distances_to_plot.columns = ['frame_ind', 'Atoms', 'Distance'] 
    
    sns.set_style("whitegrid")   
    p = sns.lineplot(data=distances_to_plot, x="frame_ind", y="Distance", hue="Atoms", linewidth=1)
    p.set_xlim(0, max(distances_to_plot['frame_ind']))
    plt.xlabel("Frame index", fontsize=12)


## other plot parameters:
p.set_ylim(min(0, math.floor(min(distances_to_plot['Distance'])-1)), math.ceil(max(distances_to_plot['Distance'])+3))
plt.ylabel("Distance (Å)", fontsize=12)
plt.tick_params(axis='x', labelsize=12)
plt.tick_params(axis='y', labelsize=12)
p.set_title('Distances in the simulation', fontsize=14)
sns.move_legend(p, "upper right")
plt.setp(p.get_legend().get_texts(), fontsize='10') # for legend text
plt.setp(p.get_legend().get_title(),fontsize='11')

## save the plot:
fig = p.get_figure()
fig.savefig(dic_input['out_file'], dpi=300) 


#plott = (
#    so.Plot(distances_to_plot, x="time_ns", y="Distance", color="Atoms")
#    .add(so.Line(linewidth=0.7))
#    .label(x="Time (ns)", y="Distance (Å)")
#    .theme(sns.axes_style("whitegrid"))
#    .save(dic_input['out_file'], dpi=300)
#)
