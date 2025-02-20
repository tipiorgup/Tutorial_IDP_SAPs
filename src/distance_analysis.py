import mdtraj as md
import numpy as np
from mpl_toolkits.mplot3d import axes3d   
import matplotlib.pyplot as plt  
from matplotlib.figure import Figure 
import re
import itertools
import pickle
import os
import peptides
from numpy import savetxt
from collections import defaultdict
try:
    import qml
    from qml.kernels import gaussian_kernel
    from qml.kernels import laplacian_kernel
    from qml.representations import get_slatm_mbtypes
    print("QML version:",qml.__version__)
except ImportError:
    print("Failed to find QML")
from MDAnalysis.analysis import polymer



#three_letter to one_letter aminoacid

three_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

def shorten(x):
    if len(x) % 3 != 0: 
        raise ValueError('Input length should be a multiple of three')
    y = ''
    for i in range(len(x) // 3):
        y += three_one[x[3 * i : 3 * i + 3]]
    return y

def chain_array(number_chains,traj):
    """Big array of structures.
    int(end.xyz.shape[1]/50) is the number of residues
    new is an array with all peptides, residues, xyz info
    Function returns a list of the names of residues, list of xyz postion for all chains ans list of the last with name of residue
    """
    xyz=traj.xyz
    n_residues=traj.n_residues
    names=np.array([residue for residue in traj.topology.residues])
    new=np.reshape(xyz, [number_chains,int(xyz.shape[1]/number_chains),3])
    new_labeled= np.reshape(np.append(np.reshape(names,[1,xyz.shape[1],1]), xyz, axis=-1),[number_chains,int(xyz.shape[1]/number_chains),4])
    return new, new_labeled


def get_sequence(number_chains,traj):
    #Gives a list of residues three_letter
    res=[]
    res=[residue for residue in traj.topology.residues][:int(traj.xyz.shape[1]/number_chains)]
    return res

def get_charges(sequence_list):
    #Info about residues, gives charges list and total charge
    three_letter=[]
    charges=[]
    res_simple=[str(ele)[:3] for ele in sequence_list]
    for i in range(len(sequence_list)):
        three_letter.append(shorten(res_simple[i]))
        charges.append(peptides.Peptide(three_letter[i]).charge(pH=7))
    sequence=''.join(three_letter)
    Total_charge=peptides.Peptide(sequence).charge(pH=7)
    return charges, Total_charge,three_letter


def possible_interactions(number_chains, number_residues):
    #All posible peptide pair interactions
    arr = list(itertools.combinations(range(0,number_chains), 2))
    arr_between_res_2= list(itertools.combinations(range(0,number_residues), 2))+list(zip(range(number_residues), range(number_residues)))
    return arr, arr_between_res_2

def all_possible_distances(arr,arr_between_res_2,xyz_file):
    mega_dist=[]
    for s in range (len(arr)):
        temp_dist=[]
        for j in range (len(arr_between_res_2)):
            first_res=xyz_file[arr[s][0]][arr_between_res_2[j][0]]
            secon_res=xyz_file[arr[s][1]][arr_between_res_2[j][1]]
            temp_dist.append(np.linalg.norm(first_res-secon_res))
        mega_dist.append(temp_dist)
    #We created a mega list with all pair intercation between peptides and residues
    #First element is the pair of chains followed by pair of residues
    return mega_dist

def filter_cutoff(distance_array,arr,cutoff):
    "Gives the positions of the distances inside the cutoff and then locates the pair of clusters, gives a dictionary of interacting chains and its keys"""
    array=np.array(distance_array)
    indexArr = np.argwhere(array< cutoff).tolist()
    lpc=np.array(indexArr).T[0]
    clusters=[]
    for e in lpc:
        clusters.append(arr[e])
    clean_clusters=list({*map(tuple, map(sorted, clusters))})
    cluster=defaultdict(list)
    for v, k in clean_clusters:
        cluster[v].append(k)
    chain_clusters_keys=[*cluster]
    return cluster,chain_clusters_keys

def filter_clusters(n,cluster_dict):
    # Lets filter the clusters grater than n
    new_keys=[]
    chain_clusters=[*cluster_dict]
    for chain_clusters in cluster_dict:
        if len(cluster_dict[chain_clusters])>int(n-1):
            new_keys.append(chain_clusters)
    interesting = {x:cluster_dict[x] for x in new_keys}
    return interesting

def show_interacting_aa(cluster_dict,arr,mega_dist,res,arr_between_res_2,cutoff):
    #Positions in arr of the list of pairs
    keys=[*cluster_dict]
    list_positions={}
    for v in keys:
        tt=[]
        for s in cluster_dict[v]:
            tt.append(arr.index((v,s)))
            list_positions[v]=tt
    #Positions of the residue pairs is arr_2        
    to_check={}
    for v in keys:
        part=[]
        for s in list_positions[v]:
            part.append(np.argwhere(np.array(mega_dist[s])< cutoff).tolist()) #Find the pair of residues
        to_check[v]=part
    #See the name of the aminoacid
    aa={}
    for v in keys:
        zipp1=[]
        zipp2=[]
        zippy=[]
        for u in range(len(to_check[v])):
            pair1=[]
            pair2=[]
            for i in list(itertools.chain(*to_check[v][u])):
                pair1.append(res[arr_between_res_2[i][0]])
                pair2.append(res[arr_between_res_2[i][1]])
            zipp1.append(pair1)
            zipp2.append(pair2)
            zipped = list(zip(pair1,pair2))
            zippy.append(zipped)    
        aa[v]=zippy
    return aa