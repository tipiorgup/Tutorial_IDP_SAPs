import numpy as np
import pandas as pd
from MDAnalysis.lib.nsgrid import FastNS
from MDAnalysis.lib.distances import distance_array
from scipy.spatial.transform import Rotation as R
import numpy as np
import pandas as pd
try:
    import pdbfixer
except Exception as e:
    print('Cannot import pdbfixer. This may affect 3SPN2 model.')
try:
    import openmm.unit as unit
except ImportError:
    import simtk.unit as unit
import sys
import os
_amino_acid_1_letter_to_3_letters_dict = dict(A='ALA', R='ARG', N='ASN', D='ASP', C='CYS', 
                                              Q='GLN', E='GLU', G='GLY', H='HIS', I='ILE', 
                                              L='LEU', K='LYS', M='MET', F='PHE', P='PRO', 
                                              S='SER', T='THR', W='TRP', Y='TYR', V='VAL')

def build_straight_CA_chain(sequence, r0=0.38):
    """
    Build a straight portein CA atom chain with given sequence. 
    
    Parameters
    ----------
    sequence : str
        Protein chain sequence (1 letter for each amino acid). 
    
    r0: float or int
        Distance in unit nm between two neighboring CA atoms. 
    
    Returns
    -------
    df_atoms : pd.DataFrame
        A pandas dataframe includes atom information. 
    
    """
    n_atoms = len(sequence)
    data = []
    for i in range(n_atoms):
        resname = _amino_acid_1_letter_to_3_letters_dict[sequence[i]]
        atom_i_dict = {'recname': 'ATOM', 'name': 'CA', 'altLoc': '', 'resname': resname, 'chainID': 'A', 
                       'iCode': '', 'occupancy': 1.0, 'tempFactor': 1.0, 'element': 'C', 'charge': ''}
        data.append(atom_i_dict)
    df_atoms = pd.DataFrame(data)
    df_atoms['serial'] = list(range(1, n_atoms + 1))
    df_atoms['resSeq'] = list(range(1, n_atoms + 1))
    df_atoms.loc[:, 'x'] = 0
    df_atoms.loc[:, 'y'] = 0
    z = r0*np.arange(n_atoms)
    z -= np.mean(z)
    df_atoms['z'] = z*10 # convert nm to angstroms
    df_atoms['z'] = df_atoms['z'].round(3)
    return df_atoms

def write_pdb(pdb_atoms, pdb_file, write_TER=False):
    """
    Write pandas dataframe to pdb file. 
    
    Parameters
    ----------
    pdb_atoms : pd.DataFrame
        A pandas dataframe includes atom information. 
    
    pdb_file : str
        Output path for the pdb file. 
    
    write_TER : bool
        Whether to write TER between two chains. 

    """
    chainID = None
    with open(pdb_file, 'w') as pdb:
        for i, atom in pdb_atoms.iterrows():
            if chainID is not None:
                if write_TER and (atom['chainID'] != chainID):
                    pdb.write('TER\n')
            chainID = atom['chainID']
            pdb_line = f'{atom.recname:<6}{int(atom.serial):>5} {atom["name"]:^4}{atom.altLoc:1}'+\
                       f'{atom.resname:<3} {atom.chainID:1}{int(atom.resSeq):>4}{atom.iCode:1}   '+\
                       f'{atom.x:>8.3f}{atom.y:>8.3f}{atom.z:>8.3f}' +\
                       f'{atom.occupancy:>6.2f}{atom.tempFactor:>6.2f}'+' ' * 10 +\
                       f'{atom.element:>2}{atom.charge:>2}'
            assert len(pdb_line) == 80, f'An item in the atom table is longer than expected ({len(pdb_line)})\n{pdb_line}'
            pdb.write(pdb_line + '\n')
        pdb.write('END\n')


"""
A python script for inserting molecules, similar to gmx insert-molecules. 
"""

def insert_molecules(new_pdb, output_pdb, n_mol, radius=0.5, existing_pdb=None, max_n_attempts=10000, 
                     box=[100, 100, 100], method='FastNS', reset_serial=True):
    """
    Insert multiple copies of given molecule into an existing pdb or a new empty box. 
    Note all the length parameter units are nm, though coordinates in pdb are in unit angstroms. 
    Currently this function only supports using orthogonal box. 
    
    Parameters
    ----------
    new_pdb : str
        Path of the new pdb file to be inserted. 
    
    output_pdb : str
        Path of the output pdb file. 
    
    n_mol : int
        Number of copies of the new molecule to be inserted. 
    
    radius : float or int
        Radius for each atom in unit nm. Insertion ensures the distance between any two atoms are larger than 2*radius under periodic boundary condition. 
    
    existing_pdb : None or str
        Path of an existing pdb file. If None, the molecules will be inserted into a new empty box, otherwise, inserted into a copy of existing_pdb. 
    
    max_n_attempts : int
        Maximal number of attempts to insert. Ensure this value no smaller than n_mol. 
    
    box : list or numpy array, shape = (3,)
        New box size in unit nm. 
    
    method : str
        Which method to use for detecting if the inserted molecule has contact with existing atoms. 
        This parameter can be 'FastNS' or 'distance_array'. 
        FastNS should be faster than distance_array, especially for large systems. 
    
    reset_serial : bool
        Whether to reset serial to 0, 1, ..., N - 1. 
        If True, the serial in the final pdb is reset as 0, 1, ..., N - 1. 
        If False, the serial remains unchanged. 
        If True and atom number > 1,000,000, the serial remains unchanged since the largest atom serial number allowed in pdb is 999,999. 
    
    """
    assert method in ['FastNS', 'distance_array'] # check method
    assert max_n_attempts >= n_mol
    print(f'Check contact with {method} method. ')
    if existing_pdb is None:
        atoms = pd.DataFrame()
    else:
        atoms = helper_functions.parse_pdb(existing_pdb)
    new_atoms = helper_functions.parse_pdb(new_pdb)
    new_coord = new_atoms[['x', 'y', 'z']].to_numpy()
    new_coord -= np.mean(new_coord, axis=0)
    count_n_mol = 0
    count_n_attempts = 0
    cutoff = float(2*10*radius) # convert nm to angstrom
    box_a, box_b, box_c = 10*box[0], 10*box[1], 10*box[2] # convert nm to angstrom
    dim = np.array([box_a, box_b, box_c, 90.0, 90.0, 90.0])
    if method == 'FastNS':
        dim = dim.astype(np.float32)
        if len(atoms.index) > 0:
            coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
            grid_search = FastNS(cutoff, coord, dim, pbc=True)
    while (count_n_mol < n_mol) and (count_n_attempts < max_n_attempts):
        # get a random rotation
        rotate = R.random()
        new_coord_i = rotate.apply(new_coord)
        # get a random translation
        translate = np.random.uniform(0, 1, 3)*np.array([box_a, box_b, box_c])
        new_coord_i += translate
        if len(atoms.index) == 0:
            new_atoms_i = new_atoms.copy()
            new_atoms_i[['x', 'y', 'z']] = new_coord_i
            atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
            count_n_mol += 1
            if method == 'FastNS':
                coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
                grid_search = FastNS(cutoff, coord, dim, pbc=True)
        else:
            flag = False
            if method == 'distance_array':
                coord = atoms[['x', 'y', 'z']].to_numpy()
                d = distance_array(coord, new_coord_i, dim)
                if np.amin(d) >= cutoff:
                    flag = True
            elif method == 'FastNS':
                results = grid_search.search(new_coord_i.astype(np.float32))
                if len(results.get_pair_distances()) == 0:
                    flag = True
            if flag:
                new_atoms_i = new_atoms.copy()
                new_atoms_i[['x', 'y', 'z']] = new_coord_i
                atoms = pd.concat([atoms, new_atoms_i], ignore_index=True)
                count_n_mol += 1
                if method == 'FastNS':
                    coord = atoms[['x', 'y', 'z']].to_numpy().astype(np.float32)
                    grid_search = FastNS(cutoff, coord, dim, pbc=True)
        count_n_attempts += 1
    
    # determine if the number of molecules were successfully added:
    if count_n_mol == n_mol:
        print(f'Successfully inserted {n_mol} molecules.')
    else:
        print(f'Could not successfully insert {n_mol} molecules in {count_n_attempts} attempts.')
        print(f'Only added {count_n_mol} molecules. Try increasing the box size or number of attempts to add more molecules.')
    if reset_serial:
        n_atoms = len(atoms.index)
        if n_atoms > 1000000:
            print(f'Too many atoms. Cannot reset serial as 0, 1, ..., N - 1. Serial remains unchanged.')
        else:
            atoms['serial'] = list(range(len(atoms.index)))
    # write the final pdb
    helper_functions.write_pdb(atoms, output_pdb)

