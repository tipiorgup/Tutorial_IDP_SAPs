import glob
import os
import mdtraj as md
import numpy as np
import pickle

def process_trajectory(path, sequence, ens_file):
    """
    Process molecular dynamics trajectory and create ensemble data.
    
    Args:
        path (str): Path to the directory containing trajectory files
        sequence (str): Sequence identifier
        ens_file (str): Path to save/load ensemble file
        
    Returns:
        dict: Ensemble data dictionary or None if processing fails
    """
    # Check if ensemble file already exists and is valid
    if os.path.isfile(ens_file):
        if os.path.getsize(ens_file) > 0:
            print('The ensemble has been done')
            return None
        else:
            # Remove empty file
            os.remove(ens_file)
            print(f"Removed empty ensemble file: {ens_file}")

    # Find DCD files
    txt_files = glob.glob(os.path.join(path, '310', '*dcd'), recursive=True)
    if not txt_files:
        print(f"No DCD files found in {path}")
        return None

    topology = os.path.join(path, 'start.pdb')
    if not os.path.isfile(topology):
        print(f"Topology file not found: {topology}")
        return None

    try:
        # Load trajectory
        trajectory = md.load(txt_files[0], top=topology)
        traj_length = len(trajectory)
        
        if traj_length < 50:
            print(f"Warning: Trajectory only has {traj_length} frames, need at least 50")
            return None

        # Process trajectory
        traj_r = trajectory[-50:]
        traj = {sequence: [traj_r[j] for j in range(50)]}

        # Create ensembles
        st = np.arange(5, 50, 10)
        e = np.arange(10, 60, 10)
        ensembles = {sequence: {}}

        for idx in range(len(st)):
            ensembles[sequence][idx] = traj[sequence][st[idx]:e[idx]]

        # Save ensembles
        try:
            with open(ens_file, 'wb') as handle:
                pickle.dump(ensembles, handle, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Successfully saved ensemble data to {ens_file}")
        except Exception as e:
            print(f"Error saving ensemble file: {e}")
            return None

        return ensembles

    except Exception as e:
        print(f"Error processing trajectory: {e}")
        return None