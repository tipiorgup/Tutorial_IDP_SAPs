# Simulation notebook 

`Simulation.ipynb` is a Google Collab notebook for running simple simulations without requiring an HPC. It utilizes the main codes in `src/` and `clean_pair_code.py`.

The interactive notebook depends on `global_vars.py` to create interactive widgets and `traj_process.py` to save final snapshots for further processing.
## Notes  
`Simulation.ipynb` is optimized for simulations with 1,000 residues.
For co-assembly simulations, input the ratio in the format `1-9`, `2-8`, etc., ensuring the total always sums to 10.

For creating the **SLATM represenation**, save locally the `stored_values.pkl` and `ensembles.pkl`.
