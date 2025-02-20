# Import main functionality to expose at package level
from .parsers import parse_arguments, load_peptide_pairs, validate_peptide_pair
from .system_setup import calculate_yukawa_parameters, setup_system
from .simulation import run_single_simulation
from .utils import build_straight_CA_chain, write_pdb, insert_molecules

# Define what should be available when using "from src import *"
__all__ = [
    'parse_arguments',
    'load_peptide_pairs',
    'validate_peptide_pair',
    'calculate_yukawa_parameters', 
    'setup_system',
    'run_single_simulation',
    'build_straight_CA_chain',
    'write_pdb',
    'insert_molecules'
]

# Package metadata
__description__ = 'Peptide simulation code using CALVADOS2 model for molecular dynamics'