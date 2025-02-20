from argparse import ArgumentParser
import pickle

def parse_arguments():
    """
    Set up command line argument parsing with detailed options for simulation control.
    
    Returns:
        ArgumentParser: Configured argument parser with simulation parameters
    """
    parser = ArgumentParser(description='Peptide simulation with customizable parameters')
    
    # Group for mutually exclusive peptide input methods
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--pickle-file', type=str,
                           help='Path to pickle file containing peptide pairs')
    input_group.add_argument('--peptide-pair', type=str,
                           help='Single peptide pair in format "XXX-YYY"')
    input_group.add_argument('--pairs-file', type=str,
                           help='Text file containing peptide pairs, one per line in format "XXX-YYY"')
    
    # Configuration for peptide selection
    parser.add_argument('--pair-index', type=int,
                       help='Index of pair to use from pickle file (if using pickle input)')
    parser.add_argument('--run-all', action='store_true',
                       help='Run simulation for all pairs in pickle/text file')
    
    # Simulation environment parameters
    parser.add_argument('--ratio', type=str, default='100-900',
                      help='Ratio of peptides in format "X-Y" (default: "100-900")')
    parser.add_argument('--ph', type=float, default=7.4,
                      help='pH value (default: 7.4)')
    parser.add_argument('--ionic-strength', type=float, default=0.0097,
                      help='Ionic strength (default: 0.0097)')
    
    # Box dimensions for simulation
    parser.add_argument('--box-size', type=float, default=15.0,
                      help='Box size in nanometers (default: 15.0)')
    parser.add_argument('--z-size', type=float, default=25.0,
                      help='Z-dimension size in nanometers (default: 25.0)')
    
    # Physical parameters
    parser.add_argument('--temperature', type=float, default=310,
                      help='Temperature in Kelvin (default: 310)')
    parser.add_argument('--steps', type=int, default=200000,
                      help='Number of simulation steps (default: 200000)')
    
    # Input/Output configuration
    parser.add_argument('--output-dir', type=str, default='.',
                      help='Output directory (default: current directory)')
    parser.add_argument('--residues-file', type=str, 
                      required=True,
                      help='Path to the residues CSV file containing amino acid properties')
    
    return parser.parse_args()

def load_peptide_pairs(args):
    """
    Load peptide pairs from specified input source (pickle, file, or command line).
    
    Args:
        args: Parsed command line arguments
        
    Returns:
        list: List of peptide pairs to simulate
        
    Raises:
        ValueError: If input parameters are invalid or missing
    """
    if args.pickle_file:
        print(f"Loading peptide pairs from pickle file: {args.pickle_file}")
        with open(args.pickle_file, 'rb') as f:
            flat_list = pickle.load(f)
            print(f"Loaded {len(flat_list)} pairs from pickle file")
            if not args.run_all:
                if args.pair_index is None:
                    raise ValueError("Must specify --pair-index when using pickle file without --run-all")
                return [flat_list[args.pair_index]]
            return flat_list
    elif args.pairs_file:
        with open(args.pairs_file, 'r') as f:
            pairs = [line.strip() for line in f if line.strip()]
            if not args.run_all:
                if args.pair_index is None:
                    raise ValueError("Must specify --pair-index when using pairs file without --run-all")
                return [pairs[args.pair_index]]
            return pairs
    else:  # Single pair from command line
        if args.peptide_pair is None:
            raise ValueError("No peptide pair specified")
        return [args.peptide_pair]

def validate_peptide_pair(pair):
    """
    Validate the format of a peptide pair string.
    
    Args:
        pair (str): Peptide pair in format "XXX-YYY"
        
    Returns:
        tuple: (peptide1, peptide2) if valid
        
    Raises:
        ValueError: If pair format is invalid
    """
    if not isinstance(pair, str) or '-' not in pair:
        raise ValueError(f"Invalid peptide pair format: {pair}. Expected format: 'XXX-YYY'")
    peptide1, peptide2 = pair.split('-')
    if not (peptide1.isalpha() and peptide2.isalpha()):
        raise ValueError(f"Invalid peptide sequence in pair: {pair}")
    return peptide1, peptide2