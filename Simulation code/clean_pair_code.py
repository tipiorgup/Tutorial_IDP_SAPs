#!/usr/bin/env python3

from src.parsers import parse_arguments, load_peptide_pairs
from src.simulation import run_single_simulation

def main():
    """
    Main entry point for peptide simulation program.
    Handles argument parsing and simulation execution.
    """
    # Parse command line arguments
    args = parse_arguments()
    
    # Load peptide pairs based on input method
    peptide_pairs = load_peptide_pairs(args)
    
    # Process each peptide pair
    for pair in peptide_pairs:
        try:
            run_single_simulation(pair, args)
        except Exception as e:
            print(f"Error processing pair {pair}: {str(e)}")
            if not args.run_all:
                raise
            continue

if __name__ == '__main__':
    main()