import os
import numpy as np
import pandas as pd
import itertools
try:
    import openmm 
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit

from .system_setup import calculate_yukawa_parameters, setup_system
from .utils import build_straight_CA_chain, write_pdb, insert_molecules

def run_single_simulation(rp, args):
    """Run simulation for a single peptide pair
    
    Args:
        rp (str): Raw peptide pair string in format "XXX-YYY"
        args: Command line arguments containing simulation parameters
    """
    try:
        peptides = rp.split('-')
        print(f"Processing pair: {rp}")
        
        # Setup directory structure
        sim_dir = os.path.join(args.output_dir, rp)
        sub_folder = os.path.join(rp, '310')
        os.makedirs(rp, exist_ok=True)
        os.makedirs(sub_folder, exist_ok=True)

        # Parse ratio for peptides
        ratio1, ratio2 = map(int, args.ratio.split('-'))
        
        # Load residues data
        residues = pd.read_csv(args.residues_file).set_index('three', drop=False)
        letters = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        
        # Generate initial structures for each peptide
        for sequence in peptides:
            n_residues = len(sequence)
            ca_pdb = os.path.join(rp, f'Peptide_{sequence}.pdb')
            ca_atoms = build_straight_CA_chain(sequence, r0=0.38)
            ca_atoms.loc[:, 'chainID'] = letters[peptides.index(sequence)]
            write_pdb(ca_atoms, ca_pdb)
        
        # Insert molecules with specified ratios
        insert_molecules(
            os.path.join(rp, f'Peptide_{peptides[0]}.pdb'),
            os.path.join(rp, 'tmp1.pdb'),
            n_mol=int(ratio1/len(peptides[0])),
            box=[15, 15, 15]
        )
        
        insert_molecules(
            os.path.join(rp, f'Peptide_{peptides[1]}.pdb'),
            os.path.join(rp, 'start.pdb'),
            n_mol=int(ratio2/len(peptides[1])),
            existing_pdb=os.path.join(rp, 'tmp1.pdb'),
            box=[15, 15, 15]
        )
        
        # Check if simulation already exists
        if os.path.isfile(os.path.join(sub_folder, f'{rp}_slab.dcd')):
            print(f"Simulation already exists for {rp}, skipping...")
            return
            
        # Setup residue information
        r = residues.copy()
        r = r.set_index('one')
        
        fasta = []
        fasta.append(list(peptides[0]))
        fasta.append(list(peptides[1]))
        
        # Modify terminal residues
        r.loc['X'] = r.loc[fasta[0][0]]
        r.loc['Z'] = r.loc[fasta[0][-1]]
        r.loc['X', 'MW'] += 2
        r.loc['Z', 'MW'] += 16
        
        fasta[0][0] = 'X'
        fasta[0][-1] = 'Z'
        fasta[1][0] = 'X'
        fasta[1][-1] = 'Z'
        
        # Setup simulation parameters
        lj_eps = 1.55 * 4.184
        types = list(np.unique(list(itertools.chain(*peptides))))
        types.extend(['X', 'Z'])
        
        MWs = [r.loc[a, 'MW'] for a in types]
        kT = 8.3145 * args.temperature * 1e-3
        
        # Setup pH dependent charges
        r.loc['H', 'q'] = 1. / (1 + 10**(args.ph - 6))
        r.loc['X', 'q'] = r.loc[fasta[0][0], 'q'] + 1.
        r.loc['Z', 'q'] = r.loc[fasta[0][-1], 'q'] - 1.
        
        # Calculate Yukawa parameters
        fepsw = lambda T: 5321/T + 233.76 - 0.9297*T + 0.1417*1e-2*T*T - 0.8292*1e-6*T**3
        epsw = fepsw(args.temperature)
        lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
        
        yukawa_eps = []
        yukawa_eps.append([r.loc[a].q*np.sqrt(lB*kT) for a in fasta[0]])
        yukawa_eps.append([r.loc[a].q*np.sqrt(lB*kT) for a in fasta[1]])
        
        yukawa_kappa = np.sqrt(8*np.pi*lB*args.ionic_strength*6.022/10)
        
        # Setup OpenMM system
        system = openmm.System()
        
        # Set box vectors
        L = 15.
        a = unit.Quantity(np.zeros([3]), unit.nanometers)
        a[0] = L * unit.nanometers
        b = unit.Quantity(np.zeros([3]), unit.nanometers)
        b[1] = L * unit.nanometers
        c = unit.Quantity(np.zeros([3]), unit.nanometers)
        c[2] = 25 * unit.nanometers
        system.setDefaultPeriodicBoxVectors(a, b, c)
        
        # Setup chains
        N = [len(peptides[0]), len(peptides[1])]
        chains = [int(ratio1/len(peptides[0])), int(ratio2/len(peptides[1]))]
        
        # Load PDB
        pdb = app.pdbfile.PDBFile(os.path.join(rp, 'start.pdb'))
        cutoff = 2
        
        # Add particles
        for _ in range(chains[0]):
            system.addParticle((r.loc[list(peptides[0])[0]].MW+2)*unit.amu)
            system.addParticle((r.loc[list(peptides[0])[-1]].MW+16)*unit.amu)
            for a in list(peptides[0])[1:-1]:
                system.addParticle(r.loc[a].MW*unit.amu)
                
        for _ in range(chains[1]):
            system.addParticle((r.loc[list(peptides[1])[0]].MW+2)*unit.amu)
            system.addParticle((r.loc[list(peptides[1])[-1]].MW+16)*unit.amu)
            for a in list(peptides[1])[1:-1]:
                system.addParticle(r.loc[a].MW*unit.amu)
        
        # Setup forces
        hb = openmm.HarmonicBondForce()
        
        # Setup nonbonded forces
        energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l))'
        ah = openmm.CustomNonbondedForce(energy_expression + '; s=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')
        
        ah.addGlobalParameter('eps', lj_eps*unit.kilojoules_per_mole)
        ah.addGlobalParameter('rc', cutoff*unit.nanometer)
        ah.addPerParticleParameter('s')
        ah.addPerParticleParameter('l')
        
        yu = openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r-shift); q=q1*q2')
        yu.addGlobalParameter('kappa', yukawa_kappa/unit.nanometer)
        yu.addGlobalParameter('shift', np.exp(-yukawa_kappa*4.0)/4.0/unit.nanometer)
        yu.addPerParticleParameter('q')
        
        # Add particles to forces
        for _ in range(chains[0]):
            for a, e in zip(list(fasta[0]), yukawa_eps[0]):
                yu.addParticle([e*unit.nanometer*unit.kilojoules_per_mole])
                ah.addParticle([r.loc[a].sigmas*unit.nanometer, r.loc[a].lambdas*unit.dimensionless])
                
        for _ in range(chains[1]):
            for a, e in zip(list(fasta[1]), yukawa_eps[1]):
                yu.addParticle([e*unit.nanometer*unit.kilojoules_per_mole])
                ah.addParticle([r.loc[a].sigmas*unit.nanometer, r.loc[a].lambdas*unit.dimensionless])
        
        # Add bonds and exclusions
        for j in range(chains[0]):
            begin = j*N[0]
            end = j*N[0]+N[0]
            for i in range(begin, end-1):
                hb.addBond(i, i+1, 0.38*unit.nanometer, 8033.0*unit.kilojoules_per_mole/(unit.nanometer**2))
                yu.addExclusion(i, i+1)
                ah.addExclusion(i, i+1)
        
        border = int(chains[0]*N[0])
        
        for j in range(chains[1]):
            begin = j*N[1]+border
            end = j*N[1]+N[1]+border
            for i in range(begin, end-1):
                hb.addBond(i, i+1, 0.38*unit.nanometer, 8033.0*unit.kilojoules_per_mole/(unit.nanometer**2))
                yu.addExclusion(i, i+1)
                ah.addExclusion(i, i+1)
        
        # Setup force groups and methods
        yu.setForceGroup(0)
        ah.setForceGroup(1)
        yu.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        ah.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
        hb.setUsesPeriodicBoundaryConditions(True)
        yu.setCutoffDistance(4*unit.nanometer)
        ah.setCutoffDistance(cutoff*unit.nanometer)
        
        # Add forces to system
        system.addForce(hb)
        system.addForce(yu)
        system.addForce(ah)
        
        # Setup integrator and simulation
        integrator = openmm.LangevinIntegrator(
            args.temperature*unit.kelvin,
            0.01/unit.picosecond,
            0.01*unit.picosecond
        )
        
        # Create and setup simulation
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        
        # Add reporters
        simulation.reporters.append(
            app.DCDReporter(
                os.path.join(sub_folder, f'{rp}_slab.dcd'),
                1000
            )
        )
        simulation.reporters.append(
            app.StateDataReporter(
                os.path.join(sub_folder, f'{rp}_slab.log'),
                1000,
                step=True,
                speed=True,
                elapsedTime=True,
                separator='\t'
            )
        )
        
        # Run simulation
        simulation.step(args.steps)
        
    except Exception as e:
        print(f"Error processing pair {rp}: {str(e)}")
        raise
        