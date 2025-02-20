import numpy as np
try:
    import openmm 
    import openmm.app as app
    import openmm.unit as unit
except ImportError:
    import simtk.openmm as mm
    import simtk.openmm.app as app
    import simtk.unit as unit
    from simtk.openmm import XmlSerializer

def calculate_yukawa_parameters(temp, ionic_strength):
    """
    Calculate parameters for the Yukawa potential based on temperature and ionic strength.
    
    Args:
        temp (float): Temperature in Kelvin
        ionic_strength (float): Ionic strength of the solution
        
    Returns:
        tuple: (kT, lB, yukawa_kappa) - thermal energy, Bjerrum length, and screening parameter
    """
    kT = 8.3145 * temp * 1e-3  # Thermal energy in kJ/mol
    
    # Calculate dielectric constant of water
    fepsw = lambda T: 5321/T + 233.76 - 0.9297*T + 0.1417*1e-2*T*T - 0.8292*1e-6*T**3
    epsw = fepsw(temp)
    
    # Calculate Bjerrum length
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    
    # Calculate Debye screening parameter
    yukawa_kappa = np.sqrt(8*np.pi*lB*ionic_strength*6.022/10)
    
    return kT, lB, yukawa_kappa

def setup_system(box_size, z_size, peptides, chains, residue_data, temp, ph, yukawa_kappa):
    """
    Initialize the OpenMM simulation system with specified parameters.
    
    Args:
        box_size (float): Size of simulation box in x and y dimensions
        z_size (float): Size of simulation box in z dimension
        peptides (list): List of peptide sequences
        chains (list): Number of chains for each peptide
        residue_data (pd.DataFrame): Amino acid properties
        temp (float): Temperature in Kelvin
        ph (float): pH of the solution
        yukawa_kappa (float): Debye screening parameter
        
    Returns:
        tuple: (system, hb, ah, yu) - OpenMM system and force objects
    """
    system = openmm.System()
    
    # Set periodic box vectors
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = box_size * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = box_size * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = z_size * unit.nanometers
    system.setDefaultPeriodicBoxVectors(a, b, c)
    
    # Define Lennard-Jones parameters
    lj_eps = 0.2 * 4.184  # kJ/mol
    cutoff = 2  # nm
    
    # Setup force objects
    hb = openmm.HarmonicBondForce()  # Bonded interactions
    
    # Custom nonbonded force for Lennard-Jones interactions
    energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l))'
    ah = openmm.CustomNonbondedForce(energy_expression + '; s=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')
    
    # Add parameters to forces
    ah.addGlobalParameter('eps', lj_eps * unit.kilojoules_per_mole)
    ah.addGlobalParameter('rc', cutoff * unit.nanometer)
    ah.addPerParticleParameter('s')  # Particle size
    ah.addPerParticleParameter('l')  # Lambda parameter
    
    # Setup Yukawa force for electrostatic interactions
    yu = openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r-shift); q=q1*q2')
    yu.addGlobalParameter('kappa', yukawa_kappa/unit.nanometer)
    yu.addGlobalParameter('shift', np.exp(-yukawa_kappa*4.0)/4.0/unit.nanometer)
    yu.addPerParticleParameter('q')  # Charge
    
    # Configure force groups and methods
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
    
    return system, hb, ah, yu