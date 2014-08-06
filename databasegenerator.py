import numpy as np
import pymatgen as pm
import pandas as pd
import os
from operator import mul

def csvbuilder(target_dir, output, column_functions):
    """
    Build a csv from structure files.

    Keyword arguments:
    target_dir -- string of path to look for structure files
    output -- filename for output csv
    column_functions -- list of function that build the csv.  Each function should take a pymatgen Structure as a parmeter and return a single numeric or categorical (string or character) value
    """
    structure_files = os.listdir(target_dir)

    df = pd.DataFrame()

    for structure_filename in structure_files:
        a = pm.read_structure(os.path.join(target_dir, structure_filename))
        next_row = len(df.index)
        if next_row > 0: 
            row = [f(a) for f in column_functions]
            df.loc[next_row] = row
        else:
            colnames = [f.func_name for f in column_functions]
            row = [f(a) for f in column_functions]
            df = pd.DataFrame([row], columns=colnames)
    df.to_csv(output)
    return df

#list of column functions you might want. Of course you can use any function that accepts a Structure object and returns a single value
def formula(a):
    c = a.composition
    return c.alphabetical_formula

def numberOfSpecies(a):
    return len(a.composition)

def density(a):
    return a.density

def fracTransitionMetal(a):
    c = a.composition
    return sum(map(mul, [x.is_transition_metal for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracRareEarth(a):
    c = a.composition
    return sum(map(mul, [x.is_rare_earth_metal for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracNobleGas(a):
    c = a.composition
    return sum(map(mul, [x.is_noble_gas for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracMetalloid(a):
    c = a.composition
    return sum(map(mul, [x.is_metalloid for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracLanthanoid(a):
    c = a.composition
    return sum(map(mul, [x.is_lanthanoid for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracHalogen(a):
    c = a.composition
    return sum(map(mul, [x.is_halogen for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracChalcogen(a):
    c = a.composition
    return sum(map(mul, [x.is_chalcogen for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracAlkaline(a):
    c = a.composition
    return sum(map(mul, [x.is_alkaline for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracAlkali(a):
    c = a.composition
    return sum(map(mul, [x.is_alkali for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def fracActinoid(a):
    c = a.composition
    return sum(map(mul, [x.is_actinoid for x in c.elements], map(c.get_atomic_fraction, c.elements)))

def electronsPerAtom(a):
    c = a.composition
    return sum(map(mul, map(c.get_atomic_fraction, c.elements),[x.common_oxidation_states[0] for x in c.elements]))

def electronegativityRange(a):
    return np.ptp([x.X for x in a.species])

def electronegativityStd(a):
    return np.std([x.X for x in a.species])

def radiiRange(a):
    return np.ptp([x.atomic_radius for x in a.species])

def radiiStd(a):
    return np.std([x.atomic_radius for x in a.species])

def rowRange(a):
    return np.ptp([x.row for x in a.species])

def rowStd(a):
    return np.std([x.row for x in a.species])

def ordered(a):
    return a.is_ordered

##the following are things we might want to predict/use in training
def crystalSystem(a):
    sf = pm.symmetry.finder.SymmetryFinder(a)
    return sf.get_crystal_system

def spaceGroup(a):
    sf = pm.symmetry.finder.SymmetryFinder(a)
    return sf.get_spacegroup_symbol()

def volumePerSite(a):
    return a.volume/a.num_sites
