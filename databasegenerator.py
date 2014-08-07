import numpy as np
import pymatgen as pm
import csv
import os
from operator import mul

def csvbuilder(target_dir, output, column_functions, nary=10):
    """
    Build a csv from structure files.

    Keyword arguments:
    target_dir -- string of path to look for structure files
    output -- filename for output csv
    column_functions -- list of function that build the csv.  Each function should take a pymatgen Structure as a parmeter and return a single numeric or categorical (string or character) value
    nary -- limit to n-ary structures (e.g. nary = 2 only binary and elemental, nary = 3 only ternary, binary, and elemental, etc.)
    """
    structure_files = os.listdir(target_dir)
    with open(output, 'wb') as csvfile:
        next_row = 0
        structurewriter = csv.writer(csvfile)
        for structure_filename in structure_files:
            a = pm.read_structure(os.path.join(target_dir, structure_filename))
            if next_row > 0 and numberOfSpecies(a) < nary: 
                row = [f(a) for f in column_functions]
                structurewriter.writerow(row)
                next_row += 1
            elif numberOfSpecies(a) < nary:
                colnames = [f.func_name for f in column_functions]
                structurewriter.writerow(colnames)
                row = [f(a) for f in column_functions]
                structurewriter.writerow(row)
                next_row += 2

#list of column functions you might want. Of course you can use any function that accepts a Structure object and returns a single value
def formula(a):
    c = a.composition
    return c.alphabetical_formula

def numberOfSpecies(a):
    return len(a.composition)

def density(a):
    return round(a.density,4)

def fracTransitionMetal(a):
    c = a.composition
    return round(sum(map(mul, [x.is_transition_metal for x in c.elements], map(c.get_atomic_fraction, c.elements))),4)

def fracRareEarth(a):
    c = a.composition
    return round(sum(map(mul, [x.is_rare_earth_metal for x in c.elements], map(c.get_atomic_fraction, c.elements))),4)

def fracNobleGas(a):
    c = a.composition
    return round(sum(map(mul, [x.is_noble_gas for x in c.elements], map(c.get_atomic_fraction, c.elements))),4)

def fracMetalloid(a):
    c = a.composition
    return round(sum(map(mul, [x.is_metalloid for x in c.elements], map(c.get_atomic_fraction, c.elements))))

def fracLanthanoid(a):
    c = a.composition
    return round(sum(map(mul, [x.is_lanthanoid for x in c.elements], map(c.get_atomic_fraction, c.elements))))

def fracHalogen(a):
    c = a.composition
    return round(sum(map(mul, [x.is_halogen for x in c.elements], map(c.get_atomic_fraction, c.elements))))

def fracChalcogen(a):
    c = a.composition
    return round(sum(map(mul, [x.is_chalcogen for x in c.elements], map(c.get_atomic_fraction, c.elements))))

def fracAlkaline(a):
    c = a.composition
    return round(sum(map(mul, [x.is_alkaline for x in c.elements], map(c.get_atomic_fraction, c.elements))))

def fracAlkali(a):
    c = a.composition
    return round(sum(map(mul, [x.is_alkali for x in c.elements], map(c.get_atomic_fraction, c.elements))))

def fracActinoid(a):
    c = a.composition
    return round(sum(map(mul, [x.is_actinoid for x in c.elements], map(c.get_atomic_fraction, c.elements))),4)

def electronsPerAtom(a):
    c = a.composition
    try:
        return round(sum(map(mul, map(c.get_atomic_fraction, c.elements),[x.common_oxidation_states[0] for x in c.elements])),4)
    except IndexError:
        print 'Unable to compute e/a for ' + c.alphabetical_formula + '.'
        return -1000

def electronegativityRange(a):
    return round(np.ptp([x.X for x in a.species]),4)

def electronegativityStd(a):
    return round(np.std([x.X for x in a.species]),4)

def radiiRange(a):
    try:
        return round(np.ptp([x.atomic_radius for x in a.species]),4)
    except TypeError:
        print 'Unable to radii for ' + a.formula + '.'
        return -1000
    
def radiiStd(a):
    try:
        return round(np.std([x.atomic_radius for x in a.species]),4)
    except TypeError:
        print 'Unable to radii for ' + a.formula + '.'
        return -1000

def rowRange(a):
    return round(np.ptp([x.row for x in a.species]),4)

def rowStd(a):
    return round(np.std([x.row for x in a.species]),4)

def ordered(a):
    return a.is_ordered

##the following are things we might want to predict/use in training
def crystalSystem(a):
    if np.std([a.lattice.a, a.lattice.b, a.lattice.c]) > 20:
        print 'Long c axis detected for ' + a.formula + '. Skipping crystal system calculation.'
        return 'NA'
    sf = pm.symmetry.finder.SymmetryFinder(a)
    return sf.get_crystal_system()

def spaceGroup(a):
    if np.std([a.lattice.a, a.lattice.b, a.lattice.c]) > 20:
        print 'Long c axis detected for ' + a.formula + '. Skipping space group calculation.'
        return 'NA'
    sf = pm.symmetry.finder.SymmetryFinder(a)
    return sf.get_spacegroup_symbol()

def volumePerSite(a):
    return round(a.volume/a.num_sites,4)
