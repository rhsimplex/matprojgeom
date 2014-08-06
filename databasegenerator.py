import numpy as np
import pymatgen as pm
import pandas as pd
import os

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

def formula(a):
    return a.formula

def numberOfSpecies(a):
    return len(a.composition)
