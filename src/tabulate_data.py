import pandas as pd
import numpy as np
from pymatgen.core import periodic_table as pt

# dictionary where keys are element symbols and values are dictionaries with data
data = pt._pt_data

# which properties to tabulate from big dictionary?
pmg_props = ['Atomic no', 'Name', 'Atomic mass', 'Atomic radius', 'Boiling point', 'Melting point', 'X']
temp_props = ["Boiling point", "Melting point"]

# set up
df_data = {p:[] for p in pmg_props}
df_data['Symbol'] = []

# pull that stuff
for k in data.keys():
    df_data['Symbol'].append(k)
    data_here = data[k]
    for prop in pmg_props:
        if prop in data_here.keys():
            prop_val = data_here[prop]
            if type(prop_val)==str and "no data" in prop_val:
                prop_val = np.nan
            elif prop in temp_props:
                spl = prop_val.split(" ")
                if len(spl)>2:
                    prop_val = float(spl[-2]) # "about" or "maybe"
                else:
                    prop_val = float(spl[0])
            elif prop=="Atomic radius":
                prop_val = float(prop_val)
            df_data[prop].append(prop_val)
        else:
            df_data[prop].append(np.nan)

df = pd.DataFrame.from_dict(df_data)

# special code to handle orbitals
orbitals = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f', '5s', '5p', '5d', '6s', '6p', '6d', '7s']
orb_data = {orb:[] for orb in orbitals}
for row in df.iterrows():
    sym = row[1]['Symbol']
    if 'Atomic orbitals' in data[sym].keys():
        data_here = data[sym]['Atomic orbitals']
        if data_here == 'no data':
            no_data = True
        else:
            no_data = False
            occ_orbs = list(data_here.keys())
            for orb in orbitals:
                if orb in occ_orbs:
                    orb_data[orb].append(data_here[orb])
                else:
                    orb_data[orb].append(np.nan)
    else:
        no_data = True
    if no_data:
        for orb in orbitals:
            orb_data[orb].append(np.nan)
for key in orb_data.keys():
    df[key] = orb_data[key]

df.sort_values(by="Atomic no", inplace=True)
df.reset_index(inplace=True, drop=True)
df = df.iloc[0:103]

# stuff that we need to instantiate Element to get: row, group, block valence

# my own valence function here until my merged PR goes live...
def valence(element):
    """
    # From full electron config obtain valence subshell
    # angular moment (L) and number of valence e- (v_e)

    """
    # the number of valence of noble gas is 0
    if element.group == 18:
        return (np.nan, 0)

    L_symbols = 'SPDFGHIKLMNOQRTUVWXYZ'
    valence = []
    full_electron_config = element.full_electronic_structure
    last_orbital = full_electron_config[-1]
    for n, l_symbol, ne in full_electron_config:
        l = L_symbols.lower().index(l_symbol)
        if ne < (2 * l + 1) * 2:
            valence.append((l, ne))
        # check for full last shell (e.g. column 2)
        elif (n, l_symbol, ne) == last_orbital and ne == (2 * l + 1) * 2 and len(valence) == 0:
            valence.append((l, ne))
    if len(valence) > 1:
        raise ValueError("Ambiguous valence")

    return valence[0]

rows = []
groups = []
valences = []
blocks = []
for row in df.iterrows():
    sym = row[1]['Symbol']
    el = pt.Element(sym)
    rows.append(el.row)
    groups.append(el.group)
    blocks.append(el.block)
    try:
        v = valence(el)[1]
    except ValueError: #ambiguous valence
        v = np.nan
    valences.append(v)

df['Row'] = rows
df['Group'] = groups
df['Valence'] = valences
df['Block'] = blocks

df.to_csv("../data/pymatgen_atom_data.csv")

# TODO: oxidation states too
# TODO: other data sources via matminer
