import numpy as np
import pandas as pd
from pymatgen.core import periodic_table as pt
from matminer.featurizers.composition import ElementProperty
from matminer.featurizers.conversions import StrToComposition

# only get ones that don't overlap with pymatgen features
deml_features = ['molar_vol', 'heat_fusion', 'heat_cap', 'first_ioniz', 'electric_pol', 'GGAU_Etot']
ep_deml = ElementProperty(data_source="deml", features=deml_features, stats=['mean'])
magpie_features = ['NsValence', 'NpValence', 'NdValence', 'NfValence', 'NValence', 'NsUnfilled', 'NpUnfilled', 'NdUnfilled', 'NfUnfilled', 'NUnfilled', 'GSvolume_pa', 'GSbandgap', 'GSmagmom']
ep_magpie = ElementProperty(data_source="magpie", features=magpie_features, stats=['mean'])

# list of elements
data = pt._pt_data.keys()

# reformat and feed into matminer for deml features
data = pd.DataFrame(data=data, columns=['symbol'])
data = StrToComposition().featurize_dataframe(data, 'symbol')
data = ep_deml.featurize_dataframe(data, "composition", ignore_errors=True)
# clean up other column names
data.rename(columns={f:f[14:] for f in data.columns[2:]}, inplace=True)

# get magpie features
data = ep_magpie.featurize_dataframe(data, "composition", ignore_errors=True)

# other column name cleanup
data.rename(columns={f:f[16:] for f in data.columns[8:]}, inplace=True)
data.rename(columns={'symbol':"Symbol"}, inplace=True)
data.drop('composition', inplace=True, axis=1)

# save
data.to_csv("matminer_atom_data.csv")