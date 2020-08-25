import numpy as np
import pandas as pd
from pymatgen.core import periodic_table as pt
from matminer.featurizers.composition import ElementProperty
from matminer.featurizers.conversions import StrToComposition

ep = ElementProperty.from_preset("matminer")

# list of elements
data = pt._pt_data.keys()
data = pd.DataFrame(data=data, columns=['symbol'])
data = StrToComposition().featurize_dataframe(data, 'symbol')
data = ep.featurize_dataframe(data, "composition", ignore_errors=True)

# only get ones that don't overlap with pymatgen features
features = ['mendeleev_no',
 'electrical_resistivity',
 'velocity_of_sound',
 'thermal_conductivity',
 'bulk_modulus',
 'coefficient_of_linear_thermal_expansion'] 