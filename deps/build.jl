using Conda

Conda.add("ase", channel = "conda-forge")
Conda.add("rdkit", channel = "conda-forge")
Conda.add("pymatgen", channel = "conda-forge")
Conda.rm("mkl")
Conda.add("mkl")
