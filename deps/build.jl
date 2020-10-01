using Conda
using Pkg

ENV["PYTHON"]=""
Pkg.build("PyCall")
using PyCall

Conda.add("ase", channel="conda-forge")
Conda.add("rdkit", channel="conda-forge")
Conda.add("pymatgen", channel="conda-forge")