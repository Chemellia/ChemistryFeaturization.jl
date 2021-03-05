# Terminology

There are a lot of seemingly similar terms used for quantities in this package that refer to disparate things. Here, we try to best define these terms.

## Feature

A quality or quantity associated with an atom that we wish to encode, such as atomic mass, row in the periodic table, etc.

## Feature vector

The encoding (typically one-hot-style) of the values of a set of features associated with a particular atom (or, in the case of Weave featurization, pair of atoms, etc.).

For example, suppose we are encoding the atomic mass (across five possible bins) and periodic table block (s, p, d, or f) of Hydrogen. \
In this case, the associated feature vector would be `[1 0 0 0 0 1 0 0 0]` - where the first five indexes correspond to atomic mass and the last four correspond to the block in the periodic table.

## Featurization

Either the process of assigning feature vectors to chemical elements, or a description of a scheme for doing so, encoded as a Vector of *AtomFeat* objects (eventually, we plan to implement BondFeat, PairFeat, etc...).

## Feature matrix

A collection of feature vectors associated with the atoms in a structure. \
The feature matrix must be of dimensions `m x n` - where *m* is the number of features, and *n* is the number of nodes.
