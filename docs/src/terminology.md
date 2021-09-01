# Terminology and Philosophy

There are a lot of seemingly similar terms used for quantities in this package that refer to disparate things (or, are used slightly differently by other people in other places). Here, we try to best define these terms as we intend them. Further down, once the terms are defined, we elaborate on why the package is designed the way it is.

```@contents
Pages = ["terminology.md"]
Depth = 3
```

## General Terms

### Feature

A quality or a quantity that is associated with an atom, which is of interest to us. A feature could be anything, like atomic mass, the row in the periodic table to which the element belongs, etc.

### Encoding

The process of translating the value of a feature from its human-readable form (such as a `Float` or a `String`) to the form fed into the machine learning model.

The complexity of the encoding logic, whether it's something as simple as an equality operation, or something much more complex, is a decision left to the user. Often, this is building a one-hot vector.

### Decoding

The inverse process to encoding.

Generally, it is recommended that there is symmetry between the encoding and decoding logic.

!!! note
    In many cases, the process isn't fully invertible. For instance, for a continuous-valued feature encoded to a one-hot vector you can't get back a precise value but rather only a range corresponding to the associated onehot bin.
    In such cases, the decoding mechanism should try and return the most meaningful and human-interpretable form.

## Design Philosophy

Points to elaborate on here...

* maintaining transparency/decodability, as well as user's choice to include (or not) particular features, and choose how they are encoded, as opposed to a "black-box" scheme with little to no customizability
* separation of concerns/modularity...as many things should be "plug-and-play" with each other as possible (e.g. swapping in different codecs to FD's, different FD's to featurizations, different featurizations to atoms objects)
* FD's and atoms objects are generic and fairly reusable across models, featurizations are less so (closer to one-to-one relationship between a featurization type and a model type, but still have flexibility to easily include/exclude different features)
