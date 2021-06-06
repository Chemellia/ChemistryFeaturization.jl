# Ideas for Modulification

Quick reference/scribble pad for all things modulify-related (could also possibly be used in revamping the docs)

## Data Structures

Major data structures that are used in `ChemistryFeaturization`

```text
|---- AbstractAtoms
    |==== AtomGraph
    |==== WeaveMol
|---- AbstractFeatureDescriptor
    |---- AbstractAtomFeatureDescriptor
        |==== ElementFeatureDescriptor
        |==== SpeciesFeatureDescriptor
    |---- AbstractPairFeatureDescriptor
        |==== PairFeatureDescriptor
        |---- BondFeatureDescriptor
            |==== BondType
            |==== InRing
            |==== IsConjugated
    |---- AbstractEnvironmentFeatureDescriptor
        |==== OrbitalFieldMatrix
|---- AbstractFeaturization
    |==== GraphNodeFeaturization
    |==== WeaveFeaturization
|---- Bonds
    |==== Single
    |==== Double
    |==== Triple
```

### Note

* `----` - Abstract Type
* `====` - Concrete Type

## Dependencies of Concrete Types

``` text
# one type depending on multiple different types
* AtomGraph
    * AbstractAtoms
    * AbstractFeaturization

* WeaveMol
    * AbstractAtoms
    * WeaveFeaturization

* GraphNodeFeaturization
    * AbstractFeaturization
    * AbstractAtomFeatureDescriptor

* WeaveFeaturization
    * AbstractFeaturization
    * AbstractAtomFeatureDescriptor
    * AbstractPairFeatureDescriptor


# the following FDs depend only on one other type (their immediate parent abstract FD type)
* AbstractAtomFeatureDescriptor
    * ElementFeatureDescriptor
    * SpeciesFeatureDescriptor


* AbstractPairFeatureDescriptor
    * PairFeatureDescriptor
    * BondFeatureDescriptor

* AbstractEnvironmentFeatureDescriptor
    * OrbitalFieldMatrix

# these types are tightly bound together (by design)... so not really consequential in identifying dependencies.
* Bonds
    * BondType

* BondFeatureDescriptor
    * BondType
    * InRing
    * IsConjugated

```

## Functions

Dependencies some functions have. Organized on the basis of files they currently reside in.

Note: Functions mentioned here do not include constructors (those have been included in the above section). Dependencies don't include the concrete types defined in the file.

```text
atoms/atomgraph.jl
  - All functions use only AtomGraph

features/elementfeature.jl
  - function (f::ElementFeatureDescriptor)(a::AbstractAtoms)

features/speciesfeature.jl
  - function (f::SpeciesFeatureDescriptor)(a::AbstractAtoms)

featurizations/graphnodefeaturization.jl
  - function decode(ag::AtomGraph) (this can probably be redefined in some other file - atoms/atomgraph.jl maybe?)

```

## Actual Ideas/Thoughts/Things to think about

* Create a new module for Atoms, FDs, and Fzns. These will make available the different atoms, fds and fzns respectively.

* Handling Abstractions -
  * Create a new module for all abstract types. All other modules (different types of atoms, fds, and fzns) will be `using` the required types from this module. This module will be kept internal, and needs to be made available only to the other modules being created (as in, this won't be `include("")`ed in `ChemistyFeaturization.jl`).
  * Each module will also export (more like channel it through?) the corresponding Abstract types they use too.
  * Note
    * This level of abstraction could also be applied to more levels. For instance, an abstractfds module could be defined (which defines a new abstract type for an FD). Only the implementation file with the actual concrete types that inherit from this FD will `include("")` this module. The `FeatureDescriptor` module will export the abstract types defined.
    * Probably worth considering how this may have any drawbacks. Feels like this may potentially cause redefinition of variables (like we had at one point), which isn't a very nice thing to look at while precompiling/testing.

* Move all `function (f::FeatureDescriptor)(a::AbstractAtoms)` similar functions into the FDs module

* Encode and decode functions need to be defined for the respective `AbstractAtoms`/`AbstractFeaturization` types in the same files itself. Any definitions involving combinations can be generalized to be in either the main module or in one of the holy trio of submodules (atoms, FDs, and Fzns) depending upon context.\
For instance,
  * `function (f::FeatureDescriptor)(a::AbstractAtoms)` would be better suited in the main module
  * `function (f::FeatureDescriptor)(a::AtomGraph)` would be better suited in the Atoms module
  * `function (f::ElementFeatureDescriptor)(a::AbstractAtoms)` would be better suited in the FDs module

Note - this however raises the question, where would something like  `function (f::ElementFeatureDescriptor)(a::AtomGraph)` be better suited.\
Should we just make this a module of its own? (feels like I'm going a little overboard with the modules haha. starting to sound a lot like "you get a module, you get a module, everybody gets a module!!!")
