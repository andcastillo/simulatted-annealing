# openchemlib-extended

[![NPM version][npm-image]][npm-url]
[![build status][travis-image]][travis-url]
[![David deps][david-image]][david-url]
[![npm download][download-image]][download-url]

Openchemlib extended

## DB

It allows to create a database of molecule and make substructure search

```js
const MoleculeDB = require(OCLE).DB;

let moleculeDB = new MoleculeDB();
moleculeDB.pushMoleculeInfo({ smiles: 'CCCC' });
moleculeDB.pushMoleculeInfo({ smiles: 'CCCC' });
moleculeDB.pushMoleculeInfo({ smiles: 'CCCCC' });

let result = moleculeDB.search('CC', {
  format: 'smiles',
  mode: 'substructure',
  flattenResult: false,
  keepMolecule: false
});
```

### Molecule.DB(options)

Class that allows to create a database of molecules in memory

```js
let moleculeDB = new MoleculeDB({ computeProperties: true });
```

### moleculeDB.search( query, options )

- query: molfile, smiles, idCode or a molecule (instance of OCL.Molecule)
- options:
  - format: 'smiles', 'molfile' or 'idcode'
  - mode: 'exact', 'substructure' or 'similarity'
  - flattenResult: after query should we keep the unique molecule hierarchy ?
  - keepMolecule: keep the molecule (OCL.Molecule) object in the result
  - computeProperties: predict properties like logP, logS, etc.

### parseSDF

```js
const MoleculeDB = require(OCLE).DB;
MoleculeDB.parseSDF(sdf, { computeProperties: true });
```

### parseCSV

```js
const MoleculeDB = require(OCLE).DB;
MoleculeDB.parseCSV(csv, { computeProperties: true });
```

### pushMoleculeInfo(moleculeInfo, data)

`moleculeInfo` is an object that may contain the following properties:

- idCode: the OCL unique code
- smiles
- molfile
- index (OCL substructure search index, if not present will be calculated)
- mw (used to sort search results, if not present will be calculated)

### pushEntry(molecule, data, moleculeInfo)

## RXN

Create an instance of the RXN object based on a text file in RXN format. This will also calculate for each reagent and product:

- SMILES
- Molecular formula
- Molecular weight
- idCode

```
var rxn = new RXN(rxnFile)
```

methods on rxn:

- addReagent(molfile)
- addProduct(molfile)
- toRXN() : create a new RXN file

## Molecule extension

### getGroupedDiastereotopicAtomIDs(options)

Returns an array containing all the different diastereotopic atoms in the molecule with the occurence of each of them.
It may be useful to first create the implicit hydrogens before using the method ``Molecule.addImplicitHydrogens()```.

options:

- atomLabel: filter to show only a specific atom (default: '')

### getAtomsInfo(options)

Returns an array containing many informations about the atoms.
It may be useful to first create the implicit hydrogens before using the method ``Molecule.addImplicitHydrogens()```.

### toDiastereotopicSVG(options)

Return a SVG containing a circle at the level of each atom having a data-atomid a diastereoisotopic unique identifier.

options:

- height (default: 300)
- width (default: 200)
- prefix : prefix for the id of each SVG element (default: ocl)

### getDiastereotopicHoseCodes(options)

Returns an array containing diastereotopic hoses codes

options:

- maxSphereSize: maxSphere for hose code calculation

### getGroupedHOSECodes(options)

Returns an extended groupedDiastereotopicID with hoses

options:

- atomLabel: filter to show only a specific atom (default: '')
- maxSphereSize: maxSphere for hose code calculation

### toVisualizerMolfile()

### getNumberOfAtoms()

### getExtendedDiastereotopicAtomIDs()

Returns an array of explicit hydrogens added molecule of diastereotopicAtomIDs.
Extra information like the diastereotopicIDs of the molecule are added

### getAllPaths(options)

Returns an array containing all-pairs shortest paths from a connectivity matrix using the FloydWarshall algorithm.
options:

- fromLabel: filter the pairs that start from this specific atom (default: '')
- toLabel: filter the pairs that end at this specific atom (default: '')
- minLength: min path length to report (default: 1)
- maxLength: max path length to report (default: 4)

### getConnectivityMatrix(options)

Returns an array of array (matrix) containing a '1' for all the connected atoms.
An atom is considered connected to itself (dialog is equal to 1). Options allows to
tune the values you want in the matrix.

Options:

- sdt: put the bond order outside the diagonal (default: false)
- mass: put the rounded mass on the diagonal (default: false)
- atomicNo: put the atomic number on the diagonal (default: false)

Example using npm:

```
var OCLE = require('openchemlib-extended');
var molecule = OCLE.Molecule.fromSmiles('c1ccccc1');
var matrix = molecule.getConnectivityMatrix({mass: true, sdt: true});
console.log(matrix);
```

### getMF()

Return an object containing the molecular formula and the molecular formula of each distinct
part of the molecule.
The molecular formula takes into account multihydrates, charges and isotopes.

## License

[BSD-3-Clause](./LICENSE)

[npm-image]: https://img.shields.io/npm/v/openchemlib-extended.svg?style=flat-square
[npm-url]: https://www.npmjs.com/package/openchemlib-extended
[travis-image]: https://img.shields.io/travis/cheminfo-js/openchemlib-extended/master.svg?style=flat-square
[travis-url]: https://travis-ci.org/cheminfo-js/openchemlib-extended
[david-image]: https://img.shields.io/david/cheminfo-js/openchemlib-extended.svg?style=flat-square
[david-url]: https://david-dm.org/cheminfo-js/openchemlib-extended
[download-image]: https://img.shields.io/npm/dm/openchemlib-extended.svg?style=flat-square
[download-url]: https://www.npmjs.com/package/openchemlib-extended
