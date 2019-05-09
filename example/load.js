const OCLE = require('openchemlib-extended');
const SD = require('spectra-data');
const FS = require('fs');



function loadData(filename) {
    var dataReaded = FS.readFileSync(filename).toString();
    return dataReaded;
}

/**
 * Load a molecule and spectrum from the given file names
 * @param {String} mol 
 * @param {String} jdx 
 */
function loadSample(mol, jdx) {
    let molfile = loadData(__dirname + '/' + mol);
    let jcamp = loadData(__dirname + '/' + jdx);

    let spectrum = SD.NMR.fromJcamp(jcamp, {});
    spectrum.sd.info = {};
    let molecule = OCLE.Molecule.fromMolfile(molfile);
    //let molecule = molMap.molecule;
    molecule.addImplicitHydrogens();
    let nH = molecule.getMolecularFormula().formula.replace(/.*H([0-9]+).*/, '$1') * 1;

    let sample = {
        index: 0,
        filename: jdx,
        general: { ocl: { id: molecule.getIDCode(), nH: nH, molfile: molfile, molecule: molecule } },
        spectra: {
            nmr: [{
                nucleus: 'H',
                experiment: '1d',
                //range: ranges,
                solvent: spectrum.getParamString('.SOLVENT NAME', spectrum.getSolventName()) ,
                spectrum: spectrum
            }]
        }
    };

    return sample;
}

module.exports = {loadData, loadSample};