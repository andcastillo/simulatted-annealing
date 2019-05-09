/**
 * Predict the chemical shifts and the coupling constants of the given sample using the cheminfo tool (does not use spinus)
 * @param {Object} sample 
 * @param {*} predictor 
 * @param {*} options 
 */
function predict1HFromSample(sample, predictor, options) {
    let molecule = sample.general.ocl.molecule;
    let h1 = predictor.proton(molecule, { OCLE: options.OCLE, group: true, ignoreLabile: false, levels: [5, 4, 3, 2, 1] });

    //Predict couplings
    let couplings = molecule.getCouplings();
    for (let coupling of couplings) {
        let fromId = coupling.fromDiaID;
        let atom = h1.filter(entry => {
            return entry.atomIDs[0] === coupling.atoms[0];
        });

        if (atom && atom.length === 1) {
            atom = atom[0];
            if (!atom.j)
                atom.j = [];

            atom.j.push({
                'assignment': [coupling.atoms[coupling.atoms.length - 1]],
                'coupling': coupling.value,
                'multiplicity': 'd',
                'distance': coupling.atoms.length - 1

            });
        }
    }

    h1.sort((a, b) => a.delta - b.delta);

    return h1;
}

module.exports = predict1HFromSample;
