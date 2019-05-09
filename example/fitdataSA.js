const functions = require("../src/functions");
const FS = require('fs');
const SD = require("spectra-data");
const predictor = require("nmr-predictor");
const OCLE = require('openchemlib-extended');
const load = require('./load');
const autoassigner = require('./autoassigner');
const simulation = require('nmr-simulation');
const SA = require('../src/index');
const predictFromSample = require('./predictor');


run('/../data/mol_0.mol', '/../data/h1_0.jdx');

async function run(mol, jdx) {
    //Load the sample from the files
    let sample = load.loadSample(mol, jdx);
    //Load  and set the db for chemical shift prediction
    let db = JSON.parse(load.loadData(__dirname + '/../data/nmrshiftdb2-1h.json'));
    predictor.setDb(db, 'proton', 'proton');

    //Run the autoassigner for the given sample. It will use the chemical shift predictions to improve the result
    let result = await autoassigner(sample, predictor, { OCLE });
    let ranges = result.ranges;
    let assigner = result.assigner;
    //let [ranges, assigner] = await autoassigner(sample, predictor, { OCLE });

    //Set the best (index 0) assignment on the peak picking of this sample
    assigner.setAssignmentOnRanges(ranges, 0);

    //Predict the chemical shift and coupling constants for the sample
    let prediction = predictFromSample(sample, predictor, { OCLE });

    /*************Fitting preparation starts here******************** */

    let spectrum = sample.spectra.nmr[0].spectrum;
    let frequency = spectrum.observeFrequencyX() * 1;

    let from = spectrum.getFirstX();
    let to = spectrum.getLastX() + 0.001;
    //We must normalize the spectrum before run this optimization
    //let factor = spectrum.getNbPoints() / spectrum.getArea(spectrum.getFirstX(), spectrum.getLastX());
    let xy = spectrum.getPointsInWindow(from, to, { outputX: true });
    let yy = xy.y;
    let xx = xy.x;
    let nbPoints = yy.length;
    let sum = 0;
    for (let i = 0; i < nbPoints; i++) {
        sum += yy[i];
    }
    console.log("nbPoints " + nbPoints)
    sum = 1 / sum;
    let dataY = new Array(nbPoints);
    let dataX = new Array(nbPoints);
    for (let i = 0; i < nbPoints; i++) {
        dataY[i] = yy[nbPoints - 1 - i] * sum;
        dataX[i] = xx[nbPoints - 1 - i];
    }


    //Copy the parameters from the assignment into the prediction
    prediction.forEach(value => {
        let oclID = value.diaIDs[0];
        for (let i = 0; i < ranges.length; i++) {
            let diaIDs = ranges[i].signal[0].diaID;
            for (let id = 0; id < diaIDs.length; id++) {
                if (diaIDs[id] + "" === oclID) {
                    value.delta = ranges[i].signal[0].delta;
                    break;
                }
            }
        }
    });

    //console.log(frequency)

    var options1h = {
        frequency: frequency,
        from: dataX[0],
        to: dataX[nbPoints - 1],
        lineWidth: 0.8,
        nbPoints: dataY.length,
        maxClusterSize: 5,
        output: 'y'
    };

    //Create a new spin system that will be used during the fitting process
    var spinSystem = simulation.SpinSystem.fromPrediction(prediction);
    //We create the cluster only once. Each simulation must use the same clustering afterwards.
    spinSystem.ensureClusterSize(options1h);
    //Move this function to SpinSystem
    setGroupsFromPrediction(spinSystem, prediction)

    // function that receives the parameters and returns
    // a function with the independent variable as a parameter
    function simulatedSpectrum(params) {
        //Set the chemical shits taking care of keeping the same chemical shift for each equivalente proton
        let index = 0;
        for (let key in spinSystem.csGroups) {
            let group = spinSystem.csGroups[key].indexes;
            for (let j = 0; j < group.length; j++) {
                spinSystem.chemicalShifts[group[j]] = params[index];
            }
            index++;
        }

        let jGroups = spinSystem.jGroups;
        let jc = spinSystem.couplingConstants;
        for (let key in jGroups) {
            let j = jGroups[key];
            let indexes = j.indexes;
            for (let k = 0; k < indexes.length; k++) {
                jc[indexes[k][0]][indexes[k][1]] = params[index];
                jc[indexes[k][1]][indexes[k][0]] = params[index];
            }
            index++;
        }

        var spectrum = simulation.simulate1D(spinSystem, options1h);
        //We must normalize the spectrum before compare
        let sum = 0;
        for (let i = 0; i < spectrum.length; i++) {
            sum += spectrum[i];
        }
        sum = 1 / sum;
        for (let i = 0; i < spectrum.length; i++) {
            spectrum[i] *= sum;
        }

        return (t) => spectrum[t];
    }

    //-------------------Everything in Herts!!!!!!!!-------------------
    let x = new Array(options1h.nbPoints);
    for (let i = 0; i < options1h.nbPoints; i++) {
        x[i] = i;
    }
    // array of points to fit
    let data = {
        x: x,
        y: dataY
    };

    //console.log(JSON.stringify(prediction))
    //console.log(spinSystem.couplingConstants);
    // array of initial parameter values
    let initialValues = []; //spinSystem.csGroups.map(value => value.delta);
    for (let key in spinSystem.csGroups) {
        initialValues.push(spinSystem.csGroups[key].delta);
    }
    let minValues = initialValues.map(value => value - 0.05);
    let maxValues = initialValues.map(value => value + 0.05);

    let jc = spinSystem.jGroups;
    for (let key in jc) {
        initialValues.push(jc[key].value);
        minValues.push(jc[key].value * 0.8);
        maxValues.push(jc[key].value * 1.0);
    }

    console.log(initialValues);
    console.log("Initial");
    //console.log(minValues);
    //console.log(maxValues);

    const options = {
        initialValues,
        minValues,
        maxValues,
        quenchingFactor: 1,
        maxIterations: 10000,
        errorTolerance: 10e-3
    };

    let fittedParams = SA(data, simulatedSpectrum, options);


    //Simulate the resulting spectrum
    console.log("Simulation finished!!!!");
    console.log(fittedParams);

    var simfx = simulatedSpectrum(fittedParams.parameterValues);
    let ySim = x.map(t => simfx(t))

    //console.log(SD.NMR.fromXY(dataX, ySim, {nucleus: "1H"}).toJcamp());
}

/**
 * Add the symmetry groups to the spin system
 * @param {Object} spinSystem 
 * @param {Array} prediction 
 */
function setGroupsFromPrediction(spinSystem, prediction) {
    let cs = {};
    let atomMap = {};
    let atomIndex = {};
    //Look up for the chemical shifts groups
    for (let i = 0; i < prediction.length; i++) {
        let spin = prediction[i];
        atomMap[spin.atomIDs[0]] = spin.diaIDs[0];
        atomIndex[spin.atomIDs[0]] = i;
        if (cs[spin.diaIDs[0]]) {
            cs[spin.diaIDs[0]].atoms.push(spin.atomIDs[0]);
            cs[spin.diaIDs[0]].indexes.push(i);

        } else {
            cs[spin.diaIDs[0]] = { atoms: [spin.atomIDs[0]], indexes: [i], delta: spin.delta };
        }
    }
    //Look up for the coupling constants groups
    let jg = {};
    for (let i = 0; i < prediction.length; i++) {
        let spin = prediction[i];
        if (spin.j) {
            let j = spin.j;
            for (let k = 0; k < j.length; k++) {
                //The coupling constan key is a composition of the oclID of the atoms + the path length among them.
                let key = spin.diaIDs[0];
                let jIndex = atomIndex[j[k].assignment];
                if (spin.diaIDs[0] < atomMap[j[k].assignment[0]]) {
                    key += '-' + atomMap[j[k].assignment[0]] + '-' + j[k].distance;
                } else {
                    key = atomMap[j[k].assignment[0]] + '-' + key + '-' + j[k].distance;
                }
                if (!jg[key]) {
                    jg[key] = { value: j[k].coupling, indexes: [[i, jIndex]] };

                } else {
                    jg[key].indexes.push([i, jIndex]);
                }
            }
        }
    }

    spinSystem.csGroups = cs;
    spinSystem.jGroups = jg;
}