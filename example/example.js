const SD = require("spectra-data");
const predictor = require("nmr-predictor");
const OCLE = require('openchemlib-extended');
const load = require('./load');
const autoassigner = require('./autoassigner');
const simulation = require('nmr-simulation');
const LM = require('ml-levenberg-marquardt');
const FFTlib = require("ml-fft");

const NORM = 10;

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

    //Predict the chemical shifts and coupling constants
    let prediction = predictParameters(predictor, sample);

    let spectrum = sample.spectra.nmr[0].spectrum;
    let frequency = spectrum.observeFrequencyX() * 1;

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

    console.log(ranges)
    //Guess the line width
    let lw = guessLineWidth(ranges, frequency);
    let from = ranges[0].to + 0.1; // 3.8;//spectrum.getFirstX();
    let to = ranges[ranges.length - 1].from - 0.1;//spectrum.getLastX() + 0.001;
    let smoothing = 4;


    let [dataX, dataY] = getTargetData(spectrum, from, to, smoothing);
    let nbPoints = dataX.length;


    let options1h = {
        frequency: frequency,
        from: dataX[0],
        to: dataX[nbPoints - 1],
        lineWidth: lw + smoothing,
        nbPoints: dataY.length,
        maxClusterSize: 5,
        output: 'y'
    };

    console.log(options1h)

    //Prepare the initial guess inside the spinSystem
    let spinSystem = simulation.SpinSystem.fromPrediction(prediction);
    spinSystem.ensureClusterSize(options1h);
    setGroupsFromPrediction(spinSystem, prediction);

    //-------------------Everything in Herts!!!!!!!!-------------------
    let x = new Array(nbPoints);
    for (let i = 0; i < nbPoints; i++) {
        x[i] = i;
    }
    // array of points to fit
    let data = { x: x, y: dataY };

    //Step 1. Fit a wide peaks
    let fittedParams = runFitting(data, spinSystem, smoothing, options1h, 0.1, 100);
    console.log("Error 1: " + fittedParams.parameterError.toFixed(2))
    //API.createData("error", fittedParams.parameterError.toFixed(2));

    //Step 2. Fit the raw data
    smoothing = 0;
    let options1h2 = {
        frequency: frequency,
        from: dataX[0],
        to: dataX[nbPoints - 1],
        lineWidth: lw + smoothing,
        nbPoints: dataY.length,
        maxClusterSize: 5,
        output: 'y'
    };
    //options1h.lineWidth = lw + smoothing;
    let spinSystem2 = setParamsOnSpinSystem(fittedParams.parameterValues, spinSystem);
    //console.log(spinSystem2);
    [dataX, dataY] = getTargetData(spectrum, from, to, smoothing);
    data = { x: x, y: dataY };
    fittedParams = runFitting(data, spinSystem2, smoothing, options1h2, 0.01, 100);


    //Display the spectrum,
    //var simfx = simulatedSpectrum(fittedParams.parameterValues, spinSystem2, options1h);
    //let ySim = x.map(t => simfx(t))

    //API.createData("sp0", SD.NMR.fromXY(dataX, ySim, { nucleus: "1H" }).sd);
    //API.createData("error", fittedParams.parameterError.toFixed(2));
    console.log("Error 2: " + fittedParams.parameterError.toFixed(2))

    return 0;
}

/**
 * This function obtains the chemical shifts and the coupling constants for a given molecule
 * @param {Object} predictor 
 * @param {*} sample 
 */
function predictParameters(predictor, sample) {
    let molecule = sample.general.ocl.molecule;
    //predictor.fetchProton(url: "https://raw.githubusercontent.com/cheminfo-js/spectra/master/packages/nmr-predictor/data/nmrshiftdb2-1h.json");
    let h1 = predictor.proton(molecule, { OCLE, group: true, ignoreLabile: false, levels: [5, 4, 3, 2, 1] });

    //Predict couplings
    let couplings = molecule.getCouplings();
    for (let coupling of couplings) {
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


/**
 * Fit the given spin system to the dataY array of points 
 */
function runFitting(data, spinSystem2, smoothing, options1h, gradientDifference, maxIterations) {
    let initialValues = []; //spinSystem.csGroups.map(value => value.delta);
    for (let key in spinSystem2.csGroups) {
        initialValues.push(spinSystem2.csGroups[key].delta * 200);
    }
    let minValues = initialValues.map(value => value - 10);
    let maxValues = initialValues.map(value => value + 10);

    let jc = spinSystem2.jGroups;
    for (let key in jc) {
        initialValues.push(jc[key].value);
        minValues.push(jc[key].value * 0.5);
        maxValues.push(jc[key].value * 1.5);
    }

    let options = {
        damping: 0.01,
        initialValues,
        minValues,
        maxValues,
        gradientDifference: gradientDifference,
        maxIterations: maxIterations,
        errorTolerance: 10e-3
    };

    /** function that receives the parameters and returns
     * a function with the independent variable as a parameter
     * @param {Array} params 
     */
    let fx = function (params) {
        //Set the chemical shits taking care of keeping the same chemical shift for each equivalente proton
        return simulatedSpectrum(params, spinSystem2, options1h);
    }

    console.log("Initial");
    console.log(initialValues);
    //console.log(options);


    let fittedParams = LM(data, fx, options);

    //Simulate the resulting spectrum
    console.log("Step finished!!!!");
    console.log(fittedParams);

    return fittedParams;
}

function simulatedSpectrum(params, sp, options1h) {
    //Set the chemical shits taking care of keeping the same chemical shift for each equivalente proton
    sp = setParamsOnSpinSystem(params, sp);
    let spectrum = simulation.simulate1D(sp, options1h);
    //We must normalize the spectrum before compare
    let sum = 0;
    for (let i = 0; i < spectrum.length; i++) {
        sum += spectrum[i];
    }
    sum = NORM / sum;
    for (let i = 0; i < spectrum.length; i++) {
        spectrum[i] *= sum;
    }
    return (t) => spectrum[t];
}

function setParamsOnSpinSystem(params, sp) {
    //Set the chemical shits taking care of keeping the same chemical shift for each equivalente proton
    let index = 0;
    for (let key in sp.csGroups) {
        let group = sp.csGroups[key].indexes;
        let value = params[index] / 200;
        if (!value)
            console.log(key + " ?? " + value + " " + index);
        sp.csGroups[key].delta = value;
        for (let j = 0; j < group.length; j++) {
            sp.chemicalShifts[group[j]] = value;
        }
        index++;
    }

    let jGroups = sp.jGroups;
    let jc = sp.couplingConstants;
    for (let key in jGroups) {
        let j = jGroups[key];
        let indexes = j.indexes;
        let value = params[index];
        j.value = value;
        if (!value)
            console.log(key + " value " + value + " " + index);
        for (let k = 0; k < indexes.length; k++) {
            jc[indexes[k][0]][indexes[k][1]] = value;
            jc[indexes[k][1]][indexes[k][0]] = value;
        }
        index++;
    }

    return sp;
}

/**
 * Guess the linewidth of the spectrum from the set peaks contained in the ranges
 * 
 */
function guessLineWidth(ranges, frequency) {
    let widths = ranges.reduce((res1, value) =>
        value.signal.reduce((res2, signal) =>
            signal.peak.reduce((res3, peak) => {
                res3.push(peak.width);
                return res3;
            }, res2), res1), []);

    return frequency * (widths.reduce((sum, v) => sum + v) / widths.length) / 1.1775;
}

/**
 * Get the data to be fitted
 */
function getTargetData(spectrum, from, to, smoothing) {
    let spectrum2 = spectrum;

    if (smoothing > 0) {
        let re = spectrum.getYData(0);
        let im = spectrum.getYData(1);
        let re0 = re.slice(); //Original real data 
        let im0 = im.slice();
        smoothSpectrum(re0, im0, smoothing * 20);
        spectrum2 = SD.NMR.fromXY(spectrum.getXData(0), re0, { nucleus: "1H" });
    }

    //API.createData("spectrum2", spectrum2.sd);

    //We must normalize the spectrum before run this optimization
    //let factor = spectrum.getNbPoints() / spectrum.getArea(spectrum.getFirstX(), spectrum.getLastX());
    let xy = spectrum2.getPointsInWindow(from, to, { outputX: true });
    let yy = xy.y;
    let xx = xy.x;
    let nbPoints = yy.length;
    let sum = 0;
    for (let i = 0; i < nbPoints; i++) {
        sum += yy[i];
    }
    sum = NORM / sum;
    let dataY = new Array(nbPoints);
    let dataX = new Array(nbPoints);
    for (let i = 0; i < nbPoints; i++) {
        dataY[i] = yy[nbPoints - 1 - i] * sum;
        dataX[i] = xx[nbPoints - 1 - i];
    }

    return [dataX, dataY];
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

/**
 *In place smoothing of the given real and imaginary arrays 
 */
function smoothSpectrum(re, im, std) {
    let nbPoints = re.length;
    FFTlib.FFT.init(nbPoints);
    FFTlib.FFT.fft(re, im);
    for (let i = 0; i < re.length; i++) {
        let xi = std * (re.length - i) / re.length;
        //let g = Math.exp(-0.5*(Math.pow(xi - mean) / std, 2)); //Gaussian
        let fn = Math.exp(-0.5 * xi);//std * 2 / (Math.pow(xi - mean, 2) + 0.25 * std * std * 4) * g + 0.01; 
        re[i] = re[i] * fn;
        im[i] = im[i] * fn;
    }
    //API.createData("apodization", fn);

    FFTlib.FFT.ifft(re, im);
}

