const functions = require("../src/functions");
const FS = require('fs');
const SD = require("spectra-data");
const nmrPredictor = require("nmr-predictor");
const simulatedAnnealing = require("../src/simulatedAnnealing");

const {
    couplingIndexesGenerator,
    rangesGenerator,
    initialParametersValuesGenerator,
    couplingConstantsSimilarity,
    fillRangeUpToDown
} = functions;

let molfile = createSpectraData(__dirname + '/../data/mol_0.mol');
let jdxData = createSpectraData(__dirname + '/../data/h1_0.jdx');

runSA(jdxData, molfile);

async function runSA(jdxData, molfile) {
    let prediction = await nmrPredictor.spinus(molfile, {group: true});
    const stepLength = 0.2;
    const minimumSpectraFrequency = 0;
    const maximumSpectraFrequency = 9;
    const frequencyRange = maximumSpectraFrequency - minimumSpectraFrequency;
    const spectraProperties = {
        frequency: 400.082470657773,//MHz
        from: minimumSpectraFrequency,//PPM
        to: maximumSpectraFrequency,//PPM
        lineWidth: 0.7,//Hz
        nbPoints: 5*1024 ,
        maxClusterSize: 8,
        output:"xy"
    };
    const rho = 0.5;
    const peaksOptions = {//important to the peak peaking****debug whit diferent spectra
        thresholdFactor: 1,
        optimize: true,
        minMaxRatio: 0.1,
        broadRatio: 0.00025,
        smoothY: true,
        widthFactor: 4,
        realTop: true,
        functionName: 'gaussian',
        broadWidth: 1,
        sgOptions: {windowSize: 9, polynomial: 3}
    };
    let experimentalData = SD.NMR.fromJcamp(jdxData);
    //experimentalData.reduceData(0, 10, spectraProperties);
    experimentalData.setMinMax(0, 1);

    let solvent = experimentalData.getSolventName();
    let solventShift =  experimentalData.getResidual(solvent)[0].shift;
    if (solvent === 'DMSO'){
        peaksOptions.minMaxRatio = 0.25;
        experimentalData.fill(3.1, 4, 0);
    }

    const experimentalSpectra = {
        'x': experimentalData.getXData(),
        'y': experimentalData.getYData()
    };

    let experimentalSignals = experimentalData.getRanges(peaksOptions).map(x => x.signal[0].delta);
    experimentalSignals = experimentalSignals.sort();
    experimentalSignals = mergeSort(experimentalSignals);

    const sortDelta = (prediction.map(x => x.delta)).sort();
    const sortPrediction =  (JSON.parse(JSON.stringify(prediction))).sort((a, b) => a.delta - b.delta);
    const range = {
        'min' : 0,
        'max' : sortDelta.length - 1,
        'signals': sortDelta.length
    };

    let condition;
    let psi = 0.15;

    if(sortPrediction.length === experimentalSignals.length) {
        for (let u = 0; u < sortPrediction.length; u++) {
            sortPrediction[u].delta = experimentalSignals[u];
        }
        psi = 0.01;
    }

    const gama = 1;
    const spinSystem = nmr.SpinSystem.fromPrediction(sortPrediction);// Create the spin system of the drawed molecule;
    const alfa = psi;
    const beta = 2;

    const couplingIndexes = couplingIndexesGenerator(sortPrediction);
    const [sortPredictionRanges, regionsRanges] = rangesGenerator(0, 10, 5);

    let originalParameters = initialParametersValuesGenerator(sortPrediction);
    let originalCouplingConstants = originalParameters.slice(sortPrediction.length);

    const [couplingConstantsSummaryIndexes, couplingsSummary] = couplingConstantsSimilarity(originalCouplingConstants);
    let sortCouplingsValues = JSON.parse(JSON.stringify(couplingsSummary));
    sortCouplingsValues = mergeSort(sortCouplingsValues);
    let initialParameters = initialParametersGenerator(sortPrediction);

    let regions = new Matrix(regionsRanges.length, 2);
    for (let i = 0; i < regionsRanges.length; i++) {
        let minimumTmpValue = regionsRanges[i][0] - gama;
        let maximumTmpValue = regionsRanges[i][regionsRanges[i].length - 1] + gama;
        if (minimumTmpValue < minimumSpectraFrequency) { minimumTmpValue = minimumSpectraFrequency}
        if (maximumTmpValue > maximumSpectraFrequency) {maximumTmpValue = maximumSpectraFrequency }

        regions[i][0] = minimumTmpValue;
        regions[i][1] = maximumTmpValue;
    }
    let minWidth = sortCouplingsValues[0] - sortCouplingsValues[0] * rho;
    let maxWidth = sortCouplingsValues[sortCouplingsValues.length - 1] + sortCouplingsValues[sortCouplingsValues.length - 1] * rho;
    let widthsRange = maxWidth - minWidth;
    let widths = fillRangeUpToDown(minWidth, maxWidth, stepLength);
    let widthsCopy = widths.map(x => x / spectraProperties.frequency);
    let maxRate = 0.5;
    let minRate = 0.01;
    let step = (maxRate - minRate) / widths.length;


    const boundsRangesLength = initialParameters[2].clone().subtract(initialParameters[1]);
//****************************************************************
    /*regions = new Matrix(regionsRanges.length, 2);
    for (let i = 0; i < regionsRanges.length; i++) {
        let minimumTmpValue = regionsRanges[i][0] - gama;
        let maximumTmpValue = regionsRanges[i][regionsRanges[i].length - 1] + gama;
        if (minimumTmpValue < minimumSpectraFrequency) { minimumTmpValue = minimumSpectraFrequency}
        if (maximumTmpValue > maximumSpectraFrequency) {maximumTmpValue = maximumSpectraFrequency }

        regions[i][0] = minimumTmpValue;
        regions[i][1] = maximumTmpValue;
    };
    minWidth = sortCouplingsValues[0] - sortCouplingsValues[0] * rho;
    maxWidth = sortCouplingsValues[sortCouplingsValues.length - 1] + sortCouplingsValues[sortCouplingsValues.length - 1] * rho;
    widthsRange = maxWidth - minWidth;
    widths = fillRangeUpToDown(minWidth,maxWidth, stepLength);
    widthsCopy = widths.map(x => x / spectraProperties.frequency);
    maxRate = 0.5;
    minRate = 0.01;
    step = (maxRate - minRate) / widths.length;
    */
    const boundsRangesLengthDecreasing = fillRangeUpToDown(minRate, maxRate, step);

    const inputs = {
        objetiveFunction: delta,
        guess :initialParameters[0],
        lowerBound :initialParameters[1],
        upperBound :initialParameters[2],
        iterationsNumber :5000,
        quenchingFactor : 1,
        toleranceValue : 1E-08
    };

    const baseMatrix = {};
    const indexesBoundsMatrix = {};
    let simulatedAnnealingOutputs = [];

    for (let i = 0; i < widths.length; i++) {
        let centralFrequency = [];
        let points = Math.round(gama * 4 * spectraProperties.frequency / widths[i]);
        for (let j = 0; j < regions.length; j++) {
            centralFrequency = centralFrequency.concat(fillRangeWithoutBounds(regions[j][0], regions[j][1], points));
        }

        let base = [];
        let boundsIndexes = [];
        let tmpOutput =[];

        for(let t = 0; t < centralFrequency.length; t++) {
            tmpOutput[t] = transformationBase(experimentalSpectra.x, centralFrequency[t], widthsCopy[i]);
            base[t] = tmpOutput[t][0];
            boundsIndexes[t] = tmpOutput[t][1];
        }
        baseMatrix.base = base;
        indexesBoundsMatrix.boundsIndexes = boundsIndexes;

        simulatedAnnealingOutputs[i] = simulatedAnnealing(inputs);
        // if ( i > 0 &&  simulatedAnnealingOutputs[i][1] - simulatedAnnealingOutputs[i - 1][1] >= 0) {
        //  inputs.guess = simulatedAnnealingOutputs[i - 1][0];
        // }

        // else {inputs.guess = simulatedAnnealingOutputs[i][0]}
        inputs.guess = simulatedAnnealingOutputs[i][0];// review
        // inputs.lowerBound = Matrix.rowVector(subtract(inputs.guess[0], multiply(boundsRangesLength[0], boundsRangesLengthDecreasing[i])));
        // inputs.upperBound = Matrix.rowVector(subtract(inputs.guess[0], multiply(boundsRangesLength[0], boundsRangesLengthDecreasing[i])));
    };

    let finalParameters = simulatedAnnealingOutputs[simulatedAnnealingOutputs.length - 1][0];

    let transformPrube = transformationBase(experimentalSpectra.x, 8, 2);
    let tPrube = {
        'x' : experimentalSpectra.x,
        'y' : transformPrube
    }
// API.createData('tPrube', tPrube);
    let prube = fillRangeUpToDown(10, 50, 1);
    let prube0 = fillRangeWithoutBounds(10, 50, 10);
    let generalFitSpectra = spectraDesigner(finalParameters);
};


// // API.createData('generalFitSpectra', generalFitSpectra);
// // API.createData('experimentalSpectra', experimentalSpectra);

function createSpectraData(filename) {
    var dataReaded = FS.readFileSync(filename).toString();
    return dataReaded;
}