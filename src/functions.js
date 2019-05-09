const Matrix = require("ml-matrix");
const SD = require("spectra-data");

module.exports.couplingIndexesGenerator = function(sortPrediction) {
    let couplingMatrixIndexes = [];
    let counter = 0;
    for (let i = 0; i < sortPrediction.length; i++){
        let iteration = 0;
        let tmp = [];
        for (let j = 0; j < sortPrediction[i].j.length; j++){
            tmp[iteration] = counter;
            counter++;
            iteration++;
        }
        couplingMatrixIndexes[i] = tmp;
    }
    return couplingMatrixIndexes;
};

module.exports.rangesGenerator = function(min, max, rangesNumber){
    let rangesLength = Math.round((max - min) / rangesNumber);
    let tmp = [];
    for(let i = 0; i < rangesNumber; i++){
        tmp[i] = sortDelta.filter(x => x > (rangesLength * i) & x < (rangesLength * (i + 1)));
    }
    let chemicalShiftsGroups = tmp.filter(x => x.length > 0);
    let sortPredictionRanges = new Matrix(chemicalShiftsGroups.length,1);
    let counter = 0;
    for(let i = 0; i < chemicalShiftsGroups.length; i++ ){
        for(let j = 0; j < chemicalShiftsGroups[i].length; j++){
            sortPredictionRanges[i][j] = counter;
            counter++;
        }
    }
    return [sortPredictionRanges, chemicalShiftsGroups];
};

module.exports.initialParametersValuesGenerator = function(sortPrediction) {
    let initialPrediction = JSON.parse(JSON.stringify(sortPrediction));
    let initialDelta = [];
    let iterations = 0;
    for (let i = range.min; i < range.max + 1; i++){
        initialDelta[iterations] = initialPrediction[i].delta;
        iterations++;
    }
    let initialCouplings = new Matrix(range.signals, 1);// esta parte se puede eliminar
    let counter = 0;
    for(let i = range.min; i < range.max + 1; i++){
        for(let j = 0; j < initialPrediction[i].j.length; j++){
            initialCouplings[counter][j]  = initialPrediction[i].j[j].coupling;
        }
        counter++;
    }
    let ungroupCouplings = [];
    for(let i = 0; i < initialCouplings.length; i++) {
        ungroupCouplings = ungroupCouplings.concat(initialCouplings[i]);
    }
    let initialSpectraParameters = initialDelta .concat(ungroupCouplings);
    return initialSpectraParameters;
};

module.exports.couplingConstantsSimilarity = function(couplingConstants) {
    let couplingConstantsSummaryIndexes = [];
    let couplingConstantsSummary = [];
    let indexesValues = [];
    let uniqueIndexesValues = [];
    for (let i = 0; i < couplingConstants.length; i++) {
        let tmp = [];
        let counter = 0;
        for (let j = i ; j < couplingConstants.length; j++) {
            if (couplingConstants[i] === couplingConstants[j]) {
                tmp[counter] = j;
                counter++
            }
            couplingConstantsSummaryIndexes[i] = tmp
        }
    }
    couplingConstantsSummaryIndexes = couplingConstantsSummaryIndexes.filter(x => x.length > 1 );
    for (let i = 0; i < couplingConstantsSummaryIndexes.length;  i++){
        couplingConstantsSummary[i] = couplingConstants[couplingConstantsSummaryIndexes[i][0]];
    }
    for (let i = 0; i < couplingConstantsSummary.length; i++) {
        for (let j = i + 1; j < couplingConstantsSummary.length; j++ ){
            if (couplingConstantsSummary[i] === couplingConstantsSummary[j])
                couplingConstantsSummary[j] = undefined;
        }
    }
    for (let i = 0; i < couplingConstantsSummaryIndexes.length; i++){
        if (couplingConstantsSummary[i] === undefined) {
            couplingConstantsSummaryIndexes[i] = undefined;
        }
    }
    couplingConstantsSummaryIndexes = couplingConstantsSummaryIndexes.filter(x => x !== undefined );
    couplingConstantsSummary = couplingConstantsSummary.filter(x => x !== undefined );
    let counter = 0;
    for (let i = 0; i < couplingConstantsSummaryIndexes.length; i++) {
        for (let j = 0; j < couplingConstantsSummaryIndexes[i].length; j++ ) {
            indexesValues[counter] = couplingConstantsSummaryIndexes[i][j];
            counter++
        }
    }
    let sortIndexesValues =  mergeSort(indexesValues);
    counter = 0;
    for (let i = 0; i < couplingConstants.length; i++) {
        sortIndexesValues.some(x => x === i) ? uniqueIndexesValues[i] = undefined : uniqueIndexesValues[i] = i ;
    }
    uniqueIndexesValues = uniqueIndexesValues.filter(x => x !== undefined);
    let tmp = [];
    for (let i = 0; i < uniqueIndexesValues.length; i++){
        tmp[i] = couplingConstants[uniqueIndexesValues[i]];
    }
    uniqueIndexesValues = Matrix.columnVector(uniqueIndexesValues)
    couplingConstantsSummary = couplingConstantsSummary.concat(tmp);
    couplingConstantsSummary = couplingConstantsSummary.filter(x => x !== undefined);
    couplingConstantsSummaryIndexes = couplingConstantsSummaryIndexes.concat(uniqueIndexesValues);

    return [couplingConstantsSummaryIndexes, couplingConstantsSummary] ;
};

module.exports.spectraDesigner = function(spectraParameters) {
    let predictionCopy = JSON.parse(JSON.stringify(sortPrediction));
    let chemicalShiftsValues = spectraParameters[0].slice(0,range.signals);
    let iteration = 0;
    for(let i = range.min; i < range.max + 1 ; i++) {
        predictionCopy[i].delta = chemicalShiftsValues[iteration];
        iteration++;
    }
    let couplingConstantsValues = spectraParameters[0].slice(range.signals);
    let ungroupCouplings = [];
    for (let i = 0; i < couplingConstantsValues.length; i++) {
        for (let j = 0; j < couplingConstantsSummaryIndexes[i].length; j++) {
            ungroupCouplings[couplingConstantsSummaryIndexes[i][j]] = couplingConstantsValues[i];
        }
    }
    for(let i = range.min; i < range.max + 1; i++){
        for(let j = 0; j < predictionCopy[i].j.length; j++){
            predictionCopy[i].j[j].coupling = ungroupCouplings[couplingIndexes[i][j]];
        }
    }
    let  simulationData = SD.NMR.fromSignals(predictionCopy, spectraProperties);
    simulationData.setMinMax(0, 1);
    // let element = simulationData.getActiveElement();
    // simulationData.suppressZone(suppressOptions.form, suppressOptions.to);
    // simulationData.setActiveElement(element);
    let  spectra = {
        'x': simulationData.getXData(),
        'y': simulationData.getYData()
    };
    return spectra;
}; // Function to make spectra from spectra Properties or spectra parameters(j and chemical shifts).

module.exports.merge = function(left, rigth) {
    let result = [];
    let i = 0;
    let j = 0;
    while (i < left.length && j < rigth.length) {
        if (left[i] < rigth[j]) {
            result.push(left[i++]);
        } else {
            result.push(rigth[j++]);
        }
    }
    return result.concat(left.slice(i)).concat(rigth.slice(j));
};

module.exports.mergeSort = function(numbers) {
    if (numbers.length < 2) {
        return numbers;
    }
    let mid = Math.floor(numbers.length/2);
    let left = numbers.slice(0,mid);
    let rigth = numbers.slice(mid);

    return merge(mergeSort(left), mergeSort(rigth));
};

module.exports.initialParametersGenerator = function(sortPrediction) {
    let initialPrediction = JSON.parse(JSON.stringify(sortPrediction));
    let ChemicalShiftsInitialError =  alfa ;
    let CouplingConstantsInitialError = beta;
    let initialDelta = [];
    let iteration = 0;
    for (let i = range.min; i < range.max + 1; i++){
        initialDelta[iteration] = initialPrediction[i].delta;
        iteration++;
    }
    let chemicalShiftsLowerBound = initialDelta.map(x => x - ChemicalShiftsInitialError);
    let chemicalShiftsUpperBound = initialDelta.map(x => x + ChemicalShiftsInitialError);
    let couplingConstantsLowerBound = couplingsSummary.map(x => x - CouplingConstantsInitialError);
    let couplingConstantsUpperBound = couplingsSummary.map(x => x + CouplingConstantsInitialError);
    let initialSpectraParameters =  Matrix.rowVector(initialDelta.concat(couplingsSummary));
    let lowerBound = Matrix.rowVector(chemicalShiftsLowerBound.concat(couplingConstantsLowerBound));
    let upperBound = Matrix.rowVector(chemicalShiftsUpperBound.concat(couplingConstantsUpperBound));
    return [initialSpectraParameters,lowerBound, upperBound];
};// Function to generate the guess values.

module.exports.bounds = function(spectraParameters, index){
    let ChemicalShiftsInitialError = alfa * boundsRangesLengthDecreasing[index];
    let CouplingConstantsInitialError = beta * boundsRangesLengthDecreasing[index];
    let chemicalShiftsValues = spectraParameters[0].slice(0,range.signals);
    let couplingConstantsValues = spectraParameters[0].slice(range.signals);
    let chemicalShiftsLowerBound = chemicalShiftsValues.map(x => x - ChemicalShiftsInitialError );
    let chemicalShiftsUpperBound = chemicalShiftsValues.map(x => x + ChemicalShiftsInitialError );
    let couplingConstantsLowerBound = couplingConstantsValues.map(x => x - CouplingConstantsInitialError);
    let couplingConstantsUpperBound = couplingConstantsValues.map(x => x + CouplingConstantsInitialError);
    let lowerBound = Matrix.rowVector(chemicalShiftsLowerBound.concat(couplingConstantsLowerBound));
    let upperBound = Matrix.rowVector(chemicalShiftsUpperBound.concat(couplingConstantsUpperBound));
    return [lowerBound, upperBound]
}

module.exports.transformationBase = function(frequencies, centralFrequency, width) {
    let lowerBoundIndexValue = frequencies.indexOf(centralFrequency - (width / 2));
    let upperBoundIndexValue = frequencies.indexOf(centralFrequency + (width / 2));
    let indexesBoundsValues = [lowerBoundIndexValue, upperBoundIndexValue];
    let yMax = 1.0;
    let newY = [];
    let indexCounter = 0;
    for (let index = lowerBoundIndexValue ; index < upperBoundIndexValue + 1; index++) {
        if (frequencies[index] < centralFrequency) {
            newY[indexCounter] = (2 * yMax / width) * (frequencies[index]- centralFrequency + width / 2);
            indexCounter++
        }
        else if (frequencies[index] > centralFrequency) {
            newY[indexCounter] = (-2 * yMax / width) * (frequencies[index] - centralFrequency- width / 2);
            indexCounter++
        }

        else if (frequencies[index] === centralFrequency) {
            newY[indexCounter] = yMax;
            indexCounter++
        }
    }
    return [newY, indexesBoundsValues];
};

module.exports.delta = function(spectraParameters){
    let expSpectra = Matrix.columnVector(experimentalSpectra.y);
    let variableTheoreticalSpectra = Matrix.columnVector(spectraDesigner(spectraParameters).y);
    let delta = 0;
    for( let h = 0; h < baseMatrix.base.length; h++ ){
        let lBound = indexesBoundsMatrix.boundsIndexes[h][0];
        let uBound = indexesBoundsMatrix.boundsIndexes[h][1];
        let integralOfExperimentalSpectra =  collect(multiply(baseMatrix.base[h],expSpectra.slice(lBound, uBound + 1 )));
        let integralOfTheoreticalSpectra  = collect(multiply(baseMatrix.base[h],variableTheoreticalSpectra.slice(lBound, uBound + 1 )));
        let differences = integralOfExperimentalSpectra - integralOfTheoreticalSpectra;
        delta = delta + Math.pow(differences, 2);
    }
    return delta;
};

module.exports.sum = function(parameterOne,parameterTwo){
    if(typeof parameterOne !== 'number' && typeof parameterTwo !== 'number' && parameterOne.length !== parameterTwo.length) {
        return 'parameterOne and parameterTwo should have the same length';
    }
    else {
        let result = [];
        if(typeof parameterOne === "number" && typeof parameterTwo === "number" ){
            result = 'is unnecesary use this function, both parameters are numbers';
        }
        else if (typeof parameterOne !== "number" && typeof parameterTwo !== "number" ){
            for(let i = 0; i < parameterOne.length; i++){
                result[i] = parameterOne[i] + parameterTwo[i];
            }
        }
        else if (typeof parameterOne !== "number" && typeof parameterTwo === "number" ){
            for(let i = 0; i < parameterOne.length; i++){
                result[i] = parameterOne[i] + parameterTwo;
            }
        }
        else if (typeof parameterOne === "number" && typeof parameterTwo !== "number" ){
            for(let i = 0; i < parameterTwo.length; i++){
                result[i] = parameterTwo[i] + parameterOne;
            }
        }
        return result;
    }
};

module.exports.subtract = function(parameterOne,parameterTwo){
    if(typeof parameterOne !== 'number' && typeof parameterTwo !== 'number' && parameterOne.length !== parameterTwo.length) {
        return 'parameterOne and parameterTwo should have the same length';
    }
    else {
        let result = [];
        if(typeof parameterOne === "number" && typeof parameterTwo == "number" ){
            result = 'is unnecesary use this function, both parameters are numbers';
        }
        else if (typeof parameterOne !== "number" && typeof parameterTwo !== "number" ){
            for(let i = 0; i < parameterOne.length; i++){
                result[i] = parameterOne[i] - parameterTwo[i];
            }
        }
        else if (typeof parameterOne !== "number" && typeof parameterTwo === "number" ){
            for(let i = 0; i < parameterOne.length; i++){
                result[i] = parameterOne[i] - parameterTwo;
            }
        }
        else if (typeof parameterOne === "number" && typeof parameterTwo !== "number" ){
            for(let i = 0; i < parameterTwo.length; i++){
                result[i] = parameterTwo[i] - parameterOne;
            }
        }
        return result;
    }
};
module.exports.multiply = function (parameterOne,parameterTwo){
    if(typeof parameterOne !== 'number' && typeof parameterTwo !== 'number' && parameterOne.length !== parameterTwo.length) {
        return 'parameterOne and parameterTwo should have the same length';
    }
    else {
        let result = [];
        if(typeof parameterOne === "number" && typeof parameterTwo == "number" ){
            result = 'is unnecesary use this function, both parameters are numbers';
        }
        else if (typeof parameterOne !== "number" && typeof parameterTwo !== "number" ){
            for(let i = 0; i < parameterOne.length; i++){
                result[i] = parameterOne[i] * parameterTwo[i];
            }
        }
        else if (typeof parameterOne !== "number" && typeof parameterTwo === "number" ){
            for(let i = 0; i < parameterOne.length; i++){
                result[i] = parameterOne[i] * parameterTwo;
            }
        }
        else if (typeof parameterOne === "number" && typeof parameterTwo !== "number" ){
            for(let i = 0; i < parameterTwo.length; i++){
                result[i] = parameterTwo[i] * parameterOne;
            }
        }
        return result;
    }
};

module.exports.collect = function(array){
    let result = 0;
    for (let i = 0; i < array.length; i++){
        result = result + array[i];
    }
    return result;
};

module.exports.fillRangeUpToDown = function(min, max, stepZise){
    let pointsNumber = Math.round((max - min )/  stepZise);
    let vector = new Array(pointsNumber);
    let counter = 0;
    for(let i = pointsNumber; i >= 0; i--) {
        vector[counter] = min + stepZise * i;
        counter++
    }
    return vector;
}
module.exports.fillRangeWithoutBounds = function(min, max, points){
    let correctPoints = points - 1;
    let stepZise = (max - min) / (points) ;
    let vector = new Array(correctPoints);
    let counter = 0;
    for(let i = 1; i <= correctPoints; i++) {
        vector[counter] = min + stepZise * i;
        counter++
    }
    return vector;
};
