const errorCalculation = require('./errorCalculation');

function simulatedAnnealing(
    data,
    parameterizedFunction,
    options = {}
) {
    let {
        maxIterations = 100,
        errorTolerance = 10e-3,
        minValues,
        temperature = _temperature,
        quenchingFactor = 1,
        maxValues,
        initialValues
    } = options;

    if (!data.x || !data.y) {
        throw new Error('The data parameter must have x and y elements');
    } else if (
        !Array.isArray(data.x) ||
        data.x.length < 2 ||
        !Array.isArray(data.y) ||
        data.y.length < 2
    ) {
        throw new Error(
            'The data parameter elements must be an array with more than 2 points'
        );
    } else if (data.x.length !== data.y.length) {
        throw new Error('The data parameter elements must have the same size');
    }

    var parameters =
        initialValues || new Array(parameterizedFunction.length).fill(1);
    let parLen = parameters.length;
    maxValues = maxValues || new Array(parLen).fill(Number.MAX_SAFE_INTEGER);
    minValues = minValues || new Array(parLen).fill(Number.MIN_SAFE_INTEGER);

    if (maxValues.length !== minValues.length) {
        throw new Error('minValues and maxValues must be the same size');
    }

    if (!Array.isArray(parameters)) {
        throw new Error('initialValues must be an array');
    }

    //const globalOptimum = parameterizedFunction(parameters);
    var error = errorCalculation(data, parameters, parameterizedFunction);
    let current = parameters;
    var converged = error <= errorTolerance;
    for (var iteration = 0; iteration < maxIterations && !converged; iteration++) {
        let temp = temperature(maxIterations, iteration);
        let candidate = candidateGenerator(minValues, maxValues, current, temp, quenchingFactor);
        let currentError = errorCalculation(data, candidate, parameterizedFunction);
        let df = currentError - error;
        if (df < 0) {
            current = candidate;
            error = currentError;
            parameters = candidate;
        }
        else if (acceptableProbability(df, temp, quenchingFactor) < Math.random()) {
            current = candidate;
            error = currentError;
            parameters = candidate;
        }
        converged = error <= errorTolerance;
        console.log(error);
    }

    return {
        parameterValues: parameters,
        parameterError: error,
        iterations: iteration
    };
}

function _temperature(maxIteration, iteration) {
    return maxIteration / (iteration * 0.1 + 1);
}

function acceptableProbability(functionDelta, temperature, quenchingFactor) {
    let probability = Math.exp(- (functionDelta * quenchingFactor) / (temperature));
    //console.log(probability);
    return probability
}

function candidateGenerator(lowerBound, upperBound, current, temperature, quienching) {
    let newCandidate = new Array(current.length);
    for (let i = 0; i < current.length; i++) {
        let range = upperBound[i] - lowerBound[i];
        let dx = (Math.random() - 0.5) * range *  Math.exp(- quienching / temperature);
        newCandidate[i] = current[i] + dx;
        if(newCandidate[i] < lowerBound[i]) {
            newCandidate[i] =   lowerBound[i] + Math.abs(dx / 2); //Be sure it does not overload again
        }
        if(newCandidate[i] > upperBound[i]) {
            newCandidate[i] =  upperBound[i] - Math.abs(dx / 2); //Be sure it does not overload again
        }
    }
    //console.log(newCandidate);
    return newCandidate;
    //infoToTest.dx = dx;
    //infoToTest.parameters = newCandidate;
    //return [newCandidate, infoToTest];
}

module.exports = simulatedAnnealing;