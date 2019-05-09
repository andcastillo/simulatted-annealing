export default function simulatedAnnealing(inputs) {
    let testInfo = { temperatureList: [], dfList: [], optimumList: [], probabilityList: [], errorList: [], currentList: [], dxList: [], candidateList: [] };// probabilityList : [], acceptValues : []
    let current = inputs.guess;
    let range = inputs.neighbour.upperBound.map((x, i) => x - inputs.neighbour.lowerBound[i]);
    let globalOptimum = inputs.goalFunction(current);
    let optimum, optimumCandidate;
    for (let iteration = 0; iteration < inputs.maxIterations; iteration++) {
      let T = temperature(inputs.maxIterations, iteration);
      let [candidate, information] = candidateGenerator(inputs.neighbour, range, current, T, inputs.quenchingFactor);
      let candidateEvaluation = inputs.goalFunction(candidate);
      let df = candidateEvaluation - globalOptimum;
      if (df < 0) {
        current = candidate;
        optimum = candidateEvaluation;
        optimumCandidate = candidate;
        globalOptimum = candidateEvaluation;
      } else if (acceptableProbability(df, T, inputs.quenchingFactor) > Math.random()) {
        current = candidate;
        optimum = candidateEvaluation;
      }
      let probability = (df < 0 ? 2 : acceptableProbability(df, T, inputs.quenchingFactor));// only for test ** the "2" values is not showed in the plot because the probability is only generate when df > 0 ***
      testInfo.temperatureList[iteration] = T;
      testInfo.currentList[iteration] = current;
      testInfo.dxList[iteration] = information.dx;
      testInfo.candidateList[iteration] = information.parameters;
      testInfo.dfList[iteration] = df;
      testInfo.probabilityList[iteration] = probability;
    }
    return [optimumCandidate, globalOptimum, testInfo];
  }
  
  function acceptableProbability(functionDelta, temperature, quenchingFactor) {
    let probability = Math.exp(-(functionDelta * quenchingFactor) / (temperature));
    return probability;
  }
  
  function candidateGenerator(neighbour, range, current, T, quienching) {
    let dx = []; 
    let infoToTest = { yList: [], dx: undefined, parameters: undefined };
    let newCandidate = [];
    for (let i = 0; i < current.length; i++) {
      let y = Math.random() * (range[i]) + neighbour.lowerBound[i];
      infoToTest.yList[i] = y;
      dx[i] = y * Math.exp(-quienching / T);
      newCandidate[i] = current[i] + dx[i];
      newCandidate[i] = (newCandidate[i] < neighbour.lowerBound[i] ? neighbour.lowerBound[i] + (neighbour.lowerBound[i] - newCandidate) : 0) + (neighbour.lowerBound[i] <= newCandidate[i] && newCandidate[i] <= neighbour.upperBound[i] ? newCandidate[i] : 0) + (neighbour.upperBound[i] < newCandidate[i] ? neighbour.upperBound[i] - (newCandidate[i] - neighbour.upperBound[i]) : 0);
    }
    infoToTest.dx = dx;
    infoToTest.parameters = newCandidate;
    return [newCandidate, infoToTest];
  }
  
  function temperature(maxIteration, iteration) {
    let T = maxIteration / (iteration * 0.1 + 1);
    return T;
  }
  