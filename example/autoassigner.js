const autoassigner = require('nmr-auto-assignment');

async function autoassignSample(sample, predictor, options) {
    let OCLE = options.OCLE;
    let spectrum = sample.spectra.nmr[0].spectrum;
    let nH = sample.general.ocl.nH;
    let ranges = spectrum.getRanges({
        nH: nH,
        realTop: false,
        thresholdFactor: 0.8,
        keepPeaks: true,
        optimize: true,
        integralType: "sum",
        compile: true
    });


    let sum = 0;
    for (let j = ranges.length - 1; j >= 0; j--) {
        if (ranges[j].from < 0 || ranges[j].from > 11.8) {
            ranges.splice(j, 1);
        } else {
            if (ranges[j].from > 2.48 && ranges[j].to < 2.55) { // && signals[j].signal[0].multiplicity === 'quint') {
                ranges.splice(j, 1);
            } else
                if (ranges[j].from > 7.10 && ranges[j].to < 7.30 && ranges[j].signal[0].multiplicity === 's') {
                    ranges.splice(j, 1);
                } else {
                    sum += ranges[j].integral;
                }
        }
    }

    // Restore the integral to nH
    for (let j = ranges.length - 1; j >= 0; j--) {
        ranges[j].integral *= nH / sum;
    }

    ranges.forEach((range, index) => {
        range.signalID = "1H_" + index;
    });

    //predictor.setDb(db, 'proton', 'proton');

    var result = await autoassigner(
        {
            general: { molfile: sample.general.ocl.molfile + "" },
            spectra: {
                nmr: [{
                    nucleus: 'H',
                    experiment: '1d',
                    range: ranges,
                    solvent: sample.spectra.nmr[0].solvent
                }]
            }
        },
        {
            minScore: 1,
            maxSolutions: 2000,
            errorCS: -1,
            predictor: predictor,
            condensed: false,
            unassigned: 0,
            ignoreLabile: false,
            levels: [5, 4, 3, 2],
            OCLE: OCLE
        }
    );

    ranges.sort((a, b) => b.from - a.from);


    return { ranges, assigner: result };
}

module.exports = autoassignSample;