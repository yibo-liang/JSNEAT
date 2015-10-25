

var ExperimentConfiguration =
        function (
                inputNum,
                inputVector,
                expFunction,
                isBooleanFunction,
                errorExpectation,
                maxGeneration,
                container
                ) {
            this.inputNum = inputNum;
            this.inputvecs = inputVector;
            this.expFunction = expFunction;
            this.isBooleanFunction = isBooleanFunction;
            this.errorExpectation = errorExpectation;
            this.maxGeneration = maxGeneration;
            //console.log(container)
            this.container = container;
            //console.log(this);
        };

function runExperiment(config) {
    //console.log(config.inputNum);
    var obj = initNEATPool(config.inputNum, 1, 50);
    var genomes = obj.genomes;
    var npool = obj.pool;
    var done = false;
    var resultGenome;
    var gen = 0;

    var correctAnswers = [];
    var minError = 99999;
    var minVec = 99999;
    var maxVec = -99999;
    for (var k = 0; k < config.inputvecs.length; k++) {
        var i = config.inputvecs[k];
        var r;
        if (typeof i === "number") {
            r = (config.expFunction(i));
        } else {
            if (i.length === 1) {
                r = (config.expFunction(i[0]));
            } else
                r = (config.expFunction(i));
        }
        if (r < minVec) {
            minVec = r;
        }
        if (r > maxVec) {
            maxVec = r;
        }
        correctAnswers.push(r);
    }
    //normalise to range (-1,1)

    if (!config.isBooleanFunction) {

        for (var k = 0; k < config.inputvecs.length; k++) {
            correctAnswers[k] = (correctAnswers[k] - minVec) / (maxVec - minVec);
        }

    }

    var bestGenome = null;
    var bestAnswers = [];
    do {
        for (var d = 0; d < npool.species.length; d++) {
            var species = npool.species[d];
            //for each genome in this species
            for (var e = 0; e < species.genomes.length; e++) {
                var genome = species.genomes[e];
                //console.log(species, e);
                genome.properties["NeuroNetwork"] = generateNN(genome.chromesomes["NEAT"]);
                var error = 0;

                var currentAnswers = [];
                for (var k = 0; k < config.inputvecs.length; k++) {
                    //console.log("input=", config.inputvecs[k]);

                    var network = genome.properties["NeuroNetwork"];
                    //console.log("input vecs= [" + config.inputvecs + "]");
                    var output = evaluateNeuroNetwork(genome, network, config.inputvecs[k], false);

                    var answer = output[0];
                    //if (answer === 0 && config.inputvecs[k][0] === -10) {
                    //console.log(config.inputvecs[k], answer);

                    // }



                    var correctAnswer = correctAnswers[k];

                    if (config.isBooleanFunction) {
                        var answer = answer > 0.5 ? true : false;
                        if (answer !== correctAnswer) {
                            error += 1;
                        }
                    } else {
                        var deltaI = correctAnswer - answer;
                        //console.log(correctAnswer, answer, deltaI);
                        error += Math.pow(deltaI, 2);
                    }
                    currentAnswers.push(answer);
                    //console.log("output=", answer, "correct=", correctAnswer);

                }

                //console.log(error);
                //error = Math.sqrt(error);
                var reverseSigmoid = function (x) {
                    return 2 * (1 - 1 / (Math.exp(-x * 0.5) + 1));
                };

                genome.fitness = reverseSigmoid(error);

                if (genome.fitness > npool.maxFitness) {
                    npool.maxFitness = genome.fitness;
                    console.log("gen=" + gen, "new max fitness = " + npool.maxFitness, "error=" + error);
                }

                if (error < minError) {
                    minError = error;
                    bestGenome = genome;
                    bestAnswers = currentAnswers;
                    /*
                     for (var i in network.neurons) {
                     if (network.neurons.hasOwnProperty(i)) {
                     var n = network.neurons[i];
                     console.log(i, n.cValue, n.incomingLinks, n.outcomingLinks);
                     }
                     }
                     console.log(network);
                     console.log("------------------");
                     if (bestAnswers[0] === 0 && config.inputvecs[0][0] === -10) {
                     console.log(genome);
                     console.log("ba=" + bestAnswers, "iv=" + config.inputvecs);
                     var network = genome.properties["NeuroNetwork"];
                     
                     
                     console.log("temp=", temp);
                     
                     
                     }
                    
                    var tempnet = generateNN(genome.chromesomes["NEAT"]);
                    var temp = evaluateNeuroNetwork(genome, tempnet, [0, 0], true);
                    var temp = evaluateNeuroNetwork(genome, tempnet, [0, 1], true);
                    var temp = evaluateNeuroNetwork(genome, tempnet, [1, 0], true);
                    var temp = evaluateNeuroNetwork(genome, tempnet, [1, 1], true); */
                    
                }
                if (error <= config.errorExpectation) {
                    done = true;
                    resultGenome = genome;
                    break;
                }

            }
            if (done === true) {
                break;
            }

        }
        nextGeneration(npool);
        gen++;
        //console.log(bestGenome.fitness);
        //console.log("next gen -------------------------------");
    } while (!done && gen < config.maxGeneration);
    console.log(bestGenome);
    //console.log("gen = " + npool.generation, npool);

    var vecToStr = function (vec) {
        var result = "";
        for (var i = 0; i < vec.length; i++) {
            result += "(" + vec[i] + ")";
            if (i < vec.length - 1) {
                result += ",";
            }
        }
        return result;
    };

    config.container.innerHTML += "<p> Experiment with function " + config.expFunction.name + "</p>";
    config.container.innerHTML += "<p> and with input " + vecToStr(config.inputvecs) + "</p>";
    config.container.innerHTML += "<p> Result: " + resultGenome + "</p>";
    config.container.innerHTML += "<p> output Vec  : " + vecToStr(bestAnswers) + "</p>";
    config.container.innerHTML += "<p> Expected Vec: " + vecToStr(correctAnswers) + "</p>";

    config.container.innerHTML += "<p> Min Error = " + minError + "</p>";
    config.container.innerHTML += "<p> Generation = " + gen + "</p>";





}

