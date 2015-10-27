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

function joinVectors(vec1, vec2) {
    if (vec1.length === vec2.length) {
        var result = [];
        for (var i = 0; i < vec1.length; i++) {
            if (vec1[i].length === 1) {
                result.push({x: vec1[i][0], y: vec2[i]});
            } else if (vec1[i].length === 2) {
                var z;
                if (vec2[i] === true) {
                    z = 1;
                } else if (vec2[i] === false) {
                    z = 0;
                } else {
                    z = vec2[i];
                }

                result.push({x: vec1[i][0], y: vec1[i][1], z: z});
            }
        }
        //console.log(result);
        return result;
    } else {
        //console.log(vec1, vec2);
        return [];
    }
}

function getDomain(data, func) {
    var min = 9999999;
    var max = -9999999;
    //console.log(data);
    for (var i = 0; i < data.length; i++) {
        var d = func(data[i]);
        if (d > max)
            max = d;
        if (d < min)
            min = d;
    }
    //console.log(min, max);
    return [min, max];
}

google.load("visualization", "1");

function render3d(rawData, container)
{
    var len = rawData.length;
    var numRows = Math.sqrt(len);
    var numCols = Math.sqrt(len);

    var tooltipStrings = new Array();
    var data = new google.visualization.DataTable();

    for (var i = 0; i < numCols; i++)
    {
        data.addColumn('number', 'col' + i);
    }
    data.addRows(numRows);
    //ar d = 360 / numRows;
    var idx = 0;
    for (var i = 0; i < numRows; i++)
    {
        for (var j = 0; j < numCols; j++)
        {
            var value = rawData[idx].z;
            data.setValue(i, j, value);

            tooltipStrings[idx] = "x:" + i + ", y:" + j + " = " + value;
            idx++;
        }
    }
    //console.log(idx);
    if (idx === 0) {
        console.log(rawData);
    }
    var containerDiv = document.getElementById(container);
    while (containerDiv.firstChild) {
        containerDiv.removeChild(containerDiv.firstChild);
    }
    var surfacePlot = new greg.ross.visualisation.SurfacePlot(containerDiv);

    // Don't fill polygons in IE. It's too slow.
    var fillPly = true;

    // Define a colour gradient.
    var colour1 = {red: 0, green: 0, blue: 255};
    var colour2 = {red: 0, green: 255, blue: 255};
    var colour3 = {red: 0, green: 255, blue: 0};
    var colour4 = {red: 255, green: 255, blue: 0};
    var colour5 = {red: 255, green: 0, blue: 0};
    var colours = [colour1, colour2, colour3, colour4, colour5];

    // Axis labels.
    var xAxisHeader = "X";
    var yAxisHeader = "Y";
    var zAxisHeader = "Z";

    var options = {xPos: 0, yPos: 0, width: 500, height: 500, fillPolygons: fillPly, colourGradient: colours,
        tooltips: tooltipStrings, xTitle: xAxisHeader, yTitle: yAxisHeader, zTitle: zAxisHeader, restrictXRotation: false};

    surfacePlot.draw(data, options);
}

function renderFunction2D(data, container, color, isIdeaGraph) {
    var vis = d3.select('#' + container),
            WIDTH = 1000,
            HEIGHT = 500,
            MARGINS = {
                top: 20,
                right: 20,
                bottom: 20,
                left: 50
            };
    if (!isIdeaGraph)
        vis.selectAll("*").remove();
    var xRange = d3.scale.linear().
            range([MARGINS.left, WIDTH - MARGINS.right]).
            domain(
                    getDomain(data, function (d) {
                        return d.x;
                    }));
    var yRange = d3.scale.linear()
            .range([HEIGHT - MARGINS.top, MARGINS.bottom])
            .domain([0, 1]);
    if (!isIdeaGraph) {
        var xAxis = d3.svg.axis()
                .scale(xRange)
                .tickSize(5)
                .tickSubdivide(true);
        var yAxis = d3.svg.axis()
                .scale(yRange)
                .tickSize(5)
                .orient('left')
                .tickSubdivide(true);


        vis.append('svg:g')
                .attr('class', 'x axis')
                .attr('transform', 'translate(0,' + (HEIGHT - MARGINS.bottom) + ')')
                .call(xAxis);

        vis.append('svg:g')
                .attr('class', 'y axis')
                .attr('transform', 'translate(' + (MARGINS.left) + ',0)')
                .call(yAxis);
    }


    var lineFunc = d3.svg.line()
            .x(function (d) {
                return xRange(d.x);
            })
            .y(function (d) {
                return yRange(d.y);
            }).interpolate('linear');

    vis.append('svg:path')
            .attr('d', lineFunc(data))
            .attr('stroke', color)
            .attr('stroke-width', 1)
            .attr('fill', 'none');
}

var ExperimentConfiguration =
        function (
                inputNum,
                inputVector,
                expFunction,
                isBooleanFunction,
                errorExpectation,
                maxGeneration,
                maxHN,
                container,
                graphType
                ) {
            this.inputNum = inputNum;
            this.inputvecs = inputVector;
            this.expFunction = expFunction;
            this.isBooleanFunction = isBooleanFunction;
            this.errorExpectation = errorExpectation;
            this.maxGeneration = maxGeneration;
            //console.log(container)
            this.maxHiddenNeurons = maxHN;
            this.container = container;
            this.graphType = graphType;
            this.population = 70;
            //console.log(this);
        };
ExperimentConfiguration.prototype.inputNumber = function (n) {
    this.inputNum = n;
    return this;
};
ExperimentConfiguration.prototype.inputVectors = function (n) {
    this.inputVector = n;
    return this;
};
ExperimentConfiguration.prototype.function = function (n) {
    this.expFunction = n;
    return this;
};
ExperimentConfiguration.prototype.isBoolean = function (n) {
    this.isBooleanFunction = n;
    return this;
};
ExperimentConfiguration.prototype.stopAtError = function (n) {
    this.errorExpectation = n;
    return this;
};
ExperimentConfiguration.prototype.MaxGeneration = function (n) {
    this.maxGeneration = n;
    return this;
};
ExperimentConfiguration.prototype.maxHiddenNeurons = function (n) {
    this.maxHN = n;
    return this;
};
ExperimentConfiguration.prototype.SVGContainer = function (n) {
    this.container = n;
    return this;
};
ExperimentConfiguration.prototype.SetGraphType = function (n) {
    this.graphType = n;
    return this;
};
ExperimentConfiguration.prototype.Population = function (n) {
    this.population = n;
    return this;
};



function experimentLoop(currentGen, config, pool, correctAnswers, finishCallback) {

    var bestGenome = null;

    var done = false;
    var bestAnswers = [];

    var resultGenome;
    var minError = 99999;
    config.container.innerHTML = "Generation " + currentGen + "/" + config.maxGeneration;

    config.container.innerHTML += "<p>best fitness = " + pool.maxFitness + "</p>";
    config.container.innerHTML += "<p>Species count = " + pool.species.length + "</p>";

    for (var d = 0; d < pool.species.length; d++) {
        var species = pool.species[d];
        //for each genome in this species
        for (var e = 0; e < species.genomes.length; e++) {
            var genome = species.genomes[e];
            //console.log(species, e);
            genome.properties["NeuroNetwork"] = generateNN(genome.chromosomes["NEAT"]);
            var error = 0;

            var currentAnswers = [];
            //var mismatch = 0;
            for (var k = 0; k < config.inputvecs.length; k++) {
                //console.log("input=", config.inputvecs[k]);

                var network = genome.properties["NeuroNetwork"];
                //console.log("input vecs= [" + config.inputvecs + "]");
                var output = evaluateNeuroNetwork(genome, network, config.inputvecs[k], false);

                var answer = output[0];
                var correctAnswer = correctAnswers[k];

                if (config.isBooleanFunction) {
                    var answer = answer > 0.5 ? true : false;
                    if (answer !== correctAnswer) {
                        var ca = correctAnswer ? 1 : 0;
                        error += Math.pow(ca - output[0], 2);
                    }
                } else {
                    var deltaI = correctAnswer - answer;
                    if (deltaI * 10 > config.errorExpectation) {
                        //mismatch++;
                    }
                    //console.log(correctAnswer, answer, deltaI);
                    error += Math.pow(deltaI, 2);
                }
                currentAnswers.push(answer);
                //console.log("output=", answer, "correct=", correctAnswer);

            }
            error = Math.sqrt(error);
            var reverseSigmoid = function (x) {
                return 2 * (1 - 1 / (Math.exp(-x * 2) + 1));
            };

            genome.fitness = reverseSigmoid(error) ;
            //genome.fitness = 1 / (error + 0.01);
            if (genome.fitness > pool.maxFitness) {
                pool.maxFitness = genome.fitness;
                //console.log("gen=" + currentGen, "new max fitness = " + pool.maxFitness, "error=" + error);



            }
            if (error < minError) {
                minError = error;
                bestGenome = genome;
                bestAnswers = currentAnswers;
                var data = joinVectors(config.inputvecs, bestAnswers);

                if (config.graphType === "2d") {
                    var container = "graph";
                    renderFunction2D(data, container, "green", false);
                } else {
                    var container = "current3d";
                    render3d(data, container);
                }
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

    if (!done && currentGen < config.maxGeneration) {
        setTimeout(function () {
            nextGeneration(pool);
            experimentLoop(currentGen + 1, config, pool, correctAnswers, finishCallback);
        }, 0);
    } else {
        finishCallback(currentGen, config, correctAnswers, bestAnswers, bestGenome, minError);
        return;
    }

}

function finishExperiment(generation, config, correctAnswers, bestAnswers, bestGenome, minError) {


    console.log(bestGenome);
    //console.log("gen = " + npool.generation, npool);



    config.container.innerHTML += "<p> Result: " + bestGenome + "</p>";
    //config.container.innerHTML += "<p> output Vec  : " + vecToStr(bestAnswers) + "</p>";
    //config.container.innerHTML += "<p> Expected Vec: " + vecToStr(correctAnswers) + "</p>";
    config.container.innerHTML += "<p> Error = " + minError + "</p>";
    //config.container.innerHTML += "<p> Generation = " + generation + "</p>";

    var table = "";
    table += "<table class='resultTable'>";

    table += "<tr><td>Input</td><td>Idea Output</td><td>Actual Output</td><td>delta</td></tr>";
    for (var i = 0; i < bestAnswers.length; i++) {
        var delta = Math.abs(correctAnswers[i] - bestAnswers[i]);
        if (config.isBooleanFunction) {
            delta = correctAnswers[i] === bestAnswers[i] ? "Same" : "Different";
        }

        table +=
                "<tr><td>" + config.inputvecs[i] +
                "</td><td>" + correctAnswers[i] +
                "</td><td>" + bestAnswers[i] +
                "</td><td>" + delta +
                "</td></tr>";
    }
    config.container.innerHTML += table + "</table>";
    showInputs();



}

function runExperiment(config) {
    //console.log(config.inputNum);

    var obj = initNEATPool(config.inputNum, 1, config.population, config.maxHiddenNeurons);

    config.container.innerHTML = "";
    config.container.innerHTML += "<p> Experiment with function " + config.expFunction.name + "</p>";
    config.container.innerHTML += "<p> and with input " + vecToStr(config.inputvecs) + "</p>";


    var correctAnswers = [];
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
    var data = joinVectors(config.inputvecs, correctAnswers);

    if (config.graphType === "2d") {
        var container = "ideagraph";
        renderFunction2D(data, container, "red");
    } else {
        var container = "idea3d";
        render3d(data, container);
    }
    experimentLoop(0, config, obj.pool, correctAnswers, finishExperiment);


}

