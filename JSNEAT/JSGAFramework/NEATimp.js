//.var MaxNodes = 30;


var useBias = true;

var isMLP = true;
var isNaiveLayered = false;
var maxNeuronPerLayer = 6;
var maxLayer = 3;
window.maxLayer = maxLayer;
var newLayerChance = 0.05;

var isRecurrent = false;


var CrossoverChance = 0.75;
var AdoptionChanche = 0.00;

var PointMutationChance = 0.3;

var LinkMutationChance = 1.5;
var NodeMutationChance = 0.6;
var BiasMutationChance = 0.20;

var DisableMutationChance = 0.4;
var EnableMutationChance = 0.2;

var MLPNewNodeChance = 0.02;
var MLPNewLayerChance = 0.02;
var MLPActivationMutateChance = 0.02;

var PerturbStep = 0.15;
var PerturbChance = 0.85;

var PruningLinkChance = 0.003;
var PruningNodeChance = 0.005;


var NeuronGenotype = function (type, index) {
    this.type = type;
    this.index = index;
};
NeuronGenotype.prototype.copy = function () {
    var result = new NeuronGenotype(this.type, this.index);
    return result;
};

function newInnovation(pool) {
    var i = pool.innovation;
    pool.innovation++;
    return i;
}

var Gene = function () {
    this.in = null;
    this.out = null;
    this.weight = 0;
    this.enable = true;
    this.innovation = 0;
    this.isOutNeuronLayerHolder = false;
    this.isMutatingActivation = false;
};

Gene.prototype.copy = function () {
    var copy = new Gene();
    copy.in = this.in;
    copy.out = this.out;
    copy.weight = this.weight;
    copy.enable = this.enable;
    copy.innovation = this.innovation;
    copy.isOutNeuronLayerHolder = this.isOutNeuronLayerHolder;
    copy.isMutatingActivation = this.isMutatingActivation;
    return copy;
};

function identity(x) {
    return x;
}

function sigmoid(x) {
    return 1 / (1 + Math.exp(-1 * x));
}

function NullActivation(x) {
    return 0;
}

function tanh(x) {
    return (Math.exp(2 * x) - 1) / (Math.exp(2 * x) + 1);
}

function cosine(x) {
    return Math.cos(x);
}

function gaussian(x) {
    return Math.exp(-(x * x) / 2);
}


function randomActivation() {
    var activationFuncs = [sigmoid, NullActivation, tanh, cosine, gaussian];
    var i = RandomIntInclusive(0, activationFuncs.length - 1);
    return activationFuncs[i];
}

var Neuron = function () {
    this.incomingLinks = [];
    this.outcomingLinks = [];
    this.cValue = 0;
    this.nValue = 0;
    this.activationFunction = randomActivation();
};
Neuron.prototype.copy = function () {
    var copy = new Neuron();
    for (var i = 0; i < this.incomingLinks.length; i++) {
        copy.incomingLinks.push(this.incomingLinks[i].copy());
    }
    for (var i = 0; i < this.outcomingLinks.length; i++) {
        copy.outcomingLinks.push(this.outcomingLinks[i].copy());
    }
    copy.cValue = this.cValue;
    copy.nValue = this.nValue;
    copy.activationFunction = this.activationFunction;
    return copy;
}


var Network = function (chromosome) {
    this.chromosome = chromosome;
    this.neurons = {};
    this.indices = [];
};
Network.prototype.copy = function () {
    var copy = new Network(null);
    for (var i in this.neurons) {
        if (this.neurons.hasOwnProperty(i)) {
            copy.neurons[i] = this.neurons[i].copy();
            copy.indices.push(i);
        }
    }
    return copy;

};



function NeuronCompare(n1, n2) {
    if (!isNaiveLayered && !isMLP) {
        var getType = function (i) {
            return i.substring(0, 1);
        };
        var getNumber = function (i) {
            return +i.substring(1);
        };

        var typeToNum = function (t) {
            if (t === "I") {
                return 0;
            }
            if (t === "B") {
                return 10000;
            }
            if (t === "H") {
                return 20000;
            }
            if (t === "O") {
                return 30000;
            }
        };
        var typea = getType(n1);
        var typeb = getType(n2);
        var na = getNumber(n1) + typeToNum(typea);
        var nb = getNumber(n2) + typeToNum(typeb);
        return na - nb;
    } else {
        var n1vec = n1.split(".");
        var n2vec = n2.split(".");
        var getType = function (vec) {
            return vec[0];
        };
        var getNumber = function (vec) {
            if (vec[0] === "B")
                return 0;
            return +vec[vec.length - 1];
        };
        var getLayer = function (vec) {
            var t = getType(vec);
            if (t === "B" || t === "I") {
                return 0;
            }
            if (t === "H") {
                return parseInt(vec[1]);
            }
            if (t === "O") {
                return 0;
            }
        };

        var typeToNum = function (t) {
            if (t === "I") {
                return 0;
            }
            if (t === "B") {
                return 1000;
            }
            if (t === "H") {
                return 20000;
            }
            if (t === "O") {
                return 30000;
            }
        };

        var typea = getType(n1vec);
        var typeb = getType(n2vec);
        var na = getNumber(n1vec) + typeToNum(typea) + getLayer(n1vec) * 100;
        var nb = getNumber(n2vec) + typeToNum(typeb) + getLayer(n2vec) * 100;
        //console.log(n1, n2);
        //console.log("na", na, "nb", nb)
        //console.log(n1vec, getNumber(n1vec), typeToNum(typea), getLayer(n1vec));
        //console.log(n2vec, getNumber(n2vec), typeToNum(typeb), getLayer(n2vec));

        return na - nb;

    }

}

function isSameLayer(n1, n2) {
    if (!isNaiveLayered) {
        var a = n1.substring(0, 1);
        var b = n2.substring(0, 1);
        return a === b;

    } else {
        var n1vec = n1.split(".");
        var n2vec = n2.split(".");
        var getType = function (vec) {
            return vec[0];
        };
        if (n1vec[0] === n2vec[0]) {
            if (n1vec[0] === "H") {
                var getLayer = function (vec) {
                    var t = getType(vec);
                    if (t === "B" || t === "I") {
                        return 0;
                    }
                    if (t === "H") {
                        return parseInt(vec[1]);
                    }
                    if (t === "O") {
                        return 0;
                    }
                };
                var l1 = getLayer(n1vec);
                var l2 = getLayer(n2vec);
                return l1 === l2;
            } else {
                return true;
            }
        } else {
            return false;
        }

    }
}

function generateNN(chromosome) {
    var network = new Network(chromosome);

    //console.log(chromosome);
    var inputNeuronNum = chromosome.genome.properties["inputs"];
    var outputNeuronNum = chromosome.genome.properties["outputs"];
    for (var i = 0; i < inputNeuronNum; i++) {
        var n = new Neuron();
        var index = "I" + i;
        if (isNaiveLayered || isMLP) {
            index = "I." + i;
        }
        network.neurons[index] = n;
        network.indices.push(index);
    }
    for (var i = 0; i < outputNeuronNum; i++) {
        var n = new Neuron();
        var index = "O" + i;
        if (isNaiveLayered || isMLP) {
            index = "O." + i;
        }
        network.neurons[index] = n;
        network.indices.push(index);
    }
    if (useBias === true) {
        var n = new Neuron();
        n.cValue = 1;
        n.nValue = 1;
        network.neurons["B"] = n;
        network.indices.push("B");
    }


    for (var i = 0; i < chromosome.genes.length; i++) {
        var gene = chromosome.genes[i];
        if (gene.enable === true) {
            if (typeof network.neurons[gene.out] === "undefined") {
                network.neurons[gene.out] = new Neuron();
                network.indices.push(gene.out);

            }
            if (gene.isOutNeuronLayerHolder === true) {
                network.neurons[gene.out].activationFunction = randomActivation();//identity;
            } else if (gene.isMutatingActivation) {
                network.neurons[gene.out].activationFunction = randomActivation();
                gene.isMutatingActivation = false;
            }
            var neuron_i = gene.out;
            var neuron = network.neurons[neuron_i];

            neuron.incomingLinks.push(gene);

        }
    }


    for (var i = 0; i < chromosome.genes.length; i++) {
        var gene = chromosome.genes[i];
        if (gene.enable === true) {
            if (typeof network.neurons[gene.in] === "undefined") {
                network.neurons[gene.in] = new Neuron();
                network.indices.push(gene.in);
            }
            var neuron_i = gene.in;
            var neuron = network.neurons[neuron_i];

            neuron.outcomingLinks.push(gene);

        }
    }

    chromosome.genome.properties["NeuroNetwork"] = network;
    //console.log("generate nn=", network);
    return network;


}

function evaluateNeuroNetwork(genome, network, inputArray, isDebuging) {
    var inputVec = [];
    for (var i = 0; i < inputArray.length; i++) {
        inputVec.push(inputArray[i]);
    }
    var inputs = genome.properties["inputs"];
    if (inputVec.length !== inputs) {
        console.log("ERROR: INPUT NUMBER IS DIFFERENT FROM NETWORK", inputVec);
        return;
    }


    for (var i = 0; i < inputs; i++) {
        if (!isNaiveLayered && !isMLP) {
            network.neurons["I" + i].cValue = inputVec[i];
        } else {
            network.neurons["I." + i].cValue = inputVec[i];
        }
    }

    //sort all neuron, from input -> bias -> hidden -> output
    var neuronList = [];
    for (var i in network.neurons) {
        if (network.neurons.hasOwnProperty(i)) {
            neuronList.push(i);
        }
    }
    neuronList.sort(NeuronCompare);
    //if (neuronList.length >= 5) {
    //    console.log(neuronList);
    //    die();
    //}

    var debugstr = "";
    for (var i in neuronList) {
        var ni = neuronList[i];
        if (ni === "B" || ni.substring(0, 1) === "I")
            continue
        var neuron = network.neurons[ni];

        debugstr += "Debug neuron ni='" + ni + "';\n";

        //console.log(i, neuron);
        var sum = 0;
        for (var j = 0; j < neuron.incomingLinks.length; j++) {
            var incoming = neuron.incomingLinks[j];
            var other = network.neurons[incoming.in];
            //console.log(incoming, other);
            sum += incoming.weight * other.cValue;
            debugstr += "value += n[" + incoming.in + "](" + other.cValue + ")*" + incoming.weight + "\n";
            //console.log("sum=",sum);
        }
        debugstr += "value=" + sum + ", " + neuron.activationFunction.name + "(value)=" + neuron.activationFunction(sum) + "\n";
        if (neuron.incomingLinks.length > 0) {

            var temp = neuron.activationFunction(sum);
            //if recurrent network is used, then save the new value for next iteration,
            //otherwise, save the new value for current iteration
            if (isRecurrent === true) {
                neuron.nValue = temp;
            } else {
                neuron.cValue = temp;
            }
            //console.log("v=",neuron.v);
        }

        //console.log("i=", i, "  v=", network.neurons[i].v);

    }

    if (typeof isDebuging !== "undefined" && isDebuging) {

        console.log(debugstr);
        //console.log(network);

    }
    if (isRecurrent === true)
        for (var i in network.neurons) {
            if (network.neurons.hasOwnProperty(i)) {
                var neuron = network.neurons[i];
                neuron.cValue = neuron.nValue;
            }
        }

    var outputArray = [];

    var outputs = genome.properties["outputs"];
    for (var i = 0; i < outputs; i++) {
        var v;
        if (!isNaiveLayered && !isMLP) {
            v = network.neurons["O" + i].cValue;
        } else {
            v = network.neurons["O." + i].cValue;
        }
        //if (v > 0) {
        //    outputArray.push(true);
        //} else {
        //    outputArray.push(false);
        //}
        outputArray.push(v);
    }
    if (isDebuging) {
        console.log("input=" + inputArray);
        console.log("output=" + outputArray);
        console.log("---------------------");
    }
    return outputArray;

}

function disjoint(geneslist1, geneslist2) {
    var l1 = {};
    var l2 = {};

    for (var i = 0; i < geneslist1.length; i++) {
        var gene = geneslist1[i];
        l1[gene.innovation] = true;
    }

    for (var i = 0; i < geneslist2.length; i++) {
        var gene = geneslist2[i];
        l2[gene.innovation] = true;
    }
    //console.log(l1, l2);
    var disjointGenes = 0;
    for (var i = 0; i < geneslist1.length; i++) {

        var gene = geneslist1[i];
        //console.log(i, gene, genes1);
        if (typeof l2[gene.innovation] === "undefined" || l2[gene.innovation] !== true) {
            disjointGenes = disjointGenes + 1;
        }

    }
    for (var i = 0; i < geneslist2.length; i++) {

        var gene = geneslist2[i];
        if (typeof l1[gene.innovation] === "undefined" || l1[gene.innovation] !== true) {
            disjointGenes = disjointGenes + 1;
        }

    }

    var n = Math.max(geneslist1.length, geneslist2.length);
    //console.log(disjointGenes,n);
    if (n === 0)
        return 0;
    return disjointGenes / n;
}


function weights(geneslist1, geneslist2) {
    //console.log("weights:",genes1, genes2);
    var i2 = {};
    for (var i = 1; i < geneslist2.length; i++) {
        var gene = geneslist2[i];
        i2[+gene.innovation] = gene;
    }
    //console.log(i2);
    var sum = 0;
    var coincident = 0;
    for (var i = 1; i < geneslist1.length; i++) {
        var gene = geneslist1[i];
        if (typeof i2[+gene.innovation] !== "undefined") {
            var gene2 = i2[+gene.innovation];
            sum = sum + Math.abs(gene.weight - gene2.weight);
            coincident = coincident + 1;
        }
    }
    //console.log("sum, coinci", sum, coincident, sum/coincident);
    if (coincident === 0) {
        return 1;
    }
    return sum / coincident;
}

var NEATcrossover = function (chromosome1, chromosome2) {


    var childchromosome = new Chromesome(
            chromosome1.genome,
            "NEAT",
            chromosome1.evaluator,
            chromosome1.chromosomeOperator,
            chromosome1.deltaFunction,
            chromosome1.deltaWeight
            );
    var innovations2 = {};
    for (var i = 0; i < chromosome2.genes.length; i++) {
        var gene = chromosome2.genes[i];
        innovations2[gene.innovation] = gene;
    }


    for (var i = 0; i < chromosome1.genes.length; i++) {
        var gene1 = chromosome1.genes[i];
        var gene2 = innovations2[gene1.innovation];
        if (typeof gene2 !== "undefined" && gene2.enable === true) {
            childchromosome.genes.push(gene2.copy());
        } else {
            childchromosome.genes.push(gene1.copy());
        }
    }

    return childchromosome;
};

function isInputNeuron(index) {
    if (index.substring(0, 1) === "I") {
        return true;
    }
    return false;
}

function isBiasNeuron(index) {
    return index === "B";
}

function randomNeuron(genome, NoInputs, NoBias) {
    var candidates = [];

    //console.log(genome.properties["NeuroNetwork"]);
    var indices = genome.properties["NeuroNetwork"].indices;
    for (var i = 0; i < indices.length; i++) {
        var index = indices[i];
        if (NoInputs && isInputNeuron(index) || NoBias && isBiasNeuron(index)) {

        } else {
            candidates.push(index);
        }
    }
    var i = RandomIntInclusive(0, candidates.length - 1);
    //console.log(candidates)
    return candidates[i];
}

function containsLink(genes, link) {
    for (var i = 0; i < genes.length; i++) {
        var gene = genes[i];
        if (gene.in === link.in && gene.out === link.out) {
            return true;
        }
    }
    return false;
}


function pointMutate(chromosome) {
    //console.log("pointMutation");
    var step = chromosome.genome.properties["PertubStep"];
    var chance = chromosome.genome.properties["PertubChance"];
    for (var i = 0; i < chromosome.genes.length; i++) {
        var gene = chromosome.genes[i];
        if (Math.random() < chance) {
            gene.weight = gene.weight + Math.random() * step * 2 - step;

        } else {
            gene.weight = Math.random() * step * 8 - 4;
        }
    }

}



function linkMutation(chromosome, forceBias) {
    //console.log("linkMutation");

    var neuron1 = randomNeuron(chromosome.genome, false, false);
    var neuron2 = randomNeuron(chromosome.genome, true, true);
    var inputs = chromosome.genome.properties["inputs"];
    if (!isNaiveLayered) {
        if (isInputNeuron(neuron1) && isInputNeuron(neuron2)) {
            // both inputs neuron
            return;
        }

        if (neuron1 === neuron2 && !isRecurrent) {
            //if input is output but recurrent is not allowed
            return;
        }

        if (isInputNeuron(neuron2)) {
            //swap 1 and 2
            var temp = neuron1;
            neuron1 = neuron2;
            neuron2 = temp;
        }
    } else {

        var getLayer = function (n) {
            var vec = n.split(".");
            var t = vec[0];
            if (t === "B" || t === "I") {
                return 0;
            }
            if (t === "H") {
                return parseInt(vec[1]);
            }
            if (t === "O") {
                return 0;
            }
        };
        if (isSameLayer(neuron1, neuron2) && getLayer(neuron1) !== "H") {
            return;
        }
        var getType = function (i) {
            return i.substring(0, 1);
        };
        var l1 = getLayer(neuron1);
        var l2 = getLayer(neuron2);
        if (l1 < l2) {
            var tl = l2;
            var l2 = l1;
            var l1 = tl;
        }
        if (getType(neuron1) === "I" || getType(neuron1) === "B") {
            if (l2 > 1) {
                return;
            }
        }
        if (getType(neuron1) === "H" && getType(neuron2) === "O") {
            var maxl = chromosome.genome.properties["hiddenLayers"];
            if (l1 < maxl) {
                return;
            }
        }

    }

    var newLink = new Gene();
    newLink.in = neuron1;
    newLink.out = neuron2;
    if (forceBias) {
        newLink.in = "B";
    }
    if (containsLink(chromosome.genes, newLink)) {
        return;
    }
    newLink.innovation = newInnovation(chromosome.genome.pool);
    //console.log("new link=", newLink)
    newLink.weight = Math.random() * 10 - 5;
    chromosome.genes.push(newLink);
}


function nodeMutation(chromosome) {
    //console.log("nodeMutation");
    var genome = chromosome.genome;
    if (genome.properties["hiddenNeurons"] >= genome.properties["maxHiddenNeurons"]) {
        return;
    }


    if (chromosome.genes.length === 0) {
        return;
    }
    if (!isNaiveLayered) {
        var gene = chromosome.genes[RandomIntInclusive(0, chromosome.genes.length - 1)];
        if (gene.enable === false) {
            return;
        }
        genome.properties["hiddenNeurons"] += 1;

        gene.enable = false;
        var gene1 = gene.copy();
        gene1.out = "H" + genome.properties["hiddenNeurons"];
        gene1.weight = 1.0;
        gene1.innovation = newInnovation(genome.pool);
        gene1.enable = true;
        chromosome.genes.push(gene1);
        var gene2 = gene.copy();
        gene2.in = "H" + genome.properties["hiddenNeurons"];
        gene2.innovation = newInnovation(genome.pool);
        gene2.enable = true;
        chromosome.genes.push(gene2);
        //console.log(gene, ",", gene1, ",", gene2);
    } else {
        var layer;
        if (genome.properties["hiddenLayers"] === 0) {
            genome.properties["hiddenLayers"]++;
            genome.properties["NeuronNumberAtLayer"][genome.properties["hiddenLayers"]] = 0;
            layer = 1;
        } else {
            layer = RandomIntInclusive(1, genome.properties["hiddenLayers"]);
            //console.log('genome.properties["hiddenLayers"]=' + genome.properties["hiddenLayers"], "layer=" + layer);
            //console.log('genome.properties["hiddenNeurons"]=' + genome.properties["hiddenNeurons"])
            if (genome.properties["hiddenNeurons"] > 1) {
                var p = newLayerChance;
                if (Math.random() < p) {
                    genome.properties["hiddenLayers"]++;
                    layer = genome.properties["hiddenLayers"];
                    genome.properties["NeuronNumberAtLayer"][layer] = 0;
                }
            }
        }
        var candidates = [];
        var getType = function (i) {
            return i.substring(0, 1);
        };
        var getLayer = function (n) {
            //console.log(n);
            var vec = n.split(".");
            var t = vec[0];
            if (t === "B" || t === "I") {
                return 0;
            }
            if (t === "H") {
                return parseInt(vec[1]);
            }
            if (t === "O") {
                return 0;
            }
        };
        //console.log(chromosome.genes);
        for (var i = 0; i < chromosome.genes.length; i++) {
            //console.log("----")
            //console.log(chromosome.genes[i], "in.layer=" + getLayer(chromosome.genes[i].in), layer - 1);

            if (chromosome.genes[i].enable === false) {
                //console.log("ignore=", chromosome.genes[i], layer)
                continue;
            }
            if (getLayer(chromosome.genes[i].in) === layer - 1) {
                //console.log("found=" + chromosome.genes[i]);
                candidates.push(chromosome.genes[i]);

            } else if (layer === 1 && (getType(chromosome.genes[i].in) === "I" || getType(chromosome.genes[i].in) === "B")) {
                candidates.push(chromosome.genes[i]);
            } else {
                //console.log("ignore=", chromosome.genes[i], layer)
            }
        }
        //console.log("***")
        //console.log(layer, chromosome.genes, chromosome, candidates);
        var index = RandomIntInclusive(0, candidates.length - 1);
        var gene = candidates[index];
        if (typeof gene === "undefined") {
            //console.log(genome, candidates, index);
            //console.log(genome.log);
            die();
        }

        //console.log(genome.properties["hiddenNeurons"], "increment hn", ", gene len=" + chromosome.genes.length);

        genome.properties["hiddenNeurons"] += 1;
        genome.log += 'genome.properties["hiddenNeurons"]=' + genome.properties["hiddenNeurons"] + "\n";
        genome.log = genome.log + "chromosome id=" + chromosome.id + "\n";
        genome.log = genome.log + "add one neuron\n";
        gene.enable = false;
        var gene1 = gene.copy();
        gene1.out = "H." + layer + "." + genome.properties["hiddenNeurons"];
        gene1.weight = 1;
        gene1.innovation = newInnovation(genome.pool);
        gene1.enable = true;
        genome.log = genome.log + "add one gene: in=" + gene1.in + ", out=" + gene1.out + ", gene len=" + chromosome.genes.length + "\n";
        chromosome.genes.push(gene1);
        var gene2 = gene.copy();
        gene2.in = "H." + layer + "." + genome.properties["hiddenNeurons"];
        gene2.weight = 1;
        gene2.innovation = newInnovation(genome.pool);
        gene2.enable = true;

        genome.log = genome.log + "add one gene: in=" + gene2.in + ", out=" + gene2.out + ", gene len=" + chromosome.genes.length + "\n";
        chromosome.genes.push(gene2);

        //console.log(genome.properties["hiddenNeurons"], "increment hn after", "len=" + chromosome.genes.length);

    }

}



function enableDisableMutate(chromosome, enable) {
    if (isNaiveLayered) {
        //does not allow enable in MLP network
        return;

    }

    var candidates = [];
    for (var i = 0; i < chromosome.genes.length; i++) {
        if (chromosome.genes[i].enable !== enable) {
            candidates.push(chromosome.genes[i]);
        }
    }
    if (candidates.length === 0) {
        return;
    }
    var gene = candidates[RandomIntInclusive(0, candidates.length - 1)];
    gene.enable = enable;
}


function deleteOneLink(chromosome, gene) {
    var len1 = chromosome.genes.length;
    for (var i = 0; i < chromosome.genes.length; i++) {
        var g = chromosome.genes[i];
        if (gene.in === g.in && g.out === gene.out)
        {

            var deleted = chromosome.genes.splice(i, 1);
            var len2 = chromosome.genes.length;
            //console.log("delete link", gene.in, gene.out, "lens=", len1, len2);
            return;
        }
    }
}

function deleteDisabledLinkMutate(chromosome) {

    var candidates = [];
    for (var i = 0; i < chromosome.genes.length; i++) {
        var g = chromosome.genes[i];
        if (g.enable === false)
            candidates.push(g);

    }

    if (candidates.length === 0)
        return;

    var i = getRandomIntInclusive(0, candidates.length - 1);
    console.log("link found..");
    deleteOneLink(chromosome, candidates[i]);
    generateNN(chromosome);


}

function MLPActivationMutate(chromosome) {
    var l = chromosome.genes.length - 1;
    var i = RandomIntInclusive(0, l);
    chromosome.genes[i].isMutatingActivation = true;
}


function MLPNewLayerMutate(chromosome) {
    if (!isMLP) {
        return;
    }
    var genome = chromosome.genome;
    var oldMaxHiddenLayerNum = genome.properties["hiddenLayers"];
    var newLayerNum = oldMaxHiddenLayerNum + 1;

    if (newLayerNum > maxLayer)
        return;

    if (genome.properties["hiddenNeurons"] + 1 > maxLayer * maxNeuronPerLayer)
        return;

    var outputs = genome.properties["outputs"];
    var newNeuronNum = outputs;


    /* debug 
     if (newLayerNum > 1) {
     console.log("----------------------------------------------")
     console.log("before mutate");
     for (var i = 0; i < chromosome.genes.length; i++) {
     var gene = chromosome.genes[i];
     //console.log(gene.in, "->", gene.out, "w=" + gene.weight);
     }
     
     var tm = generateNN(chromosome);
     var result = evaluateNeuroNetwork(chromosome.genome, tm, [1], true);
     console.log("result=", result)
     }
     --*/
    //generate new neurons' names
    var newNeuronNames = [];
    for (var i = 0; i < newNeuronNum; i++) {
        newNeuronNames.push("H." + newLayerNum + "." + (i + 1));
    }

    //reconnect all link to output layer to new Layer, 
    for (var n = 0; n < outputs; n++) {
        var outputNeuronName = "O." + (n);
        for (var i = 0; i < chromosome.genes.length; i++) {
            var gene = chromosome.genes[i];
            if (gene.out === outputNeuronName) {
                gene.out = newNeuronNames[n];
                if (newLayerNum > 1)
                    gene.isOutNeuronLayerHolder = true;
            }
            //console.log("reconnected gene", gene)
        }
    }

    if (useBias) {
        for (var i = 0; i < newNeuronNum; i++) {
            var biaslink = new Gene();
            biaslink.in = "B";
            biaslink.out = newNeuronNames[i];
            biaslink.innovation = newInnovation(genome.pool);
            biaslink.weight = Math.random() * 10 - 5;
            chromosome.genes.push(biaslink);
        }
    }

    //and connect new Layer to ouput
    for (var n = 0; n < outputs; n++) {
        for (var i = 0; i < newNeuronNum; i++) {
            var newGene = new Gene();
            newGene.in = newNeuronNames[i];
            newGene.out = "O." + n;
            newGene.innovation = newInnovation(genome.pool);
            if (n === i) {
                newGene.weight = 1;
            } else {
                newGene.weight = 0;
            }
            newGene.enable = true;
            chromosome.genes.push(newGene);
        }
    }


    genome.properties["hiddenLayers"] += 1;
    genome.properties["hiddenNeurons"] += outputs;
    genome.properties["NeuronNumberAtLayer"][newLayerNum] = outputs;
    /*
     
     if (newLayerNum > 1) {
     console.log("after mutate");
     for (var i = 0; i < chromosome.genes.length; i++) {
     var gene = chromosome.genes[i];
     //console.log(gene.in, "->", gene.out, "w=" + gene.weight);
     }
     
     var tm = generateNN(chromosome);
     var result = evaluateNeuroNetwork(chromosome.genome, tm, [1], true);
     console.log("result=", result)
     // die();
     }*/
}

function MLPNewNodeMutate(chromosome) {

    if (!isMLP) {
        return;
    }
    var genome = chromosome.genome;
    var oldMaxHiddenLayerNum = genome.properties["hiddenLayers"];
    var inputs = genome.properties["inputs"];
    var outputs = genome.properties["outputs"];



    var selectedLayer = Math.floor(Math.random() * oldMaxHiddenLayerNum) + 1;
    var newNeuronIndex = genome.properties["NeuronNumberAtLayer"][selectedLayer] + 1;
    if (newNeuronIndex > maxNeuronPerLayer) {
        return;
    }

    if (genome.properties["hiddenNeurons"] + 1 > maxLayer * maxNeuronPerLayer)
        return;

    var newNodeNmae = "H." + selectedLayer + "." + newNeuronIndex;


    if (oldMaxHiddenLayerNum === 0) {
        //does not have a hidden layer, mutate exit
        return;
    }


    var newLinksIn = [];
    if (selectedLayer === 1) {
        for (var i = 0; i < inputs; i++) {
            newLinksIn.push("I." + i);
        }

    } else {
        for (var i = 0; i < genome.properties["NeuronNumberAtLayer"][selectedLayer - 1]; i++) {
            newLinksIn.push("H." + (selectedLayer - 1) + "." + (i + 1));
        }
    }
    if (useBias) {
        newLinksIn.push("B");
    }


    var newLinksOut = [];
    if (selectedLayer === oldMaxHiddenLayerNum) {
        for (var i = 0; i < outputs; i++) {
            newLinksOut.push("O." + i);
        }
    } else {
        for (var i = 0; i < genome.properties["NeuronNumberAtLayer"][selectedLayer + 1]; i++) {
            newLinksOut.push("H." + (selectedLayer + 1) + "." + (i + 1));
        }
    }

    var randomWeight = function () {
        return Math.random() * 10 - 5;
    };

    //console.log("new node with name" + newNodeNmae, genome.properties);
    for (var i = 0; i < newLinksIn.length; i++) {
        var newGene = new Gene();
        newGene.weight = randomWeight();
        newGene.in = newLinksIn[i];
        newGene.out = newNodeNmae;
        newGene.innovation = newInnovation(chromosome.genome.pool);
        chromosome.genes.push(newGene);
    }


    for (var i = 0; i < newLinksOut.length; i++) {
        var newGene = new Gene();
        newGene.weight = randomWeight();
        newGene.in = newNodeNmae;
        newGene.out = newLinksOut[i];
        newGene.innovation = newInnovation(chromosome.genome.pool);
        chromosome.genes.push(newGene);
    }

    genome.properties["hiddenNeurons"] += 1;
    genome.properties["NeuronNumberAtLayer"][selectedLayer] += 1;
}


var NEATmutate = function (chromosome) {

    if (Math.random() < chromosome.genome.properties["PointMutationChance"]) {
        pointMutate(chromosome);

    }
    var genome = chromosome.genome;

    if (!isMLP) {
        var p = genome.properties["LinkMutationChance"];
        while (p > 0) {
            if (Math.random() < p) {
                linkMutation(chromosome, false);
            }
            p -= 1;
        }

        var p = genome.properties["BiasMutationChance"];
        while (p > 0) {
            if (Math.random() < p) {
                linkMutation(chromosome, true);
            }
            p -= 1;
        }

        var p = genome.properties["NodeMutationChance"];
        while (p > 0) {
            if (Math.random() < p) {
                nodeMutation(chromosome);
            }
            p -= 1;
        }

        var p = genome.properties["EnableMutationChance"];
        while (p > 0) {
            if (Math.random() < p) {
                enableDisableMutate(chromosome, true);
            }
            p -= 1;
        }
        var p = genome.properties["DisableMutationChance"];
        while (p > 0) {
            if (Math.random() < p) {
                enableDisableMutate(chromosome, false);
            }
            p -= 1;
        }

    } else {

        var p = genome.properties["MLPNewLayerChance"];
        while (p > 0) {
            if (Math.random() < p) {
                //console.log("mutate.MLPNewNodeMutate");
                MLPNewNodeMutate(chromosome);
            }
            p -= 1;
        }

        var p = genome.properties["MLPNewLayerChance"];
        while (p > 0) {
            if (Math.random() < p) {
                // console.log("mutate.MLPNewLayerMutate");

                MLPNewLayerMutate(chromosome);
            }
            p -= 1;
        }
    }
    generateNN(chromosome);
    return chromosome;


};


function DeltaNeurons(chromosome1, chromosome2) {
    var g1 = chromosome1.genome;
    var g2 = chromosome2.genome;
    var hn1 = g1.properties["NeuronNumberAtLayer"];
    var hn2 = g2.properties["NeuronNumberAtLayer"];
    //console.log(hn1,hn2);
    var sum = 0;
    for (var i in hn1) {
        if (hn1.hasOwnProperty(i)) {
            //console.log("i=", i)
            var n1 = hn1[i];
            var n2 = hn2[i];

            if (typeof n2 === "undefined" && n1 !== null) {
                sum += n1;
                //console.log("d layer +"+n1);
            } else {
                if (n1 !== null && n2 !== null)
                    sum += Math.abs(n1 - n2);
                // console.log("d neuron n +1")


            }
        }
    }
    // if (sum >= 1) {

    //return sum;
    //}
    var getLayerName = function (d) {
        return d.substring(0, 1);
    }
    var ns1 = g1.properties["NeuroNetwork"];
    var ns2 = g2.properties["NeuroNetwork"];
    for (var i in ns1.neurons) {
        if (ns1.neurons.hasOwnProperty(i)) {
            var n = getLayerName(i)
            if (n === "I" || n === "B") {
                continue;
            }

            if (typeof ns2.neurons[i] !== "undefined") {
                n1 = ns1.neurons[i];
                n2 = ns2.neurons[i];

                if (n1.activationFunction !== n2.activationFunction) {
                    sum += 0;
                    // console.log("d fun +1");
                }
            }
        }
    }

    /*
     console.log("=====================================")
     console.log(sum)
     console.log(ns1, ns2)
     
     */
    //die()
    if (sum > 1) {
        //console.log(ns1, ns2, sum)
        //die()
    }
    return sum;


}

var NEAToperator = new chromosomeOperator(NEATcrossover, NEATmutate);
var NEATdeltaWeight = 1;


function NEATdeltaFunction(chromosome1, chromosome2) {
    var DeltaDisjoint = 0;
    var DeltaWeights = 0.5;
    var DeltaNeuron = 1;

    var dd = DeltaDisjoint * disjoint(chromosome1.genes, chromosome2.genes);
    var dw = DeltaWeights * weights(chromosome1.genes, chromosome2.genes);
    var dn = DeltaNeuron * DeltaNeurons(chromosome1, chromosome2);
    //console.log(dd, dw, dn);
    return dd + dw + dn;
}

function baseNEATChromesome(genome, initInputs, initOutputs, maxHN) {
    genome.properties["inputs"] = initInputs;
    genome.properties["outputs"] = initOutputs;
    genome.properties["hiddenNeurons"] = 0;
    genome.properties["NeuronNumberAtLayer"] = [];

    genome.properties["hiddenLayers"] = 0;

    genome.properties["maxHiddenNeurons"] = maxHN;

    genome.properties["PertubStep"] = PerturbStep;
    genome.properties["PertubChance"] = PerturbChance;

    genome.properties["PointMutationChance"] = PointMutationChance;
    genome.properties["LinkMutationChance"] = LinkMutationChance;
    genome.properties["NodeMutationChance"] = NodeMutationChance;
    genome.properties["BiasMutationChance"] = BiasMutationChance;
    genome.properties["DisableMutationChance"] = DisableMutationChance;
    genome.properties["EnableMutationChance"] = EnableMutationChance;

    genome.properties["MLPNewNodeChance"] = MLPNewNodeChance;
    genome.properties["MLPNewLayerChance"] = MLPNewLayerChance;


    //var inputs = genome.properties["inputs"];
    //var outputs = genome.properties["outputs"];

    var result = new Chromesome(
            genome,
            "NEAT",
            generateNN,
            NEAToperator,
            NEATdeltaFunction,
            NEATdeltaWeight
            );

    var localInnovation = 0;
    if (isMLP)
        for (var i = 0; i <= initInputs; i++) {
            for (var j = 0; j < initOutputs; j++) {
                var gene = new Gene();
                var link_in;
                var link_out;
                if (i === initInputs) {
                    link_in = "B";
                } else {
                    link_in = "I." + i;
                }
                link_out = "O." + j;
                gene.in = link_in;
                gene.out = link_out;
                gene.innovation = localInnovation++;
                gene.weight = Math.random() * 6 - 3;
                result.genes.push(gene);
            }
        }

    return result;
}
;


function totalComplexity(pool) {
    //console.log(pool)
    var sum = 0;
    var population = 0;
    for (var i = 0; i < pool.species.length; i++) {
        var species = pool.species[i];
        for (var j in species.genomes) {
            if (species.genomes.hasOwnProperty(j)) {
                population++;

                if (typeof species.genomes[j].properties["NeuroNetwork"] === "undefined")
                    generateNN(species.genomes[j].chromosomes["NEAT"]);

                var network = species.genomes[j].properties["NeuroNetwork"];
                //console.log(network);
                for (var ni in network.neurons) {
                    if (network.neurons.hasOwnProperty(ni)) {
                        sum++;
                        for (var l = 0; l < network.neurons[ni].incomingLinks.length; l++) {
                            sum++;
                        }
                    }
                }
            }
        }
    }

    var result = sum / population;
    //console.log("complexity", sum, population, result);
    return result;
}

function basicNEATGenome(pool, inputn, outputn) {
    var basicGenome = new Genome(pool);
    var chromosome = baseNEATChromesome(basicGenome, inputn, outputn);
    basicGenome.chromosomes["NEAT"] = chromosome;
    basicGenome.properties["NeuroNetwork"] = generateNN(chromosome);
    //console.log(chromosome);
    chromosome.chromosomeOperator.mutateOperator(chromosome);
    return basicGenome;
}

function initNEATPool(inputn, outputn, population, maxHN) {
    var init_genomes = [];
    var pool = new Pool(population, true);

    for (var i = 0; i < population; i++) {
        var basicGenome = new Genome(pool);
        var chromosome = baseNEATChromesome(basicGenome, inputn, outputn, maxHN);
        basicGenome.chromosomes["NEAT"] = chromosome;
        basicGenome.properties["NeuroNetwork"] = generateNN(chromosome);
        pool.crossoverEnabled = !isMLP;
        pool.mutateEnabled = true;


        chromosome.chromosomeOperator.mutateOperator(chromosome);
        addToSpecies(basicGenome, pool);
        init_genomes.push(basicGenome);
    }

    pool.complexity = totalComplexity(pool);
    return {pool: pool, genomes: init_genomes};
}